## ============================================================
## Week 4 Day 2 — NAFLD Microarray Validation: limma + GSEA
## Dataset: GSE89632 (human liver, Illumina HumanHT-12 BeadChip)
## Groups: HC (n=24) vs SS (n=20) vs NASH (n=19)
## Purpose: Validate Day 1 DESeq2 findings using independent platform
## ============================================================

suppressPackageStartupMessages({
  library(GEOquery)
  library(limma)
  library(ggplot2)
  library(ggrepel)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(msigdbr)
  library(dplyr)
  library(tidyr)
})

set.seed(42)

# ============================================================
# STEP 1 — Download GSE89632 from GEO
# ============================================================
cat("\n=== STEP 1: Download GSE89632 from GEO ===\n")

gse  <- getGEO("GSE89632", GSEMatrix = TRUE, destdir = "data/")
eset <- gse[[1]]

cat("ExpressionSet loaded:", nrow(eset), "probes x", ncol(eset), "samples\n")

# ============================================================
# STEP 2 — Extract metadata and define groups
# ============================================================
cat("\n=== STEP 2: Extract metadata and sample groups ===\n")

meta <- pData(eset)

# Diagnosis is in characteristics_ch1.1: "diagnosis: HC / SS / NASH"
meta$diagnosis <- sub("diagnosis: ", "", meta$"characteristics_ch1.1")
meta$fibrosis  <- sub("fibrosis \\(stage\\): ", "", meta$"fibrosis (stage):ch1")
meta$nas_score <- meta$"nafld activity score:ch1"

cat("Sample groups:\n")
print(table(meta$diagnosis))
cat("\nSample names match expression matrix:", all(rownames(meta) == colnames(exprs(eset))), "\n")

# ============================================================
# STEP 3 — Extract expression matrix and normalize
# ============================================================
cat("\n=== STEP 3: Extract expression matrix and normalise ===\n")

ex_raw <- exprs(eset)
cat("Raw expression matrix:", nrow(ex_raw), "probes x", ncol(ex_raw), "samples\n")
cat("Value range:", round(min(ex_raw), 2), "to", round(max(ex_raw), 2), "\n")
cat("Data appears log2-transformed (range ~7-16).\n")

# Between-array quantile normalisation (aligns array distributions)
# The data is already log2 — we apply normalizeBetweenArrays to remove
# residual between-array technical variation
ex_norm <- normalizeBetweenArrays(ex_raw, method = "quantile")
cat("After quantile normalisation — range:", round(min(ex_norm), 2),
    "to", round(max(ex_norm), 2), "\n")

# ============================================================
# STEP 4 — Probe filtering: keep expressed, annotated probes
# ============================================================
cat("\n=== STEP 4: Filter probes ===\n")

fd <- fData(eset)

# Keep probes with a gene symbol
has_symbol <- fd$Symbol != "" & !is.na(fd$Symbol)
cat("Probes with gene symbol:", sum(has_symbol), "/", nrow(fd), "\n")

# Keep probes above background (expression > 7.5 in at least 15% of samples)
min_samples <- ceiling(0.15 * ncol(ex_norm))
expressed   <- rowSums(ex_norm > 7.5) >= min_samples
cat("Probes above background threshold:", sum(expressed), "\n")

keep        <- has_symbol & expressed
ex_filt     <- ex_norm[keep, ]
fd_filt     <- fd[keep, ]
cat("Probes retained for analysis:", nrow(ex_filt), "\n")

# For multi-probe genes: keep probe with highest mean expression
fd_filt$mean_expr <- rowMeans(ex_filt)
fd_filt$probe_id  <- rownames(fd_filt)

# Pick best probe per symbol
best_probes <- fd_filt %>%
  group_by(Symbol) %>%
  slice_max(order_by = mean_expr, n = 1, with_ties = FALSE) %>%
  pull(probe_id)

ex_gene  <- ex_filt[best_probes, ]
fd_gene  <- fd_filt[best_probes, ]
rownames(ex_gene) <- fd_gene$Symbol
cat("Unique genes after probe collapse:", nrow(ex_gene), "\n")

# ============================================================
# STEP 5 — limma: HC vs NASH differential expression
# ============================================================
cat("\n=== STEP 5: limma — HC vs NASH ===\n")

# Set up design matrix: HC is reference
diagnosis_f <- factor(meta$diagnosis, levels = c("HC", "SS", "NASH"))
design <- model.matrix(~ 0 + diagnosis_f)
colnames(design) <- c("HC", "SS", "NASH")

fit <- lmFit(ex_gene, design)

# Contrast: NASH vs HC
contrast_mat <- makeContrasts(NASH_vs_HC = NASH - HC,
                               SS_vs_HC   = SS   - HC,
                               levels     = design)
fit2  <- contrasts.fit(fit, contrast_mat)
fit2  <- eBayes(fit2)

# Extract full table for NASH vs HC
tt_nash <- topTable(fit2, coef = "NASH_vs_HC", number = Inf, sort.by = "none")
tt_nash$gene <- rownames(tt_nash)
tt_nash <- tt_nash %>% arrange(adj.P.Val)

cat("\nlimma summary (NASH vs HC, adj.P.Val < 0.05):\n")
sig <- tt_nash %>% filter(adj.P.Val < 0.05)
cat("  Up in NASH  (logFC > 0):", sum(sig$logFC > 0), "\n")
cat("  Down in NASH (logFC < 0):", sum(sig$logFC < 0), "\n")

write.csv(tt_nash, "results/limma_nash_vs_hc.csv", row.names = FALSE)
cat("Full limma results saved.\n")

# Also save SS vs HC for completeness
tt_ss <- topTable(fit2, coef = "SS_vs_HC", number = Inf, sort.by = "none")
tt_ss$gene <- rownames(tt_ss)
write.csv(tt_ss, "results/limma_ss_vs_hc.csv", row.names = FALSE)

# ============================================================
# STEP 6 — Volcano plot: NASH vs HC
# ============================================================
cat("\n=== STEP 6: Volcano plot ===\n")

volcano_df <- tt_nash %>%
  filter(!is.na(adj.P.Val)) %>%
  mutate(
    sig            = case_when(
      adj.P.Val < 0.05 & logFC >  0.5 ~ "Up in NASH",
      adj.P.Val < 0.05 & logFC < -0.5 ~ "Down in NASH",
      TRUE ~ "NS"),
    neg_log10_padj = -log10(adj.P.Val)
  )

top_labels <- volcano_df %>%
  filter(sig != "NS") %>%
  slice_max(order_by = neg_log10_padj, n = 25)

cat("Up in NASH   (adj.P<0.05, logFC>0.5):", sum(volcano_df$sig == "Up in NASH"), "\n")
cat("Down in NASH (adj.P<0.05, logFC<-0.5):", sum(volcano_df$sig == "Down in NASH"), "\n")

p_vol <- ggplot(volcano_df, aes(x = logFC, y = neg_log10_padj, colour = sig)) +
  geom_point(alpha = 0.45, size = 0.7) +
  geom_text_repel(data = top_labels, aes(label = gene),
                  size = 2.4, max.overlaps = 20,
                  segment.size = 0.3, segment.alpha = 0.6) +
  scale_colour_manual(values = c("Up in NASH"   = "#C62828",
                                  "Down in NASH" = "#1565C0",
                                  "NS"           = "grey70")) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed",
             colour = "black", linewidth = 0.4) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             colour = "black", linewidth = 0.4) +
  labs(title    = "Volcano Plot: NASH vs Healthy Controls",
       subtitle = "GSE89632 | 43 samples (HC n=24, NASH n=19) | limma microarray",
       x        = "log2 Fold Change (NASH / HC)",
       y        = expression(-log[10](p[adj])),
       colour   = NULL) +
  theme_bw(base_size = 12) +
  theme(plot.title    = element_text(face = "bold"),
        legend.position = "top")

ggsave("plots/volcano_nash_vs_hc.pdf", p_vol, width = 8, height = 6.5)
ggsave("plots/volcano_nash_vs_hc.png", p_vol, width = 8, height = 6.5, dpi = 150)
cat("Volcano plot saved.\n")

# ============================================================
# STEP 7 — Top 20 DEGs: NASH vs HC
# ============================================================
cat("\n=== STEP 7: Top 20 DEGs ===\n")

top20_nash <- tt_nash %>%
  filter(adj.P.Val < 0.05) %>%
  slice_max(order_by = abs(logFC), n = 20) %>%
  dplyr::select(gene, logFC, AveExpr, t, P.Value, adj.P.Val)

print(as.data.frame(top20_nash))
write.csv(top20_nash, "results/top20_degs_nash_vs_hc.csv", row.names = FALSE)

p_top20 <- top20_nash %>%
  mutate(direction = ifelse(logFC > 0, "Up in NASH", "Down in NASH"),
         gene      = reorder(gene, logFC)) %>%
  ggplot(aes(x = logFC, y = gene, fill = direction)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = c("Up in NASH"   = "#C62828",
                                "Down in NASH" = "#1565C0")) +
  labs(title    = "Top 20 DEGs: NASH vs Healthy Controls",
       subtitle = "Ranked by absolute log2 Fold Change (adj.P < 0.05)",
       x = "log2 Fold Change", y = NULL, fill = NULL) +
  theme_bw(base_size = 11) +
  theme(plot.title    = element_text(face = "bold"),
        legend.position = "top")

ggsave("plots/top20_degs_nash_vs_hc.pdf", p_top20, width = 7, height = 6)
ggsave("plots/top20_degs_nash_vs_hc.png", p_top20, width = 7, height = 6, dpi = 150)
cat("Top 20 DEG plot saved.\n")

# ============================================================
# STEP 8 — Cross-platform validation: compare with Day 1 DEGs
# ============================================================
cat("\n=== STEP 8: Cross-platform validation with Day 1 RNA-seq ===\n")

day1_file <- "../week4-day1-nafld-bulk-rnaseq/results/deseq2_nafld_vs_normal.csv"
if (file.exists(day1_file)) {
  day1 <- read.csv(day1_file)

  # Day 1 significant genes
  day1_up   <- day1 %>% filter(!is.na(padj), padj < 0.05, log2FoldChange > 1) %>%
    pull(gene_symbol) %>% na.omit()
  day1_down <- day1 %>% filter(!is.na(padj), padj < 0.05, log2FoldChange < -1) %>%
    pull(gene_symbol) %>% na.omit()

  # Day 2 significant genes
  day2_up   <- tt_nash %>% filter(adj.P.Val < 0.05, logFC > 0.5) %>% pull(gene)
  day2_down <- tt_nash %>% filter(adj.P.Val < 0.05, logFC < -0.5) %>% pull(gene)

  # Overlap
  overlap_up   <- intersect(day1_up,   day2_up)
  overlap_down <- intersect(day1_down, day2_down)

  cat("\nDay 1 (RNA-seq)    — Up:", length(day1_up),
      "| Down:", length(day1_down), "\n")
  cat("Day 2 (microarray) — Up:", length(day2_up),
      "| Down:", length(day2_down), "\n")
  cat("\nGenes UP in NAFLD/NASH in BOTH datasets:", length(overlap_up), "\n")
  if (length(overlap_up) > 0) print(overlap_up)
  cat("\nGenes DOWN in NAFLD/NASH in BOTH datasets:", length(overlap_down), "\n")
  if (length(overlap_down) > 0) print(overlap_down)

  # Save overlap
  overlap_df <- bind_rows(
    data.frame(gene = overlap_up,   direction = "Up",   stringsAsFactors = FALSE),
    data.frame(gene = overlap_down, direction = "Down", stringsAsFactors = FALSE)
  )
  write.csv(overlap_df, "results/cross_platform_overlap.csv", row.names = FALSE)

  # Scatter plot: Day 1 LFC vs Day 2 LFC for all shared genes
  all_genes <- intersect(
    day1 %>% filter(!is.na(gene_symbol)) %>% pull(gene_symbol),
    tt_nash$gene
  )
  scatter_df <- inner_join(
    day1 %>% dplyr::select(gene_symbol, log2FoldChange) %>%
      rename(gene = gene_symbol, lfc_day1 = log2FoldChange) %>%
      filter(!is.na(gene)),
    tt_nash %>% dplyr::select(gene, logFC) %>%
      rename(lfc_day2 = logFC),
    by = "gene"
  )
  cat("\nGenes in scatter comparison:", nrow(scatter_df), "\n")

  # Pearson correlation of LFCs
  r <- cor(scatter_df$lfc_day1, scatter_df$lfc_day2, use = "complete.obs")
  cat("Pearson correlation of log2FC between platforms:", round(r, 3), "\n")

  # Label shared significant genes
  scatter_df <- scatter_df %>%
    mutate(shared_sig = gene %in% c(overlap_up, overlap_down))

  p_scatter <- ggplot(scatter_df, aes(x = lfc_day1, y = lfc_day2)) +
    geom_point(aes(colour = shared_sig), alpha = 0.3, size = 0.6) +
    geom_smooth(method = "lm", colour = "black", linewidth = 0.8, se = FALSE) +
    geom_text_repel(
      data = scatter_df %>% filter(shared_sig) %>%
        slice_max(order_by = abs(lfc_day1), n = 20),
      aes(label = gene), size = 2.2, max.overlaps = 20,
      segment.size = 0.3
    ) +
    scale_colour_manual(values = c("FALSE" = "grey70", "TRUE" = "#B71C1C"),
                        labels = c("Not significant", "Shared DEG (both platforms)"),
                        guide  = guide_legend(title = NULL)) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5,
             label = paste0("r = ", round(r, 3)), size = 4, fontface = "bold") +
    labs(title    = "Cross-Platform Validation: Day 1 vs Day 2 log2FC",
         subtitle = "RNA-seq (GSE162694) vs Microarray (GSE89632) — same direction = validation",
         x        = "log2FC — RNA-seq (DESeq2, Day 1)",
         y        = "log2FC — Microarray (limma, Day 2)") +
    theme_bw(base_size = 12) +
    theme(plot.title    = element_text(face = "bold"),
          legend.position = "bottom")

  ggsave("plots/cross_platform_scatter.pdf", p_scatter, width = 8, height = 7)
  ggsave("plots/cross_platform_scatter.png", p_scatter, width = 8, height = 7, dpi = 150)
  cat("Cross-platform scatter plot saved.\n")
} else {
  cat("Day 1 results not found at:", day1_file, "— skipping validation plot.\n")
}

# ============================================================
# STEP 9 — Boxplots of key validated genes across HC / SS / NASH
# ============================================================
cat("\n=== STEP 9: Expression boxplots for key genes ===\n")

# Pick top validated genes + known NAFLD markers
key_genes <- c("CXCL8", "NR4A1", "TNFRSF12A", "MPO",
                "SPP1", "COL1A1", "ACTA2", "CCL2",
                "CYP2C19", "CYP3A4", "CYP1A2", "UGT2B7")
key_genes <- key_genes[key_genes %in% rownames(ex_gene)]

if (length(key_genes) > 0) {
  box_df <- ex_gene[key_genes, , drop = FALSE] %>%
    as.data.frame() %>%
    tibble::rownames_to_column("gene") %>%
    tidyr::pivot_longer(-gene, names_to = "sample", values_to = "expr") %>%
    left_join(meta %>% dplyr::select(geo_accession, diagnosis) %>%
                rename(sample = geo_accession),
              by = "sample") %>%
    mutate(diagnosis = factor(diagnosis, levels = c("HC", "SS", "NASH")))

  p_box <- ggplot(box_df, aes(x = diagnosis, y = expr, fill = diagnosis)) +
    geom_boxplot(outlier.size = 0.5, width = 0.6) +
    geom_jitter(width = 0.15, size = 0.4, alpha = 0.5) +
    scale_fill_manual(values = c(HC = "#2E7D32", SS = "#F57F17", NASH = "#C62828")) +
    facet_wrap(~ gene, scales = "free_y", ncol = 4) +
    labs(title    = "Expression of Key NAFLD Genes Across Groups",
         subtitle = "GSE89632 | Illumina HumanHT-12 BeadChip | quantile normalised",
         x = NULL, y = "log2 Expression", fill = "Diagnosis") +
    theme_bw(base_size = 10) +
    theme(plot.title    = element_text(face = "bold"),
          strip.text    = element_text(face = "bold"),
          legend.position = "top",
          axis.text.x   = element_text(size = 8))

  ggsave("plots/boxplots_key_genes.pdf", p_box,
         width = 3 * ceiling(length(key_genes) / 3), height = 7)
  ggsave("plots/boxplots_key_genes.png", p_box,
         width = 3 * ceiling(length(key_genes) / 3), height = 7, dpi = 150)
  cat("Key gene boxplots saved for:", paste(key_genes, collapse = ", "), "\n")
}

# ============================================================
# STEP 10 — GSEA: KEGG + Hallmark (ranked by limma t-statistic)
# ============================================================
cat("\n=== STEP 10: GSEA ===\n")

# Use t-statistic as rank (more stable than LFC for microarray)
ranked_t <- setNames(tt_nash$t, tt_nash$gene)
ranked_t <- sort(ranked_t, decreasing = TRUE)
ranked_t <- ranked_t[!duplicated(names(ranked_t))]

# Map gene symbol -> Entrez for KEGG
entrez_map <- mapIds(org.Hs.eg.db, keys = names(ranked_t),
                     column = "ENTREZID", keytype = "SYMBOL",
                     multiVals = "first")
ranked_entrez <- setNames(ranked_t, entrez_map)
ranked_entrez <- ranked_entrez[!is.na(names(ranked_entrez))]
ranked_entrez <- ranked_entrez[!duplicated(names(ranked_entrez))]
cat("Genes in KEGG ranked list:", length(ranked_entrez), "\n")

## KEGG
cat("Running KEGG GSEA...\n")
gsea_kegg <- gseKEGG(geneList     = ranked_entrez,
                      organism     = "hsa",
                      minGSSize    = 15,
                      maxGSSize    = 500,
                      pvalueCutoff = 0.25,
                      verbose      = FALSE,
                      seed         = 42)
cat("Significant KEGG pathways:", nrow(gsea_kegg), "\n")

if (nrow(gsea_kegg) > 0) {
  write.csv(as.data.frame(gsea_kegg), "results/gsea_kegg_results.csv",
            row.names = FALSE)
  p_kd <- dotplot(gsea_kegg, showCategory = 20, split = ".sign",
                  title = "KEGG GSEA: NASH vs HC (GSE89632)") +
    facet_grid(. ~ .sign) + theme_bw(base_size = 10) +
    theme(plot.title  = element_text(face = "bold"),
          axis.text.y = element_text(size = 8))
  ggsave("plots/gsea_kegg_dotplot.pdf", p_kd, width = 13, height = 8)
  ggsave("plots/gsea_kegg_dotplot.png", p_kd, width = 13, height = 8, dpi = 150)
  if (nrow(gsea_kegg) >= 5) {
    p_kr <- ridgeplot(gsea_kegg, showCategory = 15) +
      labs(title = "KEGG GSEA Ridge Plot") + theme_bw(base_size = 10) +
      theme(plot.title = element_text(face = "bold"), axis.text.y = element_text(size = 8))
    ggsave("plots/gsea_kegg_ridgeplot.pdf", p_kr, width = 11, height = 8)
    ggsave("plots/gsea_kegg_ridgeplot.png", p_kr, width = 11, height = 8, dpi = 150)
  }
  cat("KEGG GSEA plots saved.\n")
}

## Hallmark
cat("\nRunning Hallmark GSEA...\n")
h_t2g <- msigdbr(species = "Homo sapiens", collection = "H") %>%
  dplyr::select(gs_name, gene_symbol)

gsea_h <- GSEA(geneList     = ranked_t,
               TERM2GENE    = h_t2g,
               minGSSize    = 15,
               maxGSSize    = 500,
               pvalueCutoff = 0.25,
               verbose      = FALSE,
               seed         = 42)
cat("Significant Hallmark sets:", nrow(gsea_h), "\n")

if (nrow(gsea_h) > 0) {
  write.csv(as.data.frame(gsea_h), "results/gsea_hallmark_results.csv",
            row.names = FALSE)
  p_hd <- dotplot(gsea_h, showCategory = 20, split = ".sign",
                  title = "Hallmark GSEA: NASH vs HC (GSE89632)") +
    facet_grid(. ~ .sign) + theme_bw(base_size = 10) +
    theme(plot.title  = element_text(face = "bold"),
          axis.text.y = element_text(size = 7))
  ggsave("plots/gsea_hallmark_dotplot.pdf", p_hd, width = 14, height = 8)
  ggsave("plots/gsea_hallmark_dotplot.png", p_hd, width = 14, height = 8, dpi = 150)
  if (nrow(gsea_h) >= 5) {
    p_hr <- ridgeplot(gsea_h, showCategory = 15) +
      labs(title = "Hallmark GSEA Ridge Plot") + theme_bw(base_size = 10) +
      theme(plot.title = element_text(face = "bold"), axis.text.y = element_text(size = 7))
    ggsave("plots/gsea_hallmark_ridgeplot.pdf", p_hr, width = 12, height = 8)
    ggsave("plots/gsea_hallmark_ridgeplot.png", p_hr, width = 12, height = 8, dpi = 150)
  }
  cat("Hallmark GSEA plots saved.\n")
}

# ============================================================
# STEP 11 — PCA
# ============================================================
cat("\n=== STEP 11: PCA ===\n")

pca_res <- prcomp(t(ex_gene), scale. = TRUE)
pca_var <- round(100 * summary(pca_res)$importance[2, 1:2], 1)

pca_df <- data.frame(
  PC1       = pca_res$x[, 1],
  PC2       = pca_res$x[, 2],
  diagnosis = factor(meta$diagnosis, levels = c("HC", "SS", "NASH")),
  sample    = rownames(pca_res$x)
)

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, colour = diagnosis)) +
  geom_point(size = 2.5, alpha = 0.85) +
  scale_colour_manual(values = c(HC = "#2E7D32", SS = "#F57F17", NASH = "#C62828")) +
  labs(title    = "PCA: GSE89632 Microarray (Quantile Normalised)",
       subtitle = "HC n=24 | SS n=20 | NASH n=19",
       x        = paste0("PC1 (", pca_var[1], "%)"),
       y        = paste0("PC2 (", pca_var[2], "%)"),
       colour   = NULL) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), legend.position = "top")

ggsave("plots/pca_groups.pdf", p_pca, width = 7, height = 5.5)
ggsave("plots/pca_groups.png", p_pca, width = 7, height = 5.5, dpi = 150)
cat("PCA saved. PC1:", pca_var[1], "% PC2:", pca_var[2], "%\n")

# ============================================================
# Summary
# ============================================================
cat("\n=== All analyses complete ===\n")
for (f in c(list.files("plots",   full.names = TRUE),
            list.files("results", full.names = TRUE))) cat(" ", f, "\n")
