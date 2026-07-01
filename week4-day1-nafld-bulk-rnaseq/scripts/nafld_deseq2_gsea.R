## ============================================================
## Week 4 Day 1 — NAFLD Bulk RNA-seq: DESeq2 + GSEA
## Dataset: GSE162694 (human liver, normal vs NAFLD F0-F4)
## 143 samples | Ensembl gene IDs | fibrosis stage 0-4 + normal
## ============================================================

suppressPackageStartupMessages({
  library(GEOquery)
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(msigdbr)
  library(dplyr)
  library(tidyr)
  library(BiocParallel)
})

set.seed(42)
register(MulticoreParam(4))

# ============================================================
# STEP 1 — Download metadata from GEO
# ============================================================
cat("\n=== STEP 1: Download GSE162694 metadata from GEO ===\n")

gse <- getGEO("GSE162694", GSEMatrix = TRUE, destdir = "data/")
meta <- pData(gse[[1]])
cat("Loaded metadata:", nrow(meta), "samples\n")

# Extract sample label used in count matrix (second token in title field)
# title example: "nash1_F0 548nash1"  -> sample_id = "548nash1"
meta$sample_id <- sub(".* ", "", meta$title)

# Fibrosis stage from characteristics_ch1.3 = "fibrosis stage: X"
meta$fibrosis_raw <- sub("fibrosis stage: ", "", meta$"characteristics_ch1.3")

# Condition: "normal liver histology" → Normal, everything else → NAFLD
meta$condition <- ifelse(
  meta$fibrosis_raw == "normal liver histology", "Normal", "NAFLD"
)

cat("\nSample counts by condition:\n")
print(table(meta$condition))
cat("\nSample counts by fibrosis stage:\n")
print(table(meta$fibrosis_raw))

# ============================================================
# STEP 2 — Load supplementary count matrix (CSV, comma-separated)
# ============================================================
cat("\n=== STEP 2: Load count matrix ===\n")

csv_file <- "data/GSE162694/GSE162694_raw_counts.csv"

# File was already decompressed by GEOquery's getGEOSuppFiles
if (!file.exists(csv_file)) {
  getGEOSuppFiles("GSE162694", baseDir = "data/", makeDirectory = TRUE)
}

counts_raw <- read.csv(csv_file, row.names = 1, check.names = FALSE)
cat("Raw count matrix:", nrow(counts_raw), "genes x", ncol(counts_raw), "samples\n")
cat("First 4 gene IDs:", head(rownames(counts_raw), 4), "\n")
cat("First 4 sample names:", head(colnames(counts_raw), 4), "\n")

# ============================================================
# STEP 3 — Align samples between count matrix and metadata
# ============================================================
cat("\n=== STEP 3: Align metadata to count matrix ===\n")

rownames(meta) <- meta$sample_id
common_samples  <- intersect(colnames(counts_raw), meta$sample_id)
cat("Samples in common:", length(common_samples), "\n")

counts_mat <- as.matrix(counts_raw[, common_samples])
meta_sub   <- meta[common_samples, ]

storage.mode(counts_mat) <- "integer"

# ============================================================
# STEP 4 — Filter low-count genes
# ============================================================
cat("\n=== STEP 4: Filter low-count genes ===\n")
cat("Before filter:", nrow(counts_mat), "genes\n")
keep         <- rowSums(counts_mat > 10) >= 3
counts_filt  <- counts_mat[keep, ]
cat("After filter (>10 counts in >=3 samples):", nrow(counts_filt), "genes\n")

# ============================================================
# STEP 5 — DESeq2: Normal vs NAFLD (all stages combined)
# ============================================================
cat("\n=== STEP 5: DESeq2 — Normal vs NAFLD ===\n")

col_data <- data.frame(
  condition  = factor(meta_sub$condition, levels = c("Normal", "NAFLD")),
  fibrosis   = factor(meta_sub$fibrosis_raw),
  row.names  = colnames(counts_filt)
)

dds <- DESeqDataSetFromMatrix(
  countData = counts_filt,
  colData   = col_data,
  design    = ~ condition
)

dds <- DESeq(dds, parallel = TRUE)
res <- results(dds, contrast  = c("condition", "NAFLD", "Normal"),
               alpha = 0.05)
res_shrunk <- lfcShrink(dds, coef = "condition_NAFLD_vs_Normal",
                         type = "apeglm")

cat("\nDESeq2 result summary (NAFLD vs Normal):\n")
summary(res_shrunk)

# Convert to data frame, add gene symbol from Ensembl ID
res_df <- as.data.frame(res_shrunk) %>%
  tibble::rownames_to_column("ensembl_id") %>%
  arrange(padj)

# Map Ensembl IDs → gene symbols
res_df$gene_symbol <- mapIds(
  org.Hs.eg.db,
  keys     = res_df$ensembl_id,
  column   = "SYMBOL",
  keytype  = "ENSEMBL",
  multiVals = "first"
)

write.csv(res_df, "results/deseq2_nafld_vs_normal.csv", row.names = FALSE)
cat("Full DESeq2 results saved to results/deseq2_nafld_vs_normal.csv\n")

# ============================================================
# STEP 6 — Volcano plot
# ============================================================
cat("\n=== STEP 6: Volcano plot ===\n")

volcano_df <- res_df %>%
  filter(!is.na(padj), !is.na(log2FoldChange)) %>%
  mutate(
    label = ifelse(!is.na(gene_symbol), gene_symbol, ensembl_id),
    sig   = case_when(
      padj < 0.05 & log2FoldChange >  1 ~ "Up in NAFLD",
      padj < 0.05 & log2FoldChange < -1 ~ "Down in NAFLD",
      TRUE ~ "NS"
    ),
    neg_log10_padj = -log10(padj)
  )

top_labels <- volcano_df %>%
  filter(sig != "NS") %>%
  slice_max(order_by = neg_log10_padj, n = 25)

cat("Up in NAFLD  (padj<0.05, LFC>1):", sum(volcano_df$sig == "Up in NAFLD"), "\n")
cat("Down in NAFLD (padj<0.05, LFC<-1):", sum(volcano_df$sig == "Down in NAFLD"), "\n")

p_volcano <- ggplot(volcano_df,
                    aes(x = log2FoldChange, y = neg_log10_padj, colour = sig)) +
  geom_point(alpha = 0.45, size = 0.7) +
  geom_text_repel(data = top_labels, aes(label = label),
                  size = 2.4, max.overlaps = 20,
                  segment.size = 0.3, segment.alpha = 0.6) +
  scale_colour_manual(values = c("Up in NAFLD"   = "#C62828",
                                  "Down in NAFLD" = "#1565C0",
                                  "NS"            = "grey70")) +
  geom_vline(xintercept = c(-1, 1),  linetype = "dashed",
             colour = "black", linewidth = 0.4) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             colour = "black", linewidth = 0.4) +
  labs(title    = "Volcano Plot: NAFLD vs Normal Liver",
       subtitle = "GSE162694 | 143 samples | DESeq2 + apeglm shrinkage",
       x        = "log2 Fold Change (NAFLD / Normal)",
       y        = expression(-log[10](p[adj])),
       colour   = NULL) +
  theme_bw(base_size = 12) +
  theme(plot.title    = element_text(face = "bold"),
        legend.position = "top")

ggsave("plots/volcano_nafld_vs_normal.pdf", p_volcano, width = 8, height = 6.5)
ggsave("plots/volcano_nafld_vs_normal.png", p_volcano, width = 8, height = 6.5, dpi = 150)
cat("Volcano plot saved.\n")

# ============================================================
# STEP 7 — Top 20 DEG table and bar chart
# ============================================================
cat("\n=== STEP 7: Top 20 DEGs ===\n")

top20 <- res_df %>%
  filter(!is.na(padj), padj < 0.05) %>%
  slice_max(order_by = abs(log2FoldChange), n = 20) %>%
  select(gene_symbol, ensembl_id, baseMean,
         log2FoldChange, lfcSE, pvalue, padj)

cat("Top 20 DEGs by absolute LFC:\n")
print(as.data.frame(top20))
write.csv(top20, "results/top20_degs.csv", row.names = FALSE)

p_top20 <- top20 %>%
  mutate(
    label     = ifelse(!is.na(gene_symbol), gene_symbol, ensembl_id),
    direction = ifelse(log2FoldChange > 0, "Up in NAFLD", "Down in NAFLD"),
    label     = reorder(label, log2FoldChange)
  ) %>%
  ggplot(aes(x = log2FoldChange, y = label, fill = direction)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = c("Up in NAFLD"   = "#C62828",
                                "Down in NAFLD" = "#1565C0")) +
  labs(title    = "Top 20 DEGs: NAFLD vs Normal",
       subtitle = "Ranked by absolute log2 Fold Change (padj < 0.05)",
       x = "log₂ Fold Change", y = NULL, fill = NULL) +
  theme_bw(base_size = 11) +
  theme(plot.title    = element_text(face = "bold"),
        legend.position = "top")

ggsave("plots/top20_degs_barplot.pdf", p_top20, width = 7, height = 6)
ggsave("plots/top20_degs_barplot.png", p_top20, width = 7, height = 6, dpi = 150)
cat("Top 20 DEG bar chart saved.\n")

# ============================================================
# STEP 8 — PCA plot (VST-transformed)
# ============================================================
cat("\n=== STEP 8: PCA plot ===\n")
vsd      <- vst(dds, blind = TRUE)
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
pca_var  <- round(100 * attr(pca_data, "percentVar"), 1)

p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2, colour = condition)) +
  geom_point(size = 2.5, alpha = 0.85) +
  scale_colour_manual(values = c(Normal = "#2E7D32", NAFLD = "#E65100")) +
  labs(title    = "PCA — VST-transformed Counts",
       subtitle = "GSE162694: Normal vs NAFLD (all fibrosis stages)",
       x        = paste0("PC1 (", pca_var[1], "%)"),
       y        = paste0("PC2 (", pca_var[2], "%)"),
       colour   = NULL) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "top")

ggsave("plots/pca_nafld_vs_normal.pdf", p_pca, width = 7, height = 5.5)
ggsave("plots/pca_nafld_vs_normal.png", p_pca, width = 7, height = 5.5, dpi = 150)
cat("PCA plot saved.\n")

# ============================================================
# STEP 9 — GSEA: KEGG + Hallmark
# ============================================================
cat("\n=== STEP 9: GSEA ===\n")

# Build ranked gene list using Entrez IDs (for KEGG) and gene symbols (for Hallmark)
# Map Ensembl → Entrez
ranked_df <- res_df %>%
  filter(!is.na(log2FoldChange)) %>%
  mutate(
    entrez = mapIds(org.Hs.eg.db,
                    keys     = ensembl_id,
                    column   = "ENTREZID",
                    keytype  = "ENSEMBL",
                    multiVals = "first")
  )

## — KEGG GSEA (ranked by LFC, Entrez IDs) —
kegg_input <- ranked_df %>%
  filter(!is.na(entrez)) %>%
  arrange(desc(log2FoldChange))
ranked_entrez <- setNames(kegg_input$log2FoldChange, kegg_input$entrez)
ranked_entrez <- ranked_entrez[!duplicated(names(ranked_entrez))]

cat("Genes in KEGG ranked list:", length(ranked_entrez), "\n")
cat("Running KEGG GSEA...\n")

gsea_kegg <- gseKEGG(
  geneList     = ranked_entrez,
  organism     = "hsa",
  minGSSize    = 15,
  maxGSSize    = 500,
  pvalueCutoff = 0.25,
  verbose      = FALSE,
  seed         = 42
)

if (nrow(gsea_kegg) > 0) {
  cat("Significant KEGG pathways:", nrow(gsea_kegg), "\n")
  write.csv(as.data.frame(gsea_kegg), "results/gsea_kegg_results.csv", row.names = FALSE)

  p_kegg_dot <- dotplot(gsea_kegg, showCategory = 20,
                         split = ".sign",
                         title = "KEGG GSEA — NAFLD vs Normal") +
    facet_grid(. ~ .sign) +
    theme_bw(base_size = 10) +
    theme(plot.title    = element_text(face = "bold"),
          axis.text.y   = element_text(size = 8),
          strip.text    = element_text(face = "bold"))
  ggsave("plots/gsea_kegg_dotplot.pdf",  p_kegg_dot, width = 13, height = 8)
  ggsave("plots/gsea_kegg_dotplot.png",  p_kegg_dot, width = 13, height = 8, dpi = 150)

  if (nrow(gsea_kegg) >= 5) {
    p_kegg_ridge <- ridgeplot(gsea_kegg, showCategory = 15) +
      labs(title = "KEGG GSEA — Ridge Plot") +
      theme_bw(base_size = 10) +
      theme(plot.title  = element_text(face = "bold"),
            axis.text.y = element_text(size = 8))
    ggsave("plots/gsea_kegg_ridgeplot.pdf", p_kegg_ridge, width = 11, height = 8)
    ggsave("plots/gsea_kegg_ridgeplot.png", p_kegg_ridge, width = 11, height = 8, dpi = 150)
  }
  cat("KEGG GSEA plots saved.\n")
} else {
  cat("No significant KEGG pathways found.\n")
}

## — Hallmark GSEA (ranked by LFC, gene symbols) —
hallmark_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, gene_symbol)

hallmark_input <- res_df %>%
  filter(!is.na(log2FoldChange), !is.na(gene_symbol)) %>%
  arrange(desc(log2FoldChange))
ranked_sym <- setNames(hallmark_input$log2FoldChange, hallmark_input$gene_symbol)
ranked_sym <- ranked_sym[!duplicated(names(ranked_sym))]

cat("\nGenes in Hallmark ranked list:", length(ranked_sym), "\n")
cat("Running Hallmark GSEA...\n")

gsea_hallmark <- GSEA(
  geneList     = ranked_sym,
  TERM2GENE    = hallmark_t2g,
  minGSSize    = 15,
  maxGSSize    = 500,
  pvalueCutoff = 0.25,
  verbose      = FALSE,
  seed         = 42
)

if (nrow(gsea_hallmark) > 0) {
  cat("Significant Hallmark gene sets:", nrow(gsea_hallmark), "\n")
  write.csv(as.data.frame(gsea_hallmark),
            "results/gsea_hallmark_results.csv", row.names = FALSE)

  p_hm_dot <- dotplot(gsea_hallmark, showCategory = 20,
                       split = ".sign",
                       title = "Hallmark GSEA — NAFLD vs Normal") +
    facet_grid(. ~ .sign) +
    theme_bw(base_size = 10) +
    theme(plot.title    = element_text(face = "bold"),
          axis.text.y   = element_text(size = 7),
          strip.text    = element_text(face = "bold"))
  ggsave("plots/gsea_hallmark_dotplot.pdf", p_hm_dot, width = 14, height = 8)
  ggsave("plots/gsea_hallmark_dotplot.png", p_hm_dot, width = 14, height = 8, dpi = 150)

  if (nrow(gsea_hallmark) >= 5) {
    p_hm_ridge <- ridgeplot(gsea_hallmark, showCategory = 15) +
      labs(title = "Hallmark GSEA — Ridge Plot") +
      theme_bw(base_size = 10) +
      theme(plot.title  = element_text(face = "bold"),
            axis.text.y = element_text(size = 7))
    ggsave("plots/gsea_hallmark_ridgeplot.pdf", p_hm_ridge, width = 12, height = 8)
    ggsave("plots/gsea_hallmark_ridgeplot.png", p_hm_ridge, width = 12, height = 8, dpi = 150)
  }
  cat("Hallmark GSEA plots saved.\n")
} else {
  cat("No significant Hallmark gene sets found.\n")
}

# ============================================================
# Summary
# ============================================================
cat("\n============================================================\n")
cat("All analyses complete.\n")
cat("Plots saved to:   plots/\n")
cat("Results saved to: results/\n")
cat("Files generated:\n")
for (f in c(list.files("plots", full.names = TRUE),
            list.files("results", full.names = TRUE))) {
  cat("  ", f, "\n")
}
