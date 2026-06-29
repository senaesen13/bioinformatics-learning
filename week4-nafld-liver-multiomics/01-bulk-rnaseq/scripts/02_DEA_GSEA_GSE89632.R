## ============================================================
## GSE89632  –  Differential Expression + GSEA
## Comparisons: NASH vs HC  |  SS vs HC
## Platform: Illumina HumanHT-12 WG-DASL v4 (microarray)
## ============================================================

suppressPackageStartupMessages({
  library(Biobase)
  library(limma)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(msigdbr)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(RColorBrewer)
  library(pheatmap)
  library(scales)
})

# ── Paths ──────────────────────────────────────────────────────────────────────
res_dir  <- "01-bulk-rnaseq/results"
plot_dir <- "01-bulk-rnaseq/plots"
dir.create(res_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# ══════════════════════════════════════════════════════════════════════════════
# STEP 1 – Load data and inspect
# ══════════════════════════════════════════════════════════════════════════════
# The ExpressionSet was saved during the download step.
# It contains: exprs() = log2 intensity matrix, pData() = sample metadata,
#              fData() = probe annotation (ILMN IDs → gene symbols).
cat("[ STEP 1 ] Loading ExpressionSet ...\n")
eset <- readRDS(file.path(res_dir, "GSE89632_eset.rds"))
cat("  Probes :", nrow(eset), "\n")
cat("  Samples:", ncol(eset), "\n")

# Parse condition from characteristics_ch1.1  (e.g. "diagnosis: HC")
pd <- pData(eset)
pd$condition <- sub("diagnosis: ", "", pd[["characteristics_ch1.1"]])
pd$condition <- factor(pd$condition, levels = c("HC", "SS", "NASH"))
pData(eset) <- pd
cat("  Groups :", table(pd$condition), "\n\n")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 2 – Quality check: boxplot of raw expression distribution
# ══════════════════════════════════════════════════════════════════════════════
# Microarray data from GEO is typically already background-corrected and
# quantile-normalised by the submitter. We check this by plotting per-sample
# distributions – all boxes should sit at roughly the same level and spread.
cat("[ STEP 2 ] Plotting expression distributions ...\n")

ex <- exprs(eset)

png(file.path(plot_dir, "00_sample_distributions.png"),
    width = 1800, height = 600, res = 120)
par(mar = c(8, 4, 3, 1))
boxplot(ex, las = 2, cex.axis = 0.4, col = as.integer(pd$condition) + 1,
        main = "Per-sample log2 expression (colour = group)",
        ylab = "log2 intensity")
legend("topright", levels(pd$condition),
       fill = 2:(nlevels(pd$condition) + 1), bty = "n")
dev.off()
cat("  Saved: 00_sample_distributions.png\n\n")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 3 – Filter low-expression probes
# ══════════════════════════════════════════════════════════════════════════════
# We remove probes whose signal is near background in the majority of samples.
# Threshold: keep probes with median log2 intensity > 7.5 in at least one group.
# (The overall expression range is ~7.4–15.8; 7.5 sits just above the floor.)
cat("[ STEP 3 ] Filtering low-expression probes ...\n")

groups   <- pd$condition
medians  <- sapply(levels(groups), function(g) {
  rowMedians(ex[, groups == g, drop = FALSE])
})
keep     <- rowSums(medians > 7.5) >= 1
eset_f   <- eset[keep, ]
cat(sprintf("  Probes before: %d  |  after filter: %d\n\n",
            nrow(eset), nrow(eset_f)))

# ══════════════════════════════════════════════════════════════════════════════
# STEP 4 – Collapse probes → gene symbols (max-expression rule)
# ══════════════════════════════════════════════════════════════════════════════
# Multiple Illumina probes can target the same gene. We keep the probe with
# the highest average expression across all samples per gene symbol.
# Probes with no gene symbol are dropped.
cat("[ STEP 4 ] Collapsing probes to gene symbols ...\n")

fd        <- fData(eset_f)
symbols   <- fd$Symbol
valid     <- !is.na(symbols) & symbols != "" & symbols != "---"
eset_f    <- eset_f[valid, ]
symbols   <- symbols[valid]

ex_f      <- exprs(eset_f)
avg_expr  <- rowMeans(ex_f)
# For duplicates, keep highest-average probe
keep_idx  <- tapply(seq_along(symbols), symbols, function(idx) idx[which.max(avg_expr[idx])])
keep_idx  <- unlist(keep_idx)
ex_gene   <- ex_f[keep_idx, ]
rownames(ex_gene) <- symbols[keep_idx]

cat(sprintf("  Unique genes after collapsing: %d\n\n", nrow(ex_gene)))

# ══════════════════════════════════════════════════════════════════════════════
# STEP 5 – PCA to verify group separation
# ══════════════════════════════════════════════════════════════════════════════
cat("[ STEP 5 ] PCA plot ...\n")

pca      <- prcomp(t(ex_gene), scale. = TRUE)
pct      <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
pca_df   <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                        condition = pd$condition)

p_pca <- ggplot(pca_df, aes(PC1, PC2, colour = condition)) +
  geom_point(size = 3, alpha = 0.9) +
  scale_colour_manual(values = c(HC = "#2196F3", SS = "#FF9800", NASH = "#F44336")) +
  labs(title = "PCA – GSE89632 (NAFLD liver microarray)",
       x = paste0("PC1 (", pct[1], "%)"),
       y = paste0("PC2 (", pct[2], "%)")) +
  theme_bw(base_size = 13) +
  theme(legend.title = element_blank())

ggsave(file.path(plot_dir, "01_PCA.png"), p_pca, width = 6, height = 5, dpi = 150)
cat("  Saved: 01_PCA.png\n\n")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 6 – limma differential expression
# ══════════════════════════════════════════════════════════════════════════════
# limma uses a linear model + empirical Bayes shrinkage of variance estimates.
# This is the gold standard for microarray DE and also works well for
# log-normalised RNA-seq counts.
#
# Design matrix: one coefficient per group (no intercept), so contrasts are
# simple differences between group means.
cat("[ STEP 6 ] Fitting limma model ...\n")

design <- model.matrix(~ 0 + condition, data = pd)
colnames(design) <- levels(pd$condition)   # HC, SS, NASH

# Fit ordinary least-squares per gene
fit <- lmFit(ex_gene, design)

# Define contrasts
contrasts <- makeContrasts(
  NASH_vs_HC = NASH - HC,
  SS_vs_HC   = SS   - HC,
  levels     = design
)

fit2 <- contrasts.fit(fit, contrasts)

# Empirical Bayes moderation: borrows information across all genes to
# stabilise per-gene variance estimates (crucial for small n).
fit2 <- eBayes(fit2)

cat("  Model fit complete.\n\n")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 7 – Extract DEG tables
# ══════════════════════════════════════════════════════════════════════════════
# adj.P.Val = Benjamini-Hochberg FDR across all genes.
# We use |logFC| > 0.5 (i.e. ~1.4-fold) and FDR < 0.05 as significance cuts.
cat("[ STEP 7 ] Extracting DEG tables ...\n")

FC_CUT  <- 0.5
FDR_CUT <- 0.05

extract_deg <- function(coef_name) {
  tt <- topTable(fit2, coef = coef_name, number = Inf, sort.by = "P")
  tt$gene    <- rownames(tt)
  tt$sig     <- ifelse(tt$adj.P.Val < FDR_CUT & abs(tt$logFC) > FC_CUT,
                       ifelse(tt$logFC > 0, "Up", "Down"), "NS")
  tt
}

res_NASH <- extract_deg("NASH_vs_HC")
res_SS   <- extract_deg("SS_vs_HC")

# Summary counts
for (nm in c("NASH_vs_HC", "SS_vs_HC")) {
  r <- get(paste0("res_", sub("_vs_HC", "", nm)))
  cat(sprintf("  %s: %d up, %d down  (FDR<0.05, |logFC|>0.5)\n",
              nm,
              sum(r$sig == "Up"),
              sum(r$sig == "Down")))
}

# Save full tables
write.csv(res_NASH, file.path(res_dir, "DEG_NASH_vs_HC.csv"), row.names = FALSE)
write.csv(res_SS,   file.path(res_dir, "DEG_SS_vs_HC.csv"),   row.names = FALSE)
cat("  Saved: DEG_NASH_vs_HC.csv, DEG_SS_vs_HC.csv\n\n")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 8 – Volcano plots
# ══════════════════════════════════════════════════════════════════════════════
# Volcano plot: x = log2 fold-change, y = -log10(adjusted p-value).
# Each dot is a gene; significant DEGs (FDR<0.05, |FC|>0.5) are coloured.
# Top 15 genes by FDR are labelled to identify key players.
cat("[ STEP 8 ] Drawing volcano plots ...\n")

make_volcano <- function(df, title, file) {
  df$label <- ifelse(df$gene %in% head(df[df$sig != "NS", "gene"], 15),
                     df$gene, "")
  df$negLogP <- -log10(df$adj.P.Val + 1e-300)

  sig_pal <- c(Up = "#E41A1C", Down = "#377EB8", NS = "grey70")

  p <- ggplot(df, aes(logFC, negLogP, colour = sig, label = label)) +
    geom_point(size = 1, alpha = 0.6) +
    geom_hline(yintercept = -log10(FDR_CUT), linetype = "dashed", colour = "grey40") +
    geom_vline(xintercept = c(-FC_CUT, FC_CUT), linetype = "dashed", colour = "grey40") +
    geom_text_repel(size = 3, max.overlaps = 20,
                    box.padding = 0.4, segment.colour = "grey50") +
    scale_colour_manual(values = sig_pal,
                        labels = c(Up = "Up", Down = "Down", NS = "NS")) +
    labs(title = title,
         x = "log2 fold-change",
         y = "-log10(adj. p-value)",
         colour = NULL) +
    theme_bw(base_size = 13) +
    annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5, size = 3.5,
             label = sprintf("Up: %d   Down: %d", sum(df$sig=="Up"), sum(df$sig=="Down")))
  ggsave(file, p, width = 7, height = 6, dpi = 150)
}

make_volcano(res_NASH, "NASH vs Healthy Controls",
             file.path(plot_dir, "02_volcano_NASH_vs_HC.png"))
make_volcano(res_SS,   "Simple Steatosis vs Healthy Controls",
             file.path(plot_dir, "03_volcano_SS_vs_HC.png"))
cat("  Saved: 02_volcano_NASH_vs_HC.png, 03_volcano_SS_vs_HC.png\n\n")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 9 – Top DEG heatmap
# ══════════════════════════════════════════════════════════════════════════════
# Visualise the top 40 DEGs (by FDR) across all 63 samples.
# Row-scaling (z-score) highlights relative expression differences.
cat("[ STEP 9 ] Heatmap of top DEGs ...\n")

top_genes <- head(res_NASH[res_NASH$sig != "NS", "gene"], 40)
if (length(top_genes) < 5) top_genes <- head(res_NASH$gene, 40)

hm_mat   <- ex_gene[top_genes, , drop = FALSE]
hm_scale <- t(scale(t(hm_mat)))   # z-score per gene

ann_col <- data.frame(condition = pd$condition)
rownames(ann_col) <- colnames(hm_scale)
ann_colors <- list(condition = c(HC = "#2196F3", SS = "#FF9800", NASH = "#F44336"))

png(file.path(plot_dir, "04_heatmap_top40_NASH.png"),
    width = 1400, height = 1600, res = 150)
pheatmap(hm_scale,
         annotation_col  = ann_col,
         annotation_colors = ann_colors,
         cluster_rows    = TRUE,
         cluster_cols    = TRUE,
         show_colnames   = FALSE,
         color           = colorRampPalette(c("#2166AC","white","#B2182B"))(100),
         main            = "Top 40 DEGs – NASH vs HC (row z-score)",
         fontsize_row    = 8)
dev.off()
cat("  Saved: 04_heatmap_top40_NASH.png\n\n")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 10 – GSEA with clusterProfiler (MSigDB Hallmarks)
# ══════════════════════════════════════════════════════════════════════════════
# Gene Set Enrichment Analysis ranks ALL genes by their limma t-statistic,
# then tests whether curated gene sets (pathways) cluster at the top or
# bottom of that ranked list – capturing coordinated shifts even in genes
# that individually fail the FDR threshold.
#
# We use the MSigDB Hallmark collection (50 well-defined biological processes)
# plus KEGG pathways. Gene sets require Entrez IDs.
cat("[ STEP 10 ] Running GSEA ...\n")

# Helper: map gene symbols → Entrez IDs
sym2entrez <- function(symbols) {
  m <- mapIds(org.Hs.eg.db, keys = symbols, column = "ENTREZID",
              keytype = "SYMBOL", multiVals = "first")
  m[!is.na(m)]
}

run_gsea <- function(de_table, label) {
  # Rank metric: signed -log10(p) × sign(logFC) ≈ weighted t-statistic
  # (preserves direction and magnitude, robust to FC inflation)
  ranked <- de_table$t
  names(ranked) <- de_table$gene
  ranked <- sort(ranked, decreasing = TRUE)

  # Convert to Entrez for KEGG / MSigDB
  entrez_map <- sym2entrez(names(ranked))
  ranked_ez  <- ranked[names(ranked) %in% names(entrez_map)]
  names(ranked_ez) <- entrez_map[names(ranked_ez)]
  ranked_ez  <- ranked_ez[!duplicated(names(ranked_ez))]

  results <- list()

  # ── Hallmark gene sets ──────────────────────────────────────────────────────
  hallmark <- msigdbr(species = "Homo sapiens", category = "H") %>%
    dplyr::select(gs_name, entrez_gene) %>%
    dplyr::mutate(entrez_gene = as.character(entrez_gene))

  gsea_h <- GSEA(ranked_ez, TERM2GENE = hallmark,
                 minGSSize = 15, maxGSSize = 500,
                 pvalueCutoff = 0.25, eps = 0,
                 nPermSimple = 10000, verbose = FALSE)
  results$hallmark <- gsea_h

  cat(sprintf("  %s | Hallmarks: %d significant (FDR<0.25)\n",
              label, nrow(as.data.frame(gsea_h))))

  # ── KEGG pathways ───────────────────────────────────────────────────────────
  gsea_k <- gseKEGG(geneList     = ranked_ez,
                    organism     = "hsa",
                    minGSSize    = 15,
                    maxGSSize    = 500,
                    pvalueCutoff = 0.25,
                    eps          = 0,
                    nPermSimple  = 10000,
                    verbose      = FALSE)
  results$kegg <- gsea_k

  cat(sprintf("  %s | KEGG:      %d significant (FDR<0.25)\n",
              label, nrow(as.data.frame(gsea_k))))

  results
}

gsea_NASH <- run_gsea(res_NASH, "NASH_vs_HC")
gsea_SS   <- run_gsea(res_SS,   "SS_vs_HC")

# Save result tables
write.csv(as.data.frame(gsea_NASH$hallmark),
          file.path(res_dir, "GSEA_hallmark_NASH_vs_HC.csv"), row.names = FALSE)
write.csv(as.data.frame(gsea_NASH$kegg),
          file.path(res_dir, "GSEA_KEGG_NASH_vs_HC.csv"),    row.names = FALSE)
write.csv(as.data.frame(gsea_SS$hallmark),
          file.path(res_dir, "GSEA_hallmark_SS_vs_HC.csv"),  row.names = FALSE)
write.csv(as.data.frame(gsea_SS$kegg),
          file.path(res_dir, "GSEA_KEGG_SS_vs_HC.csv"),      row.names = FALSE)
cat("  Saved GSEA result tables.\n\n")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 11 – GSEA dot plots
# ══════════════════════════════════════════════════════════════════════════════
# Each dot = a gene set; x = normalised enrichment score (NES);
# dot size = gene-set size; colour = FDR. Positive NES → enriched in the
# first condition (e.g. NASH); negative → enriched in controls.
cat("[ STEP 11 ] GSEA dot plots ...\n")

make_gsea_dotplot <- function(gsea_obj, title, file, n = 20) {
  df <- as.data.frame(gsea_obj)
  if (nrow(df) == 0) { cat("  (no significant sets – skipping", file, ")\n"); return(invisible(NULL)) }
  df <- df[order(df$NES, decreasing = TRUE), ]
  df <- head(df, n)
  df$ID <- factor(df$ID, levels = df$ID[order(df$NES)])
  df$direction <- ifelse(df$NES > 0, "Up in disease", "Up in HC")

  p <- ggplot(df, aes(NES, ID, size = setSize, colour = p.adjust)) +
    geom_point() +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    scale_colour_gradient(low = "#D73027", high = "#74ADD1",
                          name = "FDR", limits = c(0, 0.25)) +
    scale_size_continuous(name = "Gene set\nsize", range = c(2, 8)) +
    labs(title = title, x = "Normalised Enrichment Score (NES)", y = NULL) +
    theme_bw(base_size = 11) +
    theme(axis.text.y = element_text(size = 8))

  ggsave(file, p,
         width  = 8,
         height = max(4, 0.35 * nrow(df) + 1.5),
         dpi    = 150,
         limitsize = FALSE)
  invisible(p)
}

make_gsea_dotplot(gsea_NASH$hallmark, "GSEA – Hallmarks (NASH vs HC)",
                  file.path(plot_dir, "05_GSEA_hallmark_NASH_vs_HC.png"))
make_gsea_dotplot(gsea_NASH$kegg,     "GSEA – KEGG (NASH vs HC)",
                  file.path(plot_dir, "06_GSEA_KEGG_NASH_vs_HC.png"))
make_gsea_dotplot(gsea_SS$hallmark,   "GSEA – Hallmarks (SS vs HC)",
                  file.path(plot_dir, "07_GSEA_hallmark_SS_vs_HC.png"))
make_gsea_dotplot(gsea_SS$kegg,       "GSEA – KEGG (SS vs HC)",
                  file.path(plot_dir, "08_GSEA_KEGG_SS_vs_HC.png"))
cat("  Saved: 05–08 GSEA dot plots.\n\n")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 12 – Enrichment plots for top pathways (NASH)
# ══════════════════════════════════════════════════════════════════════════════
# The running-sum enrichment plot (Subramanian 2005) shows how a specific gene
# set is distributed along the ranked list, with the peak ES at the "leading edge".
cat("[ STEP 12 ] Enrichment running-sum plots for top NASH pathways ...\n")

plot_top_gsea <- function(gsea_obj, n, prefix) {
  df <- as.data.frame(gsea_obj)
  if (nrow(df) == 0) return(invisible(NULL))
  df <- df[order(df$p.adjust), ]
  top_ids <- head(df$ID, n)
  for (id in top_ids) {
    safe <- gsub("[^A-Za-z0-9_]", "_", id)
    png(file.path(plot_dir, paste0(prefix, "_", safe, ".png")),
        width = 900, height = 550, res = 120)
    print(gseaplot2(gsea_obj, geneSetID = id,
                    title = id, pvalue_table = TRUE))
    dev.off()
  }
  cat(sprintf("  Saved %d running-sum plots (%s)\n", length(top_ids), prefix))
}

plot_top_gsea(gsea_NASH$hallmark, 5, "09_gseaplot_NASH_hallmark")
plot_top_gsea(gsea_NASH$kegg,     5, "10_gseaplot_NASH_KEGG")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 13 – Top DEG tables (console + CSV)
# ══════════════════════════════════════════════════════════════════════════════
cat("\n[ STEP 13 ] Top 20 DEGs – NASH vs HC\n")
top20_NASH <- head(res_NASH[res_NASH$sig != "NS",
                             c("gene","logFC","AveExpr","t","adj.P.Val","sig")], 20)
print(top20_NASH, row.names = FALSE)

cat("\n[ STEP 13 ] Top 20 DEGs – SS vs HC\n")
top20_SS <- head(res_SS[res_SS$sig != "NS",
                         c("gene","logFC","AveExpr","t","adj.P.Val","sig")], 20)
print(top20_SS, row.names = FALSE)

write.csv(top20_NASH, file.path(res_dir, "top20_DEG_NASH_vs_HC.csv"), row.names = FALSE)
write.csv(top20_SS,   file.path(res_dir, "top20_DEG_SS_vs_HC.csv"),   row.names = FALSE)

# ══════════════════════════════════════════════════════════════════════════════
# DONE
# ══════════════════════════════════════════════════════════════════════════════
cat("\n========================================\n")
cat("All outputs saved:\n")
cat("  Plots   →", plot_dir, "\n")
cat("  Results →", res_dir,  "\n")
cat("========================================\n")
