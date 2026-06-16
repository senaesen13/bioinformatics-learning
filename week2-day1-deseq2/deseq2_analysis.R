# =============================================================================
# Week 2 — Day 1 — DESeq2 Differential Expression Analysis
# =============================================================================
# Dataset: airway (Bioconductor) — human airway cells, dexamethasone vs untreated
# See NOTES.md for a full explanation of every step.
# =============================================================================


# -----------------------------------------------------------------------------
# 1. INSTALL AND LOAD LIBRARIES
# -----------------------------------------------------------------------------

# Auto-install any missing packages — safe to re-run; skips already-installed ones.

cran_pkgs <- c("ggplot2", "ggrepel", "dplyr")
missing_cran <- cran_pkgs[!cran_pkgs %in% rownames(installed.packages())]
if (length(missing_cran) > 0) {
  message("Installing missing CRAN packages: ", paste(missing_cran, collapse = ", "))
  install.packages(missing_cran)
}

bioc_pkgs <- c("DESeq2", "airway", "apeglm", "pheatmap")
missing_bioc <- bioc_pkgs[!bioc_pkgs %in% rownames(installed.packages())]
if (length(missing_bioc) > 0) {
  message("Installing missing Bioconductor packages: ", paste(missing_bioc, collapse = ", "))
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install(missing_bioc, ask = FALSE)
}

library(DESeq2)
library(airway)
library(apeglm)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(pheatmap)

# Create plots output directory if it doesn't exist
dir.create("plots", showWarnings = FALSE)


# -----------------------------------------------------------------------------
# 2. LOAD THE AIRWAY DATASET
# -----------------------------------------------------------------------------

data("airway")
se <- airway

dim(assay(se))   # 63,677 genes × 8 samples
colData(se)      # key columns: "dex" (trt/untrt) and "cell" (cell line)


# -----------------------------------------------------------------------------
# 3. CREATE DESeqDataSet + SET REFERENCE LEVEL
# -----------------------------------------------------------------------------

# design = ~ cell + dex: control for cell-line differences, test the dex effect
dds <- DESeqDataSet(se, design = ~ cell + dex)

# Make "untrt" the reference so LFC = log2(treated / untreated)
dds$dex <- relevel(dds$dex, ref = "untrt")


# -----------------------------------------------------------------------------
# 4. PRE-FILTER LOW-COUNT GENES
# -----------------------------------------------------------------------------

# Remove genes with fewer than 10 total reads across all samples
keep <- rowSums(counts(dds)) >= 10
dds  <- dds[keep, ]

dim(dds)


# -----------------------------------------------------------------------------
# 5. SIZE FACTOR ESTIMATION (NORMALISATION)
# -----------------------------------------------------------------------------

# Median-of-ratios method: corrects for differences in sequencing depth
dds <- estimateSizeFactors(dds)

sizeFactors(dds)

normalised_counts <- counts(dds, normalized = TRUE)
head(normalised_counts)


# -----------------------------------------------------------------------------
# 6. DISPERSION ESTIMATION
# -----------------------------------------------------------------------------

# Empirical Bayes shrinkage: shrinks noisy per-gene dispersions toward a
# fitted trend, borrowing strength across genes
dds <- estimateDispersions(dds)

png("plots/dispersion_estimates.png", width = 800, height = 600, res = 120)
plotDispEsts(dds, main = "Dispersion Estimates\nblack=raw, red=trend, blue=shrunken")
dev.off()


# -----------------------------------------------------------------------------
# 7. WALD TEST
# -----------------------------------------------------------------------------

# Fits a GLM per gene and tests whether the dex LFC is significantly ≠ 0
# (If estimateSizeFactors / estimateDispersions were already run, DESeq()
#  skips those steps and just runs the GLM + Wald test.)
dds <- DESeq(dds)

resultsNames(dds)


# -----------------------------------------------------------------------------
# 8. EXTRACT RESULTS
# -----------------------------------------------------------------------------

res <- results(dds,
               name  = "dex_trt_vs_untrt",
               alpha = 0.05)

summary(res)

head(res[order(res$padj), ])


# -----------------------------------------------------------------------------
# 9. LFC SHRINKAGE (apeglm)
# -----------------------------------------------------------------------------

# Shrinks noisy fold changes for low-count genes toward zero.
# P-values are unchanged — only the LFC estimates are corrected.
# Use shrunken LFCs for all plots and gene rankings.
res_shrunk <- lfcShrink(dds,
                        coef = "dex_trt_vs_untrt",
                        type = "apeglm")

head(cbind(raw_lfc      = res$log2FoldChange,
           shrunken_lfc = res_shrunk$log2FoldChange,
           padj         = res$padj))


# -----------------------------------------------------------------------------
# 10. MA PLOT
# -----------------------------------------------------------------------------

png("plots/ma_plot.png", width = 1200, height = 600, res = 120)
par(mfrow = c(1, 2))
plotMA(res,        ylim = c(-5, 5), main = "MA Plot — Raw LFC")
plotMA(res_shrunk, ylim = c(-5, 5), main = "MA Plot — Shrunken LFC")
par(mfrow = c(1, 1))
dev.off()


# -----------------------------------------------------------------------------
# 11. VOLCANO PLOT
# -----------------------------------------------------------------------------

res_df <- as.data.frame(res_shrunk) %>%
  tibble::rownames_to_column("gene_id") %>%
  filter(!is.na(padj)) %>%
  mutate(
    sig = case_when(
      padj < 0.05 & log2FoldChange >  1 ~ "Up",
      padj < 0.05 & log2FoldChange < -1 ~ "Down",
      TRUE                              ~ "NS"
    ),
    neg_log10_padj = pmin(-log10(padj), 50)
  )

top_genes <- res_df %>%
  filter(sig != "NS") %>%
  slice_min(padj, n = 20)

ggplot(res_df, aes(x = log2FoldChange, y = neg_log10_padj, colour = sig)) +
  geom_point(alpha = 0.5, size = 0.8) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey40") +
  geom_vline(xintercept = c(-1, 1),     linetype = "dashed", colour = "grey40") +
  geom_text_repel(
    data           = top_genes,
    aes(label      = gene_id),
    size           = 2.5,
    max.overlaps   = 20,
    segment.colour = "grey60"
  ) +
  scale_colour_manual(values = c(Up = "#E64B35", Down = "#4DBBD5", NS = "grey70")) +
  labs(
    title    = "Volcano Plot — Dexamethasone vs Untreated",
    subtitle = "Dashed lines: |LFC| = 1 and padj = 0.05",
    x        = "log2 Fold Change (shrunken)",
    y        = "-log10 Adjusted P-value",
    colour   = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "top")

ggsave("plots/volcano_plot.png", width = 8, height = 6, dpi = 150)


# -----------------------------------------------------------------------------
# 12. VARIANCE STABILISING TRANSFORMATION (VST)
# -----------------------------------------------------------------------------

# Equalises variance across the expression range before PCA and heatmaps.
# blind = TRUE: ignore the design so we see raw data structure.
vst_data <- vst(dds, blind = TRUE)


# -----------------------------------------------------------------------------
# 13. PCA
# -----------------------------------------------------------------------------

pca_data    <- plotPCA(vst_data, intgroup = c("dex", "cell"), returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))

ggplot(pca_data, aes(x = PC1, y = PC2, colour = dex, shape = cell)) +
  geom_point(size = 4, alpha = 0.9) +
  geom_text_repel(aes(label = name), size = 3, show.legend = FALSE) +
  scale_colour_manual(values = c(trt = "#E64B35", untrt = "#4DBBD5")) +
  labs(
    title  = "PCA — VST-transformed counts",
    x      = paste0("PC1: ", percent_var[1], "% variance"),
    y      = paste0("PC2: ", percent_var[2], "% variance"),
    colour = "Treatment",
    shape  = "Cell line"
  ) +
  theme_bw(base_size = 12)

ggsave("plots/pca_plot.png", width = 7, height = 5, dpi = 150)


# -----------------------------------------------------------------------------
# 14. SAMPLE-TO-SAMPLE DISTANCE HEATMAP
# -----------------------------------------------------------------------------

sample_dists  <- dist(t(assay(vst_data)))
sample_dist_m <- as.matrix(sample_dists)
annotation    <- as.data.frame(colData(dds)[, c("cell", "dex")])

pheatmap(sample_dist_m,
         clustering_distance_rows = sample_dists,
         clustering_distance_cols = sample_dists,
         annotation_col           = annotation,
         main                     = "Sample-to-Sample Distance",
         color                    = colorRampPalette(c("#2166AC", "white"))(50),
         filename                 = "plots/sample_distance_heatmap.png",
         width = 7, height = 6)


# -----------------------------------------------------------------------------
# 15. TOP DE GENES HEATMAP
# -----------------------------------------------------------------------------

top40 <- res_df %>%
  filter(sig != "NS") %>%
  slice_min(padj, n = 40) %>%
  pull(gene_id)

# Z-score each gene so all rows are on the same scale
mat_top    <- assay(vst_data)[top40, ]
mat_scaled <- t(scale(t(mat_top)))

pheatmap(mat_scaled,
         annotation_col = annotation,
         show_rownames  = TRUE,
         show_colnames  = TRUE,
         main           = "Top 40 DE Genes — Row-scaled VST counts",
         fontsize_row   = 7,
         color          = colorRampPalette(c("#4DBBD5", "white", "#E64B35"))(100),
         filename       = "plots/top40_genes_heatmap.png",
         width = 8, height = 9)


# -----------------------------------------------------------------------------
# 16. EXPORT RESULTS
# -----------------------------------------------------------------------------

res_export <- res_df %>%
  arrange(padj) %>%
  mutate(across(where(is.numeric), ~ round(.x, 6)))

write.csv(res_export, "deseq2_results.csv", row.names = FALSE)

sig_genes <- res_export %>% filter(sig != "NS")
write.csv(sig_genes, "deseq2_significant_genes.csv", row.names = FALSE)

cat("Total genes tested:       ", nrow(res_export), "\n")
cat("Significant (FDR<5%, |LFC|>1):", nrow(sig_genes), "\n")
cat("  Upregulated:  ", sum(sig_genes$sig == "Up"), "\n")
cat("  Downregulated:", sum(sig_genes$sig == "Down"), "\n")
