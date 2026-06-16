# =============================================================================
# Week 2 — Day 3 — DESeq2: Mouse Myocardial Infarction (KCL Dataset)
# =============================================================================
# Dataset: KCL Module 2 Transcriptomics (Yang Hong, sysmedicine/KCLModule2_2022)
# Kallisto abundances → tximport → DESeq2 (MI vs sham, mouse heart)
# Run this script from: week2-day3-kcl-mouse-deseq2/
# See NOTES.md for experimental design and key differences from the airway analysis.
# =============================================================================


# -----------------------------------------------------------------------------
# 1. INSTALL AND LOAD LIBRARIES
# -----------------------------------------------------------------------------

cran_pkgs <- c("ggplot2", "ggrepel", "dplyr")
missing_cran <- cran_pkgs[!cran_pkgs %in% rownames(installed.packages())]
if (length(missing_cran) > 0) install.packages(missing_cran)

bioc_pkgs <- c("tximport", "DESeq2", "apeglm", "pheatmap")
missing_bioc <- bioc_pkgs[!bioc_pkgs %in% rownames(installed.packages())]
if (length(missing_bioc) > 0) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install(missing_bioc, ask = FALSE)
}

library(tximport)
library(DESeq2)
library(apeglm)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(pheatmap)

dir.create("plots", showWarnings = FALSE)

# Suppress Rplots.pdf — all output is saved explicitly to plots/
if (!interactive()) pdf(NULL)


# -----------------------------------------------------------------------------
# 2. LOAD tx2gene MAPPING
# -----------------------------------------------------------------------------

# proteinCodingGenes maps transcript IDs → gene IDs (Ensembl) for mouse.
# It was pre-fetched from Ensembl BioMart in the original KCL workshop.
# Columns: ensembl_transcript_id, ensembl_gene_id, external_gene_name, description

kcl_dir <- "KCLModule2_2022/Transcriptomics"
load(file.path(kcl_dir, "proteinCodingGenes.Rda"))

head(proteinCodingGenes)
# tximport needs columns 1 (transcript ID) and 2 (gene ID) for tx2gene


# -----------------------------------------------------------------------------
# 3. LOAD KALLISTO ABUNDANCES WITH tximport
# -----------------------------------------------------------------------------

# 8 samples in sorted SRR order map to these biological names:
sample_names <- c(
  "sham_1dMI_1", "sham_1dMI_2",   # healthy, day 1
  "MI_1dMI_1",   "MI_1dMI_2",     # heart attack, day 1
  "sham_3dMI_1", "sham_3dMI_2",   # healthy, day 3
  "MI_3dMI_1",   "MI_3dMI_2"      # heart attack, day 3
)

srr_ids <- paste0("SRR606840", 2:9)   # SRR6068402 … SRR6068409
files   <- file.path(kcl_dir, "abundances", srr_ids, "abundance.h5")
names(files) <- sample_names

# tximport aggregates transcript-level estimates to gene level using tx2gene.
# ignoreTxVersion strips version suffixes (e.g. ENSMUST00000001.5 → ENSMUST00000001)
txi <- tximport(files,
                type            = "kallisto",
                txOut           = FALSE,
                tx2gene         = proteinCodingGenes,
                ignoreTxVersion = TRUE)

mat_count <- txi$counts      # raw counts matrix
mat_tpm   <- txi$abundance   # TPM matrix (for QC/heatmap)

cat("Count matrix dimensions:", dim(mat_count), "\n")   # ~21,906 genes × 8 samples


# -----------------------------------------------------------------------------
# 4. BUILD DESIGN MATRIX
# -----------------------------------------------------------------------------

mat_design <- data.frame(
  row.names = sample_names,
  Treatment = factor(c("sham", "sham", "MI",   "MI",   "sham", "sham", "MI",   "MI"),
                     levels = c("sham", "MI")),   # sham = reference
  Time      = factor(c("1dMI", "1dMI", "1dMI", "1dMI", "3dMI", "3dMI", "3dMI", "3dMI"),
                     levels = c("1dMI", "3dMI"))  # 1dMI = reference
)

mat_design


# -----------------------------------------------------------------------------
# 5. PRE-FILTER LOW-COUNT GENES
# -----------------------------------------------------------------------------

keep      <- rowSums(mat_count) > 10
mat_count <- mat_count[keep, ]
mat_tpm   <- mat_tpm[keep, ]

cat("Genes after filtering:", nrow(mat_count), "\n")


# -----------------------------------------------------------------------------
# 6. CREATE DESeqDataSet
# -----------------------------------------------------------------------------

# Two-factor design: test Treatment effect while controlling for Time.
# Using DESeqDataSetFromMatrix because we already have a count matrix from tximport,
# rather than importing a SummarizedExperiment directly.

dds <- DESeqDataSetFromMatrix(
  countData = round(mat_count),   # tximport returns non-integer estimates; round for DESeq2
  colData   = mat_design,
  design    = ~ Treatment + Time
)

dds


# -----------------------------------------------------------------------------
# 7. RUN DESeq2
# -----------------------------------------------------------------------------

dds <- DESeq(dds, fitType = "parametric")

resultsNames(dds)
# Expected: "Intercept", "Treatment_MI_vs_sham", "Time_3dMI_vs_1dMI"


# -----------------------------------------------------------------------------
# 8. DISPERSION PLOT
# -----------------------------------------------------------------------------

png("plots/dispersion_estimates.png", width = 800, height = 600, res = 120)
plotDispEsts(dds, main = "Dispersion Estimates — Mouse MI")
dev.off()


# -----------------------------------------------------------------------------
# 9. EXTRACT RESULTS: MI vs sham
# -----------------------------------------------------------------------------

# alpha = 0.01: stricter threshold used for independent filtering (this dataset
# has a small dynamic range; 0.01 matches the original KCL workshop cutoff)
res <- results(dds,
               contrast = c("Treatment", "MI", "sham"),
               alpha    = 0.01)

summary(res)
head(res[order(res$padj), ])


# -----------------------------------------------------------------------------
# 10. LFC SHRINKAGE (apeglm)
# -----------------------------------------------------------------------------

# apeglm requires a named coefficient, not a contrast vector.
# "Treatment_MI_vs_sham" is available because sham is set as the reference level.
res_shrunk <- lfcShrink(dds,
                        coef = "Treatment_MI_vs_sham",
                        type = "apeglm")


# -----------------------------------------------------------------------------
# 11. MA PLOT
# -----------------------------------------------------------------------------

png("plots/ma_plot.png", width = 1200, height = 600, res = 120)
par(mfrow = c(1, 2))
plotMA(res,        ylim = c(-6, 6), main = "MA Plot — Raw LFC (MI vs sham)")
plotMA(res_shrunk, ylim = c(-6, 6), main = "MA Plot — Shrunken LFC (MI vs sham)")
par(mfrow = c(1, 1))
dev.off()


# -----------------------------------------------------------------------------
# 12. VOLCANO PLOT
# -----------------------------------------------------------------------------

res_df <- as.data.frame(res_shrunk) %>%
  tibble::rownames_to_column("gene_id") %>%
  filter(!is.na(padj)) %>%
  mutate(
    sig = case_when(
      padj < 0.01 & log2FoldChange >  1 ~ "Up",
      padj < 0.01 & log2FoldChange < -1 ~ "Down",
      TRUE                              ~ "NS"
    ),
    neg_log10_padj = pmin(-log10(padj), 50)
  )

top_genes <- res_df %>%
  filter(sig != "NS") %>%
  slice_min(padj, n = 20)

# Add mouse gene symbols (external_gene_name) for labels
gene_symbols <- proteinCodingGenes[
  match(res_df$gene_id, proteinCodingGenes$ensembl_gene_id),
  "external_gene_name"
]
res_df$symbol <- ifelse(!is.na(gene_symbols) & gene_symbols != "", gene_symbols, res_df$gene_id)

top_genes <- res_df %>%
  filter(sig != "NS") %>%
  slice_min(padj, n = 20)

ggplot(res_df, aes(x = log2FoldChange, y = neg_log10_padj, colour = sig)) +
  geom_point(alpha = 0.5, size = 0.8) +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", colour = "grey40") +
  geom_vline(xintercept = c(-1, 1),     linetype = "dashed", colour = "grey40") +
  geom_text_repel(
    data           = top_genes,
    aes(label      = symbol),
    size           = 2.5,
    max.overlaps   = 20,
    segment.colour = "grey60"
  ) +
  scale_colour_manual(values = c(Up = "#E64B35", Down = "#4DBBD5", NS = "grey70")) +
  labs(
    title    = "Volcano Plot — MI vs Sham (mouse heart)",
    subtitle = "Dashed lines: |LFC| = 1 and padj = 0.01",
    x        = "log2 Fold Change (shrunken, MI / sham)",
    y        = "-log10 Adjusted P-value",
    colour   = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "top")

ggsave("plots/volcano_plot.png", width = 8, height = 6, dpi = 150)


# -----------------------------------------------------------------------------
# 13. VST + PCA
# -----------------------------------------------------------------------------

vst_data <- vst(dds, blind = TRUE)

# PCA coloured by Treatment
pca_trt     <- plotPCA(vst_data, intgroup = "Treatment", returnData = TRUE)
percent_var <- round(100 * attr(pca_trt, "percentVar"))

ggplot(pca_trt, aes(x = PC1, y = PC2, colour = Treatment, label = name)) +
  geom_point(size = 4, alpha = 0.9) +
  geom_text_repel(size = 3, show.legend = FALSE) +
  scale_colour_manual(values = c(sham = "#4DBBD5", MI = "#E64B35")) +
  labs(
    title  = "PCA — VST counts (coloured by Treatment)",
    x      = paste0("PC1: ", percent_var[1], "% variance"),
    y      = paste0("PC2: ", percent_var[2], "% variance")
  ) +
  theme_bw(base_size = 12)

ggsave("plots/pca_treatment.png", width = 7, height = 5, dpi = 150)

# PCA coloured by Time
pca_time <- plotPCA(vst_data, intgroup = "Time", returnData = TRUE)

ggplot(pca_time, aes(x = PC1, y = PC2, colour = Time, label = name)) +
  geom_point(size = 4, alpha = 0.9) +
  geom_text_repel(size = 3, show.legend = FALSE) +
  scale_colour_manual(values = c(`1dMI` = "#F39B7F", `3dMI` = "#8491B4")) +
  labs(
    title  = "PCA — VST counts (coloured by Time)",
    x      = paste0("PC1: ", percent_var[1], "% variance"),
    y      = paste0("PC2: ", percent_var[2], "% variance")
  ) +
  theme_bw(base_size = 12)

ggsave("plots/pca_time.png", width = 7, height = 5, dpi = 150)


# -----------------------------------------------------------------------------
# 14. SAMPLE-TO-SAMPLE DISTANCE HEATMAP
# -----------------------------------------------------------------------------

sample_dists  <- dist(t(assay(vst_data)))
sample_dist_m <- as.matrix(sample_dists)
annotation    <- as.data.frame(colData(dds)[, c("Treatment", "Time")])

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

sig_genes <- res_df %>%
  filter(sig != "NS") %>%
  slice_min(padj, n = 40) %>%   # slice_min returns all rows if fewer than 40 exist
  pull(gene_id)

mat_top    <- assay(vst_data)[sig_genes, ]
mat_scaled <- t(scale(t(mat_top)))

# Replace Ensembl IDs with gene symbols for row labels
rownames(mat_scaled) <- proteinCodingGenes[
  match(rownames(mat_scaled), proteinCodingGenes$ensembl_gene_id),
  "external_gene_name"
]

pheatmap(mat_scaled,
         annotation_col = annotation,
         show_rownames  = TRUE,
         show_colnames  = TRUE,
         main           = "Top DE Genes — Row-scaled VST counts (MI vs sham)",
         fontsize_row   = 8,
         color          = colorRampPalette(c("#4DBBD5", "white", "#E64B35"))(100),
         filename       = "plots/top_genes_heatmap.png",
         width = 8, height = 9)


# -----------------------------------------------------------------------------
# 16. EXPORT RESULTS
# -----------------------------------------------------------------------------

res_export <- res_df %>%
  arrange(padj) %>%
  mutate(across(where(is.numeric), ~ round(.x, 6)))

write.csv(res_export, "deseq2_results_MI_vs_sham.csv", row.names = FALSE)

sig_export <- res_export %>% filter(sig != "NS")
write.csv(sig_export, "deseq2_significant_genes.csv", row.names = FALSE)

cat("\nTotal genes tested:            ", nrow(res_export), "\n")
cat("Significant (FDR<1%, |LFC|>1): ", nrow(sig_export), "\n")
cat("  Upregulated in MI:  ", sum(sig_export$sig == "Up"),   "\n")
cat("  Downregulated in MI:", sum(sig_export$sig == "Down"), "\n")
