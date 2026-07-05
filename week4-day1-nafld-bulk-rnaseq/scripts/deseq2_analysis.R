## ============================================================
## Week 4 Day 1 — NAFLD Bulk RNA-seq: Corrected Discovery Pipeline
## Dataset: GSE162694 (Suppli et al. 2021, Hepatology)
## Methodology: protein-coding filter → mean-count filter →
##   DESeq2 (~condition, Normal vs NAFLD) → padj<0.01, |LFC|>2
## ============================================================

suppressPackageStartupMessages({
  library(GEOquery)
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(BiocParallel)
})

set.seed(42)
register(MulticoreParam(4))

# Run this script from the week4-day1-nafld-bulk-rnaseq/ project root
# e.g.: cd week4-day1-nafld-bulk-rnaseq && Rscript scripts/deseq2_analysis.R
if (interactive()) {
  tryCatch({
    proj_dir <- dirname(dirname(rstudioapi::getSourceEditorContext()$path))
    setwd(proj_dir)
  }, error = function(e) invisible(NULL))
}
cat("Working directory:", getwd(), "\n")

dir.create("results", showWarnings = FALSE)
dir.create("plots",   showWarnings = FALSE)

# ============================================================
# STEP 1 — Load metadata
# ============================================================
cat("\n=== STEP 1: Load GSE162694 metadata ===\n")

gse <- getGEO("GSE162694", GSEMatrix = TRUE, destdir = "data/")
meta <- pData(gse[[1]])

# sample_id is the second token in the title (e.g. "nash1_F0 548nash1" → "548nash1")
meta$sample_id    <- sub(".* ", "", meta$title)
meta$fibrosis_raw <- sub("fibrosis stage: ", "", meta$"characteristics_ch1.3")
meta$condition    <- ifelse(
  meta$fibrosis_raw == "normal liver histology", "Normal", "NAFLD"
)

cat("Sample counts by condition:\n"); print(table(meta$condition))
cat("Sample counts by fibrosis stage:\n"); print(table(meta$fibrosis_raw))

# ============================================================
# STEP 2 — Load count matrix
# ============================================================
cat("\n=== STEP 2: Load raw count matrix ===\n")

csv_file <- "data/GSE162694/GSE162694_raw_counts.csv"
if (!file.exists(csv_file)) {
  getGEOSuppFiles("GSE162694", baseDir = "data/", makeDirectory = TRUE)
}
counts_raw <- read.csv(csv_file, row.names = 1, check.names = FALSE)
cat("Starting count matrix:", nrow(counts_raw), "genes x", ncol(counts_raw), "samples\n")

# ============================================================
# STEP 3 — Align metadata to count matrix
# ============================================================
cat("\n=== STEP 3: Align samples ===\n")

rownames(meta)  <- meta$sample_id
common_samples  <- intersect(colnames(counts_raw), meta$sample_id)
cat("Samples in common:", length(common_samples), "\n")

counts_mat <- as.matrix(counts_raw[, common_samples])
meta_sub   <- meta[common_samples, ]
storage.mode(counts_mat) <- "integer"

# ============================================================
# STEP 4 — PROTEIN-CODING FILTER (before any count filtering)
# ============================================================
cat("\n=== STEP 4: Protein-coding gene filter ===\n")
cat("Genes before protein-coding filter:", nrow(counts_mat), "\n")

protein_coding_genes <- NULL
biotype_method <- NULL

# Method 1: biomaRt (most current; requires internet)
tryCatch({
  suppressPackageStartupMessages(library(biomaRt))
  cat("Trying biomaRt (useEnsembl)...\n")
  mart <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")
  gene_info <- getBM(
    attributes = c("ensembl_gene_id", "gene_biotype"),
    filters    = "ensembl_gene_id",
    values     = rownames(counts_mat),
    mart       = mart
  )
  protein_coding_genes <- gene_info$ensembl_gene_id[
    gene_info$gene_biotype == "protein_coding"
  ]
  biotype_method <- "biomaRt (Ensembl current)"
  cat("biomaRt query successful.\n")
}, error = function(e) {
  cat("biomaRt failed:", conditionMessage(e), "\n")
})

# Method 2: EnsDb.Hsapiens.v86 (offline fallback, Ensembl 86/GRCh38)
if (is.null(protein_coding_genes)) {
  tryCatch({
    suppressPackageStartupMessages({
      library(EnsDb.Hsapiens.v86)
      library(ensembldb)
    })
    cat("Trying EnsDb.Hsapiens.v86...\n")
    edb_genes <- genes(EnsDb.Hsapiens.v86, return.type = "DataFrame")
    protein_coding_genes <- edb_genes$gene_id[
      edb_genes$gene_biotype == "protein_coding"
    ]
    biotype_method <- "EnsDb.Hsapiens.v86 (Ensembl 86)"
    cat("EnsDb fallback successful.\n")
  }, error = function(e) {
    cat("EnsDb.Hsapiens.v86 failed:", conditionMessage(e), "\n")
  })
}

# Method 3: AnnotationHub (online, cached)
if (is.null(protein_coding_genes)) {
  tryCatch({
    suppressPackageStartupMessages({
      library(AnnotationHub)
      library(ensembldb)
    })
    cat("Trying AnnotationHub...\n")
    ah  <- AnnotationHub()
    edb <- ah[["AH116291"]]  # EnsDb.Hsapiens.v111
    edb_genes <- genes(edb, return.type = "DataFrame")
    protein_coding_genes <- edb_genes$gene_id[
      edb_genes$gene_biotype == "protein_coding"
    ]
    biotype_method <- "AnnotationHub EnsDb.Hsapiens.v111"
    cat("AnnotationHub fallback successful.\n")
  }, error = function(e) {
    cat("AnnotationHub failed:", conditionMessage(e), "\n")
    stop("All biotype annotation methods failed. Cannot continue.")
  })
}

cat("\n--- Biotype annotation method used:", biotype_method, "---\n")
n_pc_in_matrix <- sum(rownames(counts_mat) %in% protein_coding_genes)
n_removed      <- nrow(counts_mat) - n_pc_in_matrix
cat("Protein-coding genes found in count matrix:", n_pc_in_matrix, "\n")
cat("Non-coding / other biotype genes removed:  ", n_removed, "\n")

counts_pc <- counts_mat[rownames(counts_mat) %in% protein_coding_genes, ]
cat("Genes after protein-coding filter:", nrow(counts_pc), "\n")

# Save list for reuse in validation cohort
saveRDS(protein_coding_genes, "results/protein_coding_gene_list.rds")

# ============================================================
# STEP 5 — MEAN-COUNT FILTER (rowMeans >= 10)
# ============================================================
cat("\n=== STEP 5: Mean-count filter (rowMeans >= 10) ===\n")
cat("Genes before mean-count filter:", nrow(counts_pc), "\n")

keep_mean   <- rowMeans(counts_pc) >= 10
counts_filt <- counts_pc[keep_mean, ]
cat("Genes after  mean-count filter:", nrow(counts_filt), "\n")
cat("Genes removed by mean-count filter:", sum(!keep_mean), "\n")

cat("\n--- Filtering summary ---\n")
cat("  Starting genes:                 ", nrow(counts_mat), "\n")
cat("  After protein-coding filter:    ", nrow(counts_pc), "\n")
cat("  After mean-count filter:        ", nrow(counts_filt), "\n")

# ============================================================
# STEP 6 — DESeq2: Normal vs NAFLD
# ============================================================
cat("\n=== STEP 6: DESeq2 — Normal vs NAFLD ===\n")

col_data <- data.frame(
  condition = factor(meta_sub$condition, levels = c("Normal", "NAFLD")),
  row.names = colnames(counts_filt)
)
cat("Samples in DESeq2:\n"); print(table(col_data$condition))

dds <- DESeqDataSetFromMatrix(
  countData = counts_filt,
  colData   = col_data,
  design    = ~ condition
)
dds <- DESeq(dds, parallel = TRUE)

res_shrunk <- lfcShrink(dds,
  coef = "condition_NAFLD_vs_Normal",
  type = "apeglm"
)
cat("\nDESeq2 summary (NAFLD vs Normal, apeglm shrinkage):\n")
summary(res_shrunk)

# ============================================================
# STEP 7 — Map Ensembl IDs to gene symbols
# ============================================================
cat("\n=== STEP 7: Map gene symbols ===\n")
suppressPackageStartupMessages(library(org.Hs.eg.db))

res_df <- as.data.frame(res_shrunk) %>%
  tibble::rownames_to_column("ensembl_id") %>%
  mutate(
    gene_symbol = mapIds(
      org.Hs.eg.db,
      keys      = ensembl_id,
      column    = "SYMBOL",
      keytype   = "ENSEMBL",
      multiVals = "first"
    )
  ) %>%
  arrange(padj)

write.csv(res_df, "results/deseq2_results.csv", row.names = FALSE)
cat("Full DESeq2 results saved: results/deseq2_results.csv\n")

# ============================================================
# STEP 8 — Significant gene list (padj < 0.01, |LFC| > 2)
# ============================================================
cat("\n=== STEP 8: Significant genes (padj < 0.01, |LFC| > 2) ===\n")

sig_genes <- res_df %>%
  filter(!is.na(padj), !is.na(log2FoldChange)) %>%
  filter(padj < 0.01, abs(log2FoldChange) > 2)

n_up   <- sum(sig_genes$log2FoldChange > 0)
n_down <- sum(sig_genes$log2FoldChange < 0)
n_tot  <- nrow(sig_genes)

cat("Significant genes (padj<0.01, |LFC|>2):\n")
cat("  Up in NAFLD:  ", n_up, "\n")
cat("  Down in NAFLD:", n_down, "\n")
cat("  Total:        ", n_tot, "(expected ~1500)\n")

write.csv(sig_genes, "results/significant_genes.csv", row.names = FALSE)
cat("Significant gene list saved: results/significant_genes.csv\n")

# ============================================================
# STEP 9 — TREM2 / SPP1 / GPNMB check
# ============================================================
cat("\n=== STEP 9: TREM2 / SPP1 / GPNMB in significant list ===\n")

marker_genes <- c("TREM2", "SPP1", "GPNMB")
for (g in marker_genes) {
  row <- res_df %>% filter(gene_symbol == g)
  if (nrow(row) == 0) {
    cat(g, ": not found in results\n")
  } else {
    row <- row[1, ]
    in_sig <- !is.na(row$padj) && row$padj < 0.01 && abs(row$log2FoldChange) > 2
    cat(sprintf("%s: log2FC = %.3f | padj = %.2e | significant = %s\n",
                g, row$log2FoldChange, row$padj, in_sig))
  }
}

# ============================================================
# STEP 10 — PCA plot (VST-transformed)
# ============================================================
cat("\n=== STEP 10: PCA plot ===\n")

vsd      <- vst(dds, blind = TRUE)
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
pca_var  <- round(100 * attr(pca_data, "percentVar"), 1)

p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2, colour = condition)) +
  geom_point(size = 2.5, alpha = 0.85) +
  scale_colour_manual(values = c(Normal = "#2E7D32", NAFLD = "#E65100")) +
  labs(
    title    = "PCA — VST-Normalised Counts",
    subtitle = "GSE162694 | protein-coding genes only | Normal vs NAFLD",
    x        = paste0("PC1 (", pca_var[1], "%)"),
    y        = paste0("PC2 (", pca_var[2], "%)"),
    colour   = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), legend.position = "top")

ggsave("plots/pca.png", p_pca, width = 7, height = 5.5, dpi = 150)
cat("PCA plot saved: plots/pca.png\n")

# ============================================================
# STEP 11 — Volcano plot (padj<0.01, |LFC|>2 thresholds)
# ============================================================
cat("\n=== STEP 11: Volcano plot ===\n")

highlight_genes <- c("TREM2", "SPP1", "GPNMB")

volcano_df <- res_df %>%
  filter(!is.na(padj), !is.na(log2FoldChange)) %>%
  mutate(
    label          = ifelse(!is.na(gene_symbol), gene_symbol, ensembl_id),
    neg_log10_padj = -log10(padj + 1e-300),
    sig = case_when(
      padj < 0.01 & log2FoldChange >  2 ~ "Up in NAFLD",
      padj < 0.01 & log2FoldChange < -2 ~ "Down in NAFLD",
      TRUE ~ "NS"
    )
  )

# Label: highlighted markers + top 20 by significance
top_sig_labels <- volcano_df %>%
  filter(sig != "NS") %>%
  slice_max(order_by = neg_log10_padj, n = 20)

highlight_df <- volcano_df %>% filter(label %in% highlight_genes)
label_df     <- bind_rows(top_sig_labels, highlight_df) %>% distinct(ensembl_id, .keep_all = TRUE)

p_volcano <- ggplot(volcano_df, aes(x = log2FoldChange, y = neg_log10_padj, colour = sig)) +
  geom_point(alpha = 0.4, size = 0.6) +
  geom_point(data = highlight_df, size = 2.5, shape = 18, colour = "#FFD600") +
  geom_text_repel(
    data        = label_df,
    aes(label   = label),
    size        = 2.5,
    max.overlaps = 25,
    segment.size = 0.3,
    segment.alpha = 0.6,
    colour      = "black"
  ) +
  scale_colour_manual(values = c(
    "Up in NAFLD"   = "#C62828",
    "Down in NAFLD" = "#1565C0",
    "NS"            = "grey70"
  )) +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed",
             colour = "black", linewidth = 0.4) +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed",
             colour = "black", linewidth = 0.4) +
  labs(
    title    = "Volcano Plot: NAFLD vs Normal Liver",
    subtitle = sprintf("GSE162694 | padj<0.01 & |LFC|>2 | %d up, %d down | TREM2/SPP1/GPNMB = yellow diamonds",
                       n_up, n_down),
    x        = "log2 Fold Change (NAFLD / Normal)",
    y        = expression(-log[10](p[adj])),
    colour   = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), legend.position = "top")

ggsave("plots/volcano.png", p_volcano, width = 9, height = 7, dpi = 150)
cat("Volcano plot saved: plots/volcano.png\n")

# ============================================================
# Summary
# ============================================================
cat("\n============================================================\n")
cat("DISCOVERY COHORT COMPLETE — GSE162694\n")
cat("============================================================\n")
cat("  Biotype annotation method:", biotype_method, "\n")
cat("  Starting genes:          ", nrow(counts_mat),   "\n")
cat("  After protein-coding:    ", nrow(counts_pc),    "\n")
cat("  After mean-count filter: ", nrow(counts_filt),  "\n")
cat("  Sig genes (padj<0.01, |LFC|>2):\n")
cat("    Up in NAFLD:  ", n_up,   "\n")
cat("    Down in NAFLD:", n_down, "\n")
cat("    Total:        ", n_tot,  "\n")
cat("Results: results/deseq2_results.csv\n")
cat("         results/significant_genes.csv\n")
cat("Plots:   plots/pca.png | plots/volcano.png\n")
