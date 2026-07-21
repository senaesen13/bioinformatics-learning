## ============================================================
## Coexpression Analysis: 139 shared NAFLD genes (Dataset 1 × Dataset 2)
## Expression data: GSE162694 (Dataset 1, n=143, VST-normalised)
## Run from: week6-day2-nafld-coexpression-analysis/
## ============================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(GEOquery)
  library(dplyr)
  library(pheatmap)
  library(ggplot2)
})

set.seed(42)

if (interactive()) {
  tryCatch({
    proj_dir <- dirname(dirname(rstudioapi::getSourceEditorContext()$path))
    setwd(proj_dir)
  }, error = function(e) invisible(NULL))
}
cat("Working directory:", getwd(), "\n")
dir.create("results", showWarnings = FALSE)
dir.create("plots",   showWarnings = FALSE)

d1_dir <- "../week4-day1-nafld-bulk-rnaseq"
d2_dir <- "../week4-day2-nafld-validation-rnaseq"

# ============================================================
# STEP 1 — Compute the 139-gene overlap
# ============================================================
cat("\n=== STEP 1: Compute 139-gene overlap (Day 1 ∩ Day 2) ===\n")

sig_d1 <- read.csv(file.path(d1_dir, "results/significant_genes.csv"),
                   stringsAsFactors = FALSE)
sig_d2 <- read.csv(file.path(d2_dir, "results/significant_genes.csv"),
                   stringsAsFactors = FALSE)

overlap_genes <- intersect(sig_d1$gene_symbol, sig_d2$gene_symbol)
overlap_genes <- overlap_genes[!is.na(overlap_genes) & overlap_genes != ""]
cat("Day 1 sig genes:", nrow(sig_d1), "\n")
cat("Day 2 sig genes:", nrow(sig_d2), "\n")
cat("Overlap genes:  ", length(overlap_genes), "\n")

# ============================================================
# STEP 2 — Load Dataset 1 raw counts + metadata
# ============================================================
cat("\n=== STEP 2: Load GSE162694 counts and metadata ===\n")

gse <- getGEO(filename = file.path(d1_dir, "data/GSE162694_series_matrix.txt.gz"))
meta <- pData(gse)
meta$sample_id    <- sub(".* ", "", meta$title)
meta$fibrosis_raw <- sub("fibrosis stage: ", "", meta$"characteristics_ch1.3")
meta$condition    <- ifelse(meta$fibrosis_raw == "normal liver histology", "Normal", "NAFLD")
rownames(meta)    <- meta$sample_id
cat("Samples:", nrow(meta), "\n")
print(table(meta$condition))

counts_raw <- read.csv(
  file.path(d1_dir, "data/GSE162694/GSE162694_raw_counts.csv"),
  row.names = 1, check.names = FALSE
)
cat("Count matrix:", nrow(counts_raw), "genes x", ncol(counts_raw), "samples\n")

# ============================================================
# STEP 3 — Filter and align (reuse protein-coding list from Day 1)
# ============================================================
cat("\n=== STEP 3: Filter to protein-coding genes, align samples ===\n")

protein_coding_genes <- readRDS(file.path(d1_dir, "results/protein_coding_gene_list.rds"))
counts_pc   <- counts_raw[rownames(counts_raw) %in% protein_coding_genes, ]
keep_mean   <- rowMeans(counts_pc) >= 10
counts_filt <- counts_pc[keep_mean, ]
cat("After protein-coding + mean-count filter:", nrow(counts_filt), "genes\n")

common_samples  <- intersect(colnames(counts_filt), rownames(meta))
counts_mat      <- as.matrix(counts_filt[, common_samples])
meta_sub        <- meta[common_samples, ]
storage.mode(counts_mat) <- "integer"
cat("Samples in common:", length(common_samples), "\n")

# ============================================================
# STEP 4 — VST normalisation
# ============================================================
cat("\n=== STEP 4: VST normalisation ===\n")

col_data <- data.frame(
  condition = factor(meta_sub$condition, levels = c("Normal", "NAFLD")),
  row.names = colnames(counts_mat)
)

dds <- DESeqDataSetFromMatrix(
  countData = counts_mat,
  colData   = col_data,
  design    = ~ condition
)
dds   <- estimateSizeFactors(dds)
vsd   <- vst(dds, blind = TRUE)
vst_all <- assay(vsd)          # all ~10k genes × 143 samples
cat("VST matrix:", nrow(vst_all), "genes x", ncol(vst_all), "samples\n")

# ============================================================
# STEP 5 — Map Ensembl IDs → gene symbols, subset to 139 genes
# ============================================================
cat("\n=== STEP 5: Map Ensembl IDs and subset to overlap genes ===\n")

suppressPackageStartupMessages(library(org.Hs.eg.db))

sym_map <- mapIds(
  org.Hs.eg.db,
  keys    = rownames(vst_all),
  column  = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)
sym_df <- data.frame(
  ensembl_id  = rownames(vst_all),
  gene_symbol = sym_map,
  stringsAsFactors = FALSE
) %>%
  filter(!is.na(gene_symbol), gene_symbol != "")

# Deduplicate: keep one Ensembl ID per gene symbol (lowest row index = first mapped)
sym_df_unique <- sym_df %>%
  filter(gene_symbol %in% overlap_genes) %>%
  group_by(gene_symbol) %>%
  slice(1) %>%
  ungroup()

cat("Overlap genes found in VST matrix:", nrow(sym_df_unique), "\n")

vst_139 <- vst_all[sym_df_unique$ensembl_id, ]
rownames(vst_139) <- sym_df_unique$gene_symbol

cat("VST subset matrix:", nrow(vst_139), "genes x", ncol(vst_139), "samples\n")

# Save VST matrix for reference
write.csv(as.data.frame(vst_139),
          "results/vst_matrix_139genes.csv", row.names = TRUE)
cat("Saved: results/vst_matrix_139genes.csv\n")

# ============================================================
# STEP 6 — Pairwise Pearson correlation matrix
# ============================================================
cat("\n=== STEP 6: Pearson correlation matrix (genes × genes) ===\n")

# Correlation is computed gene × gene (transpose so genes are columns)
cor_mat <- cor(t(vst_139), method = "pearson")
cat("Correlation matrix:", nrow(cor_mat), "x", ncol(cor_mat), "\n")

write.csv(as.data.frame(cor_mat),
          "results/correlation_matrix.csv", row.names = TRUE)
cat("Saved: results/correlation_matrix.csv\n")

# ============================================================
# STEP 7 — Hierarchical clustering and module assignment
# ============================================================
cat("\n=== STEP 7: Hierarchical clustering → modules ===\n")

# Distance = 1 - correlation (so highly correlated genes cluster together)
dist_mat <- as.dist(1 - cor_mat)
hc       <- hclust(dist_mat, method = "average")

# k = 7 selected by silhouette analysis (best avg silhouette = 0.311 across k=2..8)
# Gives two main modules (n=70, n=62) + 5 singleton/small outlier clusters
k <- 7
modules  <- cutree(hc, k = k)
module_df <- data.frame(
  gene_symbol = names(modules),
  module      = paste0("M", modules),
  stringsAsFactors = FALSE
) %>%
  arrange(module, gene_symbol)

cat("\nModule sizes (k =", k, "):\n")
print(table(module_df$module))

# Key gene module report
key_genes <- c("TREM2", "SPP1", "GPNMB", "COL1A1")
cat("\nKey gene modules:\n")
for (g in key_genes) {
  row <- module_df %>% filter(gene_symbol == g)
  if (nrow(row) > 0) {
    cat(sprintf("  %-8s → %s\n", g, row$module[1]))
  } else {
    cat(sprintf("  %-8s → NOT FOUND in overlap\n", g))
  }
}

# Print full module membership
cat("\nFull module membership:\n")
for (m in sort(unique(module_df$module))) {
  genes_in <- module_df$gene_symbol[module_df$module == m]
  cat(sprintf("\n%s (%d genes):\n", m, length(genes_in)))
  cat(paste(sort(genes_in), collapse=", "), "\n")
}

write.csv(module_df, "results/module_assignments.csv", row.names = FALSE)
cat("\nSaved: results/module_assignments.csv\n")

# ============================================================
# STEP 8 — Heatmap (genes ordered by module)
# ============================================================
cat("\n=== STEP 8: Correlation heatmap ===\n")

# Gene annotation: which module each gene belongs to
module_colors <- c(
  M1 = "#1565C0",   # deep blue  (fibrosis/ECM/stellate)
  M2 = "#B71C1C",   # deep red   (macrophage/immune — contains TREM2, SPP1)
  M3 = "#2E7D32",   # deep green
  M4 = "#F57F17",   # amber
  M5 = "#6A1B9A",   # purple
  M6 = "#00695C",   # teal
  M7 = "#4E342E"    # brown
)

annotation_row <- data.frame(
  Module = module_df$module[match(rownames(cor_mat), module_df$gene_symbol)],
  row.names = rownames(cor_mat)
)

ann_colors <- list(Module = module_colors[sort(unique(module_df$module))])

# Order rows/cols by hclust order
gene_order <- hc$labels[hc$order]

pheatmap(
  cor_mat[gene_order, gene_order],
  cluster_rows    = FALSE,
  cluster_cols    = FALSE,
  annotation_row  = annotation_row[gene_order, , drop = FALSE],
  annotation_col  = annotation_row[gene_order, , drop = FALSE],
  annotation_colors = ann_colors,
  color           = colorRampPalette(c("#053061", "#2166AC", "#F7F7F7",
                                       "#D6604D", "#67001F"))(100),
  breaks          = seq(-1, 1, length.out = 101),
  show_rownames   = TRUE,
  show_colnames   = FALSE,
  fontsize_row    = 5,
  border_color    = NA,
  main            = "Coexpression: 139 Shared NAFLD Genes (GSE162694, n=143)",
  filename        = "plots/coexpression_heatmap.png",
  width           = 12,
  height          = 11
)
cat("Saved: plots/coexpression_heatmap.png\n")

# ============================================================
# STEP 9 — Module correlation summary (inter/intra module stats)
# ============================================================
cat("\n=== STEP 9: Module correlation summary ===\n")

for (m in sort(unique(module_df$module))) {
  genes_m <- module_df$gene_symbol[module_df$module == m]
  if (length(genes_m) < 2) next
  sub_cor <- cor_mat[genes_m, genes_m]
  upper   <- sub_cor[upper.tri(sub_cor)]
  cat(sprintf("%s: n=%d genes | mean intra-module r = %.3f | range [%.3f, %.3f]\n",
              m, length(genes_m), mean(upper), min(upper), max(upper)))
}

cat("\n============================================================\n")
cat("COEXPRESSION ANALYSIS COMPLETE\n")
cat("============================================================\n")
