## ============================================================
## Week 4 Day 3 — NAFLD Bulk RNA-seq: GSE130970 (third cohort)
## Dataset: GSE130970 (Salmon/tximport counts, Entrez gene IDs)
## Methodology: protein-coding filter → mean-count filter →
##   DESeq2 (~condition, Normal vs NAFLD) → padj<0.05, |MLE LFC|>1
## Run from: week4-day3-nafld-gse130970-rnaseq/
## ============================================================

suppressPackageStartupMessages({
  library(GEOquery)
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(org.Hs.eg.db)
  library(BiocParallel)
})

set.seed(42)
register(MulticoreParam(4))

if (interactive()) {
  tryCatch({
    proj_dir <- dirname(dirname(rstudioapi::getSourceEditorContext()$path))
    setwd(proj_dir)
  }, error = function(e) invisible(NULL))
}
cat("Working directory:", getwd(), "\n")

dir.create("data",    showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)
dir.create("plots",   showWarnings = FALSE)

# ============================================================
# STEP 1 — Download metadata via GEOquery
# ============================================================
cat("\n=== STEP 1: GSE130970 metadata ===\n")

gse <- getGEO("GSE130970", GSEMatrix = TRUE, destdir = "data/", AnnotGPL = FALSE)
meta <- pData(gse[[1]])

cat("\nAll metadata column names:\n")
print(names(meta))

# Print all characteristics columns to identify condition
char_cols <- grep("characteristics", names(meta), value = TRUE)
cat("\nCharacteristics columns:\n")
for (col in char_cols) {
  cat(sprintf("\n[%s]:\n", col))
  print(sort(unique(as.character(meta[[col]]))))
}

# GSE130970 uses NAFLD activity score (NAS) to stratify disease:
# NAS = 0  → no active NAFLD histology → treat as Normal / control
# NAS > 0  → active NAFLD / NASH at varying severity
# The direct metadata column "nafld activity score:ch1" (or the
# characteristics_ch1.6 column) holds "nafld activity score: N" values.

cat("\n--- Parsing condition from NAFLD activity score ---\n")

# Try the direct metadata column first, then fall back to characteristics
if ("nafld activity score:ch1" %in% names(meta)) {
  nas_raw <- as.character(meta[["nafld activity score:ch1"]])
  cat("Using direct metadata column 'nafld activity score:ch1'\n")
} else {
  # Locate characteristics column that mentions nafld activity score
  nas_col <- NULL
  for (col in char_cols) {
    if (any(grepl("nafld activity score", tolower(as.character(meta[[col]]))))) {
      nas_col <- col; break
    }
  }
  if (is.null(nas_col)) stop("Cannot find NAS column in metadata.")
  nas_raw <- as.character(meta[[nas_col]])
  cat("Using characteristics column:", nas_col, "\n")
}

nas_values <- suppressWarnings(
  as.numeric(sub(".*:\\s*", "", nas_raw))
)
cat("NAS distribution:\n"); print(table(nas_values))

meta$condition <- ifelse(nas_values == 0, "Normal", "NAFLD")
cat("\nFinal condition breakdown:\n"); print(table(meta$condition))

if (sum(meta$condition == "Normal") == 0) {
  stop("No NAS = 0 samples found — check NAS parsing above.")
}

# Also store fibrosis stage for reference
if ("fibrosis stage:ch1" %in% names(meta)) {
  fib_raw <- as.character(meta[["fibrosis stage:ch1"]])
  meta$fibrosis_stage <- suppressWarnings(as.numeric(sub(".*:\\s*", "", fib_raw)))
  cat("Fibrosis stage breakdown:\n"); print(table(meta$fibrosis_stage))
}

# ============================================================
# STEP 2 — Download supplementary count file
# ============================================================
cat("\n=== STEP 2: Download count matrix ===\n")

count_gz <- "data/GSE130970_all_sample_salmon_tximport_counts_entrez_gene_ID.csv.gz"
if (!file.exists(count_gz)) {
  cat("Downloading supplementary count file...\n")
  supp <- getGEOSuppFiles("GSE130970", baseDir = "data/", makeDirectory = FALSE)
  cat("Downloaded files:\n"); print(rownames(supp))
  # Also try direct URL
  if (!file.exists(count_gz)) {
    url <- paste0("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE130nnn/",
                  "GSE130970/suppl/",
                  "GSE130970_all_sample_salmon_tximport_counts_entrez_gene_ID.csv.gz")
    download.file(url, destfile = count_gz, mode = "wb")
  }
}
cat("Count file exists:", file.exists(count_gz), "\n")

# ============================================================
# STEP 3 — Load count matrix
# ============================================================
cat("\n=== STEP 3: Load count matrix ===\n")

counts_raw <- read.csv(gzfile(count_gz), row.names = 1, check.names = FALSE)
cat("Raw count matrix dimensions:", nrow(counts_raw), "genes x", ncol(counts_raw), "samples\n")
cat("First few row names (gene IDs):\n"); print(head(rownames(counts_raw)))
cat("First few column names:\n"); print(head(colnames(counts_raw)))

# Round to integers (tximport estimated counts may be non-integer)
counts_raw <- round(counts_raw)

# ============================================================
# STEP 4 — Align samples with metadata
# ============================================================
cat("\n=== STEP 4: Align samples ===\n")

rownames(meta) <- meta$geo_accession

cat("Sample titles in metadata (first 10):\n")
print(head(meta$title, 10))
cat("Count matrix columns (first 10):\n")
print(head(colnames(counts_raw), 10))

# Strategy 1: direct column name vs geo_accession
common_samples <- intersect(colnames(counts_raw), rownames(meta))
cat("Strategy 1 — direct geo_accession match:", length(common_samples), "\n")

# Strategy 2: numeric prefix of count column matches leading number in title
# count col "440349.1.X_1" → prefix "440349"
# meta title "440349 ..." → leading number "440349"
if (length(common_samples) < ncol(counts_raw) * 0.5) {
  cat("Strategy 2 — title-prefix numeric matching...\n")
  count_col_prefix <- sub("[^0-9].*", "", colnames(counts_raw))

  # Extract all numeric substrings from title and take first
  meta_title_num <- gsub("^[^0-9]*([0-9]+).*", "\\1", meta$title)
  cat("Count col prefixes (first 6):\n"); print(head(count_col_prefix, 6))
  cat("Meta title numbers (first 6):\n");  print(head(meta_title_num, 6))

  match_idx <- match(count_col_prefix, meta_title_num)
  n_matched <- sum(!is.na(match_idx))
  cat("Matched:", n_matched, "of", ncol(counts_raw), "\n")

  if (n_matched >= ncol(counts_raw) * 0.8) {
    keep      <- !is.na(match_idx)
    counts_raw <- counts_raw[, keep, drop = FALSE]
    colnames(counts_raw) <- rownames(meta)[match_idx[keep]]
    common_samples <- colnames(counts_raw)
    cat("Aligned", length(common_samples), "samples via title prefix.\n")
  }
}

# Strategy 3: positional — only if both tables have the same number of samples
if (length(common_samples) < ncol(counts_raw) * 0.5 &&
    ncol(counts_raw) == nrow(meta)) {
  cat("Strategy 3 — positional assignment (same-order assumption).\n")
  cat("WARNING: this assumes count matrix and series matrix have the same sample order.\n")
  colnames(counts_raw) <- rownames(meta)
  common_samples <- rownames(meta)
}

if (length(common_samples) == 0) {
  stop("Cannot match count matrix columns to metadata after all strategies.")
}

cat("Final matched samples:", length(common_samples), "\n")
counts_mat <- as.matrix(counts_raw[, common_samples])
meta_sub   <- meta[common_samples, ]
storage.mode(counts_mat) <- "integer"

cat("Condition check after alignment:\n"); print(table(meta_sub$condition))

# ============================================================
# STEP 5 — PROTEIN-CODING FILTER (Entrez ID → GENETYPE)
# ============================================================
cat("\n=== STEP 5: Protein-coding gene filter ===\n")
cat("Genes before protein-coding filter:", nrow(counts_mat), "\n")

# Map Entrez IDs to gene type and symbol
entrez_keys <- rownames(counts_mat)
gene_annot <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = entrez_keys,
  columns = c("SYMBOL", "GENETYPE"),
  keytype = "ENTREZID"
) %>%
  filter(!duplicated(ENTREZID))

cat("Annotation retrieved for", nrow(gene_annot), "of", length(entrez_keys), "Entrez IDs\n")
cat("Gene types present:\n"); print(table(gene_annot$GENETYPE))

protein_coding_entrez <- gene_annot$ENTREZID[
  !is.na(gene_annot$GENETYPE) & gene_annot$GENETYPE == "protein-coding"
]
cat("Protein-coding Entrez IDs:", length(protein_coding_entrez), "\n")

counts_pc  <- counts_mat[rownames(counts_mat) %in% protein_coding_entrez, ]
n_removed  <- nrow(counts_mat) - nrow(counts_pc)
cat("After protein-coding filter:", nrow(counts_pc),
    "(removed:", n_removed, ")\n")

# ============================================================
# STEP 6 — MEAN-COUNT FILTER (rowMeans >= 10)
# ============================================================
cat("\n=== STEP 6: Mean-count filter (rowMeans >= 10) ===\n")
keep_mean   <- rowMeans(counts_pc) >= 10
counts_filt <- counts_pc[keep_mean, ]
cat("After mean-count filter:", nrow(counts_filt),
    "(removed:", sum(!keep_mean), ")\n")

cat("\n--- Filtering summary ---\n")
cat("  Starting genes:              ", nrow(counts_mat),  "\n")
cat("  After protein-coding filter: ", nrow(counts_pc),   "\n")
cat("  After mean-count filter:     ", nrow(counts_filt), "\n")

# ============================================================
# STEP 7 — DESeq2: Normal vs NAFLD
# ============================================================
cat("\n=== STEP 7: DESeq2 — Normal vs NAFLD ===\n")

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

# MLE (un-shrunk) for threshold filtering
res_mle <- as.data.frame(
  results(dds, contrast = c("condition", "NAFLD", "Normal"), alpha = 0.05)
) %>%
  tibble::rownames_to_column("entrez_id") %>%
  rename(lfc_mle = log2FoldChange, lfcSE_mle = lfcSE,
         stat_mle = stat, pvalue_mle = pvalue, padj_mle = padj)

# apeglm shrinkage for reporting and GSEA ranking
res_ape <- as.data.frame(
  lfcShrink(dds, coef = "condition_NAFLD_vs_Normal", type = "apeglm")
) %>%
  tibble::rownames_to_column("entrez_id") %>%
  dplyr::select(entrez_id, lfc_apeglm = log2FoldChange, lfcSE_apeglm = lfcSE)

cat("\nDESeq2 MLE summary (NAFLD vs Normal):\n")
summary(results(dds, contrast = c("condition", "NAFLD", "Normal"), alpha = 0.05))

# ============================================================
# STEP 8 — Map Entrez IDs to gene symbols
# ============================================================
cat("\n=== STEP 8: Map gene symbols ===\n")

res_df <- res_mle %>%
  left_join(res_ape, by = "entrez_id") %>%
  mutate(
    gene_symbol = mapIds(
      org.Hs.eg.db,
      keys      = entrez_id,
      column    = "SYMBOL",
      keytype   = "ENTREZID",
      multiVals = "first"
    )
  ) %>%
  arrange(padj_mle)

write.csv(res_df, "results/gse130970_results.csv", row.names = FALSE)
cat("Full DESeq2 results saved: results/gse130970_results.csv\n")

# ============================================================
# STEP 9 — Significant genes: padj < 0.05, |MLE LFC| > 1
# ============================================================
cat("\n=== STEP 9: Significant genes (padj<0.05, |MLE LFC|>1) ===\n")

sig_genes <- res_df %>%
  filter(!is.na(padj_mle), !is.na(lfc_mle)) %>%
  filter(padj_mle < 0.05, abs(lfc_mle) > 1) %>%
  arrange(padj_mle)

n_up   <- sum(sig_genes$lfc_mle > 0)
n_down <- sum(sig_genes$lfc_mle < 0)
n_tot  <- nrow(sig_genes)

cat("  Up in NAFLD:  ", n_up, "\n")
cat("  Down in NAFLD:", n_down, "\n")
cat("  Total:        ", n_tot, "\n")

write.csv(sig_genes, "results/significant_genes.csv", row.names = FALSE)
cat("Saved: results/significant_genes.csv\n")

# ============================================================
# STEP 10 — TREM2 / SPP1 / GPNMB check
# ============================================================
cat("\n=== STEP 10: TREM2 / SPP1 / GPNMB ===\n")

for (g in c("TREM2", "SPP1", "GPNMB")) {
  row <- res_df %>% filter(gene_symbol == g)
  if (nrow(row) == 0) { cat(g, ": not found in results\n"); next }
  r <- row[1, ]
  in_sig <- !is.na(r$padj_mle) && r$padj_mle < 0.05 && abs(r$lfc_mle) > 1
  cat(sprintf("%s: MLE log2FC=%+.3f | apeglm=%+.3f | padj=%.2e | sig=%s\n",
              g, r$lfc_mle, r$lfc_apeglm, r$padj_mle, in_sig))
}

# ============================================================
# STEP 11 — PCA plot (VST-transformed)
# ============================================================
cat("\n=== STEP 11: PCA plot ===\n")

vsd      <- vst(dds, blind = TRUE)
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
pca_var  <- round(100 * attr(pca_data, "percentVar"), 1)

p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2, colour = condition)) +
  geom_point(size = 2.5, alpha = 0.85) +
  scale_colour_manual(values = c(Normal = "#2E7D32", NAFLD = "#E65100")) +
  labs(
    title    = "PCA — VST-Normalised Counts",
    subtitle = "GSE130970 | protein-coding genes only | Normal vs NAFLD",
    x        = paste0("PC1 (", pca_var[1], "%)"),
    y        = paste0("PC2 (", pca_var[2], "%)"),
    colour   = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), legend.position = "top")

ggsave("plots/pca.png", p_pca, width = 7, height = 5.5, dpi = 150)
cat("PCA plot saved: plots/pca.png\n")

# ============================================================
# STEP 12 — Volcano plot
# ============================================================
cat("\n=== STEP 12: Volcano plot ===\n")

highlight_genes <- c("TREM2", "SPP1", "GPNMB")

volcano_df <- res_df %>%
  filter(!is.na(padj_mle), !is.na(lfc_mle)) %>%
  mutate(
    label          = ifelse(!is.na(gene_symbol), gene_symbol, entrez_id),
    neg_log10_padj = -log10(padj_mle + 1e-300),
    sig = case_when(
      padj_mle < 0.05 & lfc_mle >  1 ~ "Up in NAFLD",
      padj_mle < 0.05 & lfc_mle < -1 ~ "Down in NAFLD",
      TRUE ~ "NS"
    )
  )

top_sig_labels <- volcano_df %>%
  filter(sig != "NS") %>%
  slice_max(order_by = neg_log10_padj, n = 20)

highlight_df <- volcano_df %>% filter(label %in% highlight_genes)
label_df     <- bind_rows(top_sig_labels, highlight_df) %>%
  distinct(entrez_id, .keep_all = TRUE)

p_volcano <- ggplot(volcano_df, aes(x = lfc_mle, y = neg_log10_padj, colour = sig)) +
  geom_point(alpha = 0.4, size = 0.6) +
  geom_point(data = highlight_df, size = 2.5, shape = 18, colour = "#FFD600") +
  geom_text_repel(
    data         = label_df,
    aes(label    = label),
    size         = 2.5,
    max.overlaps = 25,
    segment.size = 0.3,
    segment.alpha = 0.6,
    colour       = "black"
  ) +
  scale_colour_manual(values = c(
    "Up in NAFLD"   = "#C62828",
    "Down in NAFLD" = "#1565C0",
    "NS"            = "grey70"
  )) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed",
             colour = "black", linewidth = 0.4) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             colour = "black", linewidth = 0.4) +
  labs(
    title    = "Volcano Plot: NAFLD vs Normal Liver",
    subtitle = sprintf(
      "GSE130970 | padj<0.05 & |MLE LFC|>1 | %d up, %d down | TREM2/SPP1/GPNMB = yellow diamonds",
      n_up, n_down),
    x        = "MLE log2 Fold Change (NAFLD / Normal)",
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
cat("DISCOVERY — GSE130970 COMPLETE\n")
cat("============================================================\n")
cat("  Starting genes:              ", nrow(counts_mat),  "\n")
cat("  After protein-coding filter: ", nrow(counts_pc),   "\n")
cat("  After mean-count filter:     ", nrow(counts_filt), "\n")
cat("  Normal samples:              ", sum(col_data$condition == "Normal"), "\n")
cat("  NAFLD samples:               ", sum(col_data$condition == "NAFLD"),  "\n")
cat("  Sig genes (padj<0.05, |MLE LFC|>1):\n")
cat("    Up in NAFLD:  ", n_up,   "\n")
cat("    Down in NAFLD:", n_down, "\n")
cat("    Total:        ", n_tot,  "\n")
