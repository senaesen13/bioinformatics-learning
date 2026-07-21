## ============================================================
## Week 4 Day 1 — NAFLD Bulk RNA-seq: DESeq2 + GSEA  (_ak version)
## Dataset: GSE162694  (Suppli et al. 2021, Hepatology)
## Human liver biopsy | Normal (n=31) vs NAFLD F0-F4 (n=112) | 143 samples
##
## Changes vs Sena's original:
##   1. Download also fetches GPL21290 platform table (gene biotype)
##   2. Low-expression filter: CPM>1 in ≥10% of samples (n≥15) — replaces
##      the arbitrary >10 counts in ≥3 samples threshold
##   3. Protein-coding filter: retains only genes annotated as
##      "protein_coding" in the Ensembl biotype column from GPL21290 /
##      biomaRt fallback; lncRNA, pseudogene, etc. are excluded
##   4. GSEA ranking metric: signed –log10(pvalue) × sign(LFC) instead of
##      raw LFC — preserves direction while giving more dynamic range for
##      genes with tiny but consistent effects; biologically more sensitive
##   5. GSEA pvalueCutoff set to 0.05 (FDR) instead of 0.25 — 0.25 is
##      the clusterProfiler default but is far too lenient for publication
##   6. GSEA nPermSimple = 10000 (was default 1000) for stable p-values
##   7. Additional gene sets: Reactome (C2:CP:REACTOME) from MSigDB
##   8. PCA coloured by both condition and fibrosis stage (two panels)
##   9. Sample-level QC: flagged outliers via Mahalanobis distance on PCA
##  10. Biotype filter stats and filter log printed at each step
##
## NOTE (updated after Sena's Week 4 Day 1 rebuild, see deseq2_analysis.R /
## gsea_analysis.R): she independently adopted protein-coding-first filtering
## and a stricter significance threshold (padj<0.01, |LFC|>2) with her own
## MLE-vs-apeglm rationale — see NOTES.md. Points 4-9 above (signed -log10(p)
## GSEA ranking, Reactome gene sets, fibrosis-stage PCA panel, Mahalanobis
## outlier flagging) remain unique to this _ak script and are kept here as a
## standalone comparison, not merged into her pipeline.
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
  library(patchwork)
  library(edgeR)        # for cpm()
})

set.seed(42)
register(MulticoreParam(4))

# Run this script from the week4-day1-nafld-bulk-rnaseq/ project root
# e.g.: cd week4-day1-nafld-bulk-rnaseq && Rscript scripts/nafld_deseq2_gsea_ak.R
if (interactive()) {
  tryCatch({
    proj_dir <- dirname(dirname(rstudioapi::getSourceEditorContext()$path))
    setwd(proj_dir)
  }, error = function(e) invisible(NULL))
}
cat("Working directory:", getwd(), "\n")

OUTDIR_PLOTS   <- "plots_ak"
OUTDIR_RESULTS <- "results_ak"
dir.create(OUTDIR_PLOTS,   showWarnings = FALSE, recursive = TRUE)
dir.create(OUTDIR_RESULTS, showWarnings = FALSE, recursive = TRUE)

## ============================================================
## STEP 1 — Download GSE metadata + GPL platform annotation
## ============================================================
cat("\n=== STEP 1: Download GSE162694 metadata + GPL21290 platform table ===\n")

gse <- getGEO("GSE162694", GSEMatrix = TRUE, destdir = "data/")
meta <- pData(gse[[1]])
cat("Loaded metadata:", nrow(meta), "samples\n")

# Sample ID: second token of title field
# e.g. "nash1_F0 548nash1" -> "548nash1"
meta$sample_id <- sub(".* ", "", meta$title)

# Fibrosis stage
meta$fibrosis_raw <- sub("fibrosis stage: ", "", meta$"characteristics_ch1.3")

# Binary condition
meta$condition <- ifelse(
  meta$fibrosis_raw == "normal liver histology", "Normal", "NAFLD"
)

# Ordered fibrosis factor (for PCA colouring)
fibrosis_levels <- c("normal liver histology", "F0", "F1", "F2", "F3", "F4")
meta$fibrosis_ordered <- factor(meta$fibrosis_raw,
                                levels = fibrosis_levels,
                                ordered = TRUE)

cat("\nSample counts by condition:\n")
print(table(meta$condition))
cat("\nSample counts by fibrosis stage:\n")
print(table(meta$fibrosis_raw))

# Download GPL21290 platform annotation (gene biotype column)
cat("\nDownloading GPL21290 platform annotation...\n")
gpl <- getGEO("GPL21290", destdir = "data/")
gpl_table <- Table(gpl)
cat("GPL21290 table columns:", paste(colnames(gpl_table), collapse = ", "), "\n")
cat("GPL21290 rows:", nrow(gpl_table), "\n")

# Identify the Ensembl ID column and biotype column
# GPL21290 typically has columns: ID (ENSG...), gene_biotype or similar
ensembl_col  <- intersect(c("ID", "ensembl_gene_id", "Ensembl_Gene_ID"),
                           colnames(gpl_table))[1]
biotype_col  <- intersect(c("gene_biotype", "Gene_Biotype", "gene_type", "GENE_BIOTYPE"),
                           colnames(gpl_table))[1]
symbol_col   <- intersect(c("gene_name", "Gene_Symbol", "GENE_SYMBOL", "Symbol"),
                           colnames(gpl_table))[1]

if (is.na(ensembl_col) || is.na(biotype_col)) {
  stop("GPL21290 table is missing an expected Ensembl ID or biotype column ",
       "(found: ", paste(colnames(gpl_table), collapse = ", "),
       "). The biomaRt fallback below only triggers on empty coding_ids, ",
       "not on missing columns, so this must fail fast instead of silently ",
       "producing NA-filled annotations.")
}

cat("Using columns — Ensembl:", ensembl_col,
    "| Biotype:", biotype_col,
    "| Symbol:", symbol_col, "\n")

# Build gene annotation lookup
gene_annot <- gpl_table %>%
  dplyr::select(ensembl_id = all_of(ensembl_col),
                biotype    = all_of(biotype_col),
                symbol_gpl = all_of(symbol_col)) %>%
  dplyr::distinct(ensembl_id, .keep_all = TRUE)

cat("Biotype distribution (top 10):\n")
print(sort(table(gene_annot$biotype), decreasing = TRUE)[1:10])

## ============================================================
## STEP 2 — Load supplementary count matrix
## ============================================================
cat("\n=== STEP 2: Load count matrix ===\n")

csv_file <- "data/GSE162694/GSE162694_raw_counts.csv"
if (!file.exists(csv_file)) {
  getGEOSuppFiles("GSE162694", baseDir = "data/", makeDirectory = TRUE)
  # gunzip if still compressed
  gz_file <- paste0(csv_file, ".gz")
  if (file.exists(gz_file)) R.utils::gunzip(gz_file, remove = FALSE)
}

counts_raw <- read.csv(csv_file, row.names = 1, check.names = FALSE)
cat("Raw count matrix:", nrow(counts_raw), "genes x", ncol(counts_raw), "samples\n")
cat("First 4 gene IDs:", head(rownames(counts_raw), 4), "\n")

## ============================================================
## STEP 3 — Align samples
## ============================================================
cat("\n=== STEP 3: Align metadata to count matrix ===\n")

rownames(meta)  <- meta$sample_id
common_samples  <- intersect(colnames(counts_raw), meta$sample_id)
cat("Samples in common:", length(common_samples), "\n")

counts_mat <- as.matrix(counts_raw[, common_samples])
meta_sub   <- meta[common_samples, ]
storage.mode(counts_mat) <- "integer"

## ============================================================
## STEP 4 — Protein-coding filter (from GPL21290)
## ============================================================
cat("\n=== STEP 4: Protein-coding gene filter ===\n")

cat("Genes before biotype filter:", nrow(counts_mat), "\n")

# Genes in the count matrix
mat_genes <- rownames(counts_mat)

# Annotate with biotype
coding_ids <- gene_annot %>%
  filter(ensembl_id %in% mat_genes, biotype == "protein_coding") %>%
  pull(ensembl_id)

# Fallback: if GPL had no biotype info, use biomaRt
if (length(coding_ids) == 0) {
  cat("WARNING: No biotype found in GPL table — falling back to biomaRt\n")
  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    stop("Install biomaRt: BiocManager::install('biomaRt')")
  }
  mart <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  bm <- biomaRt::getBM(
    attributes = c("ensembl_gene_id", "gene_biotype", "hgnc_symbol"),
    filters    = "ensembl_gene_id",
    values     = mat_genes,
    mart       = mart
  )
  write.csv(bm, file.path(OUTDIR_RESULTS, "biomart_gene_annotation.csv"),
            row.names = FALSE)
  coding_ids <- bm %>%
    filter(gene_biotype == "protein_coding") %>%
    pull(ensembl_gene_id) %>%
    unique()
}

cat("Protein-coding genes in matrix:", length(coding_ids), "\n")
counts_coding <- counts_mat[rownames(counts_mat) %in% coding_ids, ]
cat("After protein-coding filter:", nrow(counts_coding), "genes\n")
cat("Removed:", nrow(counts_mat) - nrow(counts_coding),
    "non-coding / pseudogene / lncRNA / etc.\n")

## ============================================================
## STEP 5 — Low-expression filter (CPM-based)
## ============================================================
cat("\n=== STEP 5: Low-expression filter (CPM > 1 in >= 10% samples) ===\n")

# edgeR::cpm normalises per sample library size — much more principled than
# a raw count threshold, which is biased by sequencing depth.
# Threshold: CPM > 1 in at least 15 samples (10% of 143).
# This corresponds roughly to "a gene must be detectably expressed in the
# smallest biological group (n=8 for F3), so groups of ≥10% of the cohort."
cpm_mat       <- edgeR::cpm(counts_coding)
min_samples   <- max(round(0.10 * ncol(counts_coding)), 3)
keep_cpm      <- rowSums(cpm_mat > 1) >= min_samples
counts_filt   <- counts_coding[keep_cpm, ]

cat("CPM threshold: >1 CPM in >=", min_samples, "samples\n")
cat("Before filter:", nrow(counts_coding), "protein-coding genes\n")
cat("After filter: ", nrow(counts_filt), "genes retained\n")
cat("Removed (low expression):", nrow(counts_coding) - nrow(counts_filt), "\n")

## ============================================================
## STEP 6 — DESeq2: Normal vs NAFLD
## ============================================================
cat("\n=== STEP 6: DESeq2 — Normal vs NAFLD ===\n")

col_data <- data.frame(
  condition = factor(meta_sub$condition, levels = c("Normal", "NAFLD")),
  fibrosis  = factor(meta_sub$fibrosis_raw),
  row.names = colnames(counts_filt)
)

dds <- DESeqDataSetFromMatrix(
  countData = counts_filt,
  colData   = col_data,
  design    = ~ condition
)

dds <- DESeq(dds, parallel = TRUE)

res        <- results(dds, contrast = c("condition", "NAFLD", "Normal"), alpha = 0.05)
res_shrunk <- lfcShrink(dds, coef = "condition_NAFLD_vs_Normal", type = "apeglm")

cat("\nDESeq2 summary (NAFLD vs Normal, protein-coding only):\n")
summary(res_shrunk)

# Build result data frame with gene symbols
# Merge with GPL annotation (which has gene symbols) + org.Hs.eg.db for Entrez
res_df <- as.data.frame(res_shrunk) %>%
  tibble::rownames_to_column("ensembl_id") %>%
  left_join(gene_annot %>% dplyr::select(ensembl_id, symbol_gpl, biotype),
            by = "ensembl_id") %>%
  mutate(
    gene_symbol_orgdb = mapIds(org.Hs.eg.db,
                               keys      = ensembl_id,
                               column    = "SYMBOL",
                               keytype   = "ENSEMBL",
                               multiVals = "first"),
    # Prefer GPL symbol; fall back to org.Hs.eg.db
    gene_symbol = ifelse(!is.na(symbol_gpl) & symbol_gpl != "",
                         symbol_gpl, gene_symbol_orgdb),
    entrez = mapIds(org.Hs.eg.db,
                    keys      = ensembl_id,
                    column    = "ENTREZID",
                    keytype   = "ENSEMBL",
                    multiVals = "first")
  ) %>%
  arrange(padj)

write.csv(res_df,
          file.path(OUTDIR_RESULTS, "deseq2_nafld_vs_normal.csv"),
          row.names = FALSE)
cat("Full DESeq2 results saved.\n")

## ============================================================
## STEP 7 — Volcano plot
## ============================================================
cat("\n=== STEP 7: Volcano plot ===\n")

volcano_df <- res_df %>%
  filter(!is.na(padj), !is.na(log2FoldChange)) %>%
  mutate(
    label          = ifelse(!is.na(gene_symbol), gene_symbol, ensembl_id),
    sig            = case_when(
      padj < 0.05 & log2FoldChange >  1 ~ "Up in NAFLD",
      padj < 0.05 & log2FoldChange < -1 ~ "Down in NAFLD",
      TRUE ~ "NS"
    ),
    neg_log10_padj = -log10(padj)
  )

top_labels <- volcano_df %>%
  filter(sig != "NS") %>%
  slice_max(order_by = neg_log10_padj, n = 30)

cat("Up in NAFLD   (padj<0.05, LFC>1):", sum(volcano_df$sig == "Up in NAFLD"), "\n")
cat("Down in NAFLD (padj<0.05, LFC<-1):", sum(volcano_df$sig == "Down in NAFLD"), "\n")

p_volcano <- ggplot(volcano_df,
                    aes(x = log2FoldChange, y = neg_log10_padj, colour = sig)) +
  geom_point(alpha = 0.4, size = 0.65) +
  geom_text_repel(data = top_labels, aes(label = label),
                  size = 2.3, max.overlaps = 25,
                  segment.size = 0.25, segment.alpha = 0.6) +
  scale_colour_manual(values = c("Up in NAFLD"   = "#C62828",
                                  "Down in NAFLD" = "#1565C0",
                                  "NS"            = "grey72")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed",
             colour = "black", linewidth = 0.35) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             colour = "black", linewidth = 0.35) +
  labs(title    = "Volcano Plot: NAFLD vs Normal Liver (protein-coding genes)",
       subtitle = "GSE162694 | 143 samples | DESeq2 + apeglm | FDR < 0.05 & |LFC| > 1",
       x        = "log₂ Fold Change (NAFLD / Normal)",
       y        = expression(-log[10](p[adj])),
       colour   = NULL) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), legend.position = "top")

ggsave(file.path(OUTDIR_PLOTS, "volcano_nafld_vs_normal.pdf"),
       p_volcano, width = 8, height = 6.5)
ggsave(file.path(OUTDIR_PLOTS, "volcano_nafld_vs_normal.png"),
       p_volcano, width = 8, height = 6.5, dpi = 150)
cat("Volcano plot saved.\n")

## ============================================================
## STEP 8 — Top 20 DEGs
## ============================================================
cat("\n=== STEP 8: Top 20 DEGs ===\n")

top20 <- res_df %>%
  filter(!is.na(padj), padj < 0.05, !is.na(gene_symbol)) %>%
  slice_max(order_by = abs(log2FoldChange), n = 20) %>%
  dplyr::select(gene_symbol, ensembl_id, biotype, baseMean,
                log2FoldChange, lfcSE, pvalue, padj)

print(as.data.frame(top20))
write.csv(top20, file.path(OUTDIR_RESULTS, "top20_degs.csv"), row.names = FALSE)

p_top20 <- top20 %>%
  mutate(
    direction = ifelse(log2FoldChange > 0, "Up in NAFLD", "Down in NAFLD"),
    gene_symbol = reorder(gene_symbol, log2FoldChange)
  ) %>%
  ggplot(aes(x = log2FoldChange, y = gene_symbol, fill = direction)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = c("Up in NAFLD"   = "#C62828",
                                "Down in NAFLD" = "#1565C0")) +
  labs(title    = "Top 20 DEGs: NAFLD vs Normal (protein-coding)",
       subtitle = "Ranked by absolute log₂ Fold Change | padj < 0.05",
       x = "log₂ Fold Change", y = NULL, fill = NULL) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"), legend.position = "top")

ggsave(file.path(OUTDIR_PLOTS, "top20_degs_barplot.pdf"), p_top20, width = 7, height = 6)
ggsave(file.path(OUTDIR_PLOTS, "top20_degs_barplot.png"), p_top20, width = 7, height = 6, dpi = 150)
cat("Top 20 DEG bar chart saved.\n")

## ============================================================
## STEP 9 — PCA (two panels: condition + fibrosis stage)
## ============================================================
cat("\n=== STEP 9: PCA (VST) — condition + fibrosis stage panels ===\n")

vsd      <- vst(dds, blind = TRUE)

# Condition panel
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
pca_var  <- round(100 * attr(pca_data, "percentVar"), 1)

pca_data$fibrosis <- meta_sub[rownames(pca_data), "fibrosis_raw"]

# Outlier detection via Mahalanobis distance on PC1+PC2
pca_xy   <- as.matrix(pca_data[, c("PC1", "PC2")])
mah_dist <- mahalanobis(pca_xy, colMeans(pca_xy), cov(pca_xy))
pca_data$outlier <- mah_dist > qchisq(0.975, df = 2)
n_out <- sum(pca_data$outlier)
cat("PCA outliers (Mahalanobis p < 0.025):", n_out, "\n")
if (n_out > 0) {
  cat("Outlier samples:\n")
  print(rownames(pca_data)[pca_data$outlier])
}

p_pca_cond <- ggplot(pca_data,
                     aes(x = PC1, y = PC2, colour = condition,
                         shape = outlier)) +
  geom_point(size = 2.5, alpha = 0.85) +
  scale_colour_manual(values = c(Normal = "#2E7D32", NAFLD = "#E65100")) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 4),
                     labels = c("Normal", "Outlier"),
                     name   = NULL) +
  labs(title    = "PCA — VST counts (by condition)",
       subtitle  = "GSE162694 | protein-coding genes",
       x = paste0("PC1 (", pca_var[1], "%)"),
       y = paste0("PC2 (", pca_var[2], "%)"),
       colour = NULL) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), legend.position = "top")

fibrosis_colours <- c(
  "normal liver histology" = "#2E7D32",
  "F0" = "#AED581", "F1" = "#FFF176",
  "F2" = "#FFB300", "F3" = "#E64A19", "F4" = "#880E4F"
)

p_pca_fib <- ggplot(pca_data,
                    aes(x = PC1, y = PC2, colour = fibrosis)) +
  geom_point(size = 2.5, alpha = 0.85) +
  scale_colour_manual(values = fibrosis_colours, name = "Fibrosis stage") +
  labs(title    = "PCA — VST counts (by fibrosis stage)",
       subtitle  = "GSE162694 | protein-coding genes",
       x = paste0("PC1 (", pca_var[1], "%)"),
       y = paste0("PC2 (", pca_var[2], "%)")) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), legend.position = "right")

p_pca_combined <- p_pca_cond + p_pca_fib
ggsave(file.path(OUTDIR_PLOTS, "pca_combined.pdf"),
       p_pca_combined, width = 14, height = 5.5)
ggsave(file.path(OUTDIR_PLOTS, "pca_combined.png"),
       p_pca_combined, width = 14, height = 5.5, dpi = 150)
cat("PCA plots saved.\n")

## ============================================================
## STEP 10 — GSEA (Hallmark + KEGG + Reactome)
##
## Ranking metric: signed –log10(pvalue) × sign(LFC)
##   Rationale: raw LFC (Sena's approach) collapses genes with similar
##   fold-change but very different statistical confidence.
##   Signed –log10(p) gives higher weight to consistently expressed genes
##   across many samples.  Standard in many NAFLD/clinical RNA-seq papers.
##   Both metrics are scientifically defensible; we use the more sensitive one.
##
## pvalueCutoff = 0.05 (FDR-corrected q-value) — standard for publication.
##   Sena used 0.25 which is the clusterProfiler default but inflates results.
##
## nPermSimple = 10000 — 10× more permutations for stable p-values
##   (especially important with small gene sets).
## ============================================================
cat("\n=== STEP 10: GSEA ===\n")

# Build signed ranking metric
ranked_df <- res_df %>%
  filter(!is.na(pvalue), !is.na(log2FoldChange)) %>%
  mutate(
    stat = sign(log2FoldChange) * (-log10(pvalue))
  )

## — Helper to build a deduplicated named ranked vector —
make_ranked_sym <- function(df) {
  df <- df %>%
    filter(!is.na(gene_symbol)) %>%
    arrange(desc(stat))
  v <- setNames(df$stat, df$gene_symbol)
  v[!duplicated(names(v))]
}

make_ranked_entrez <- function(df) {
  df <- df %>%
    filter(!is.na(entrez)) %>%
    arrange(desc(stat))
  v <- setNames(df$stat, df$entrez)
  v[!duplicated(names(v))]
}

ranked_sym    <- make_ranked_sym(ranked_df)
ranked_entrez <- make_ranked_entrez(ranked_df)

cat("Genes in symbol-ranked list:", length(ranked_sym), "\n")
cat("Genes in Entrez-ranked list:", length(ranked_entrez), "\n")

## — 10a: Hallmark —
cat("\nRunning Hallmark GSEA...\n")
hallmark_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, gene_symbol)

gsea_hallmark <- GSEA(
  geneList     = ranked_sym,
  TERM2GENE    = hallmark_t2g,
  minGSSize    = 15,
  maxGSSize    = 500,
  pvalueCutoff = 0.05,
  eps          = 0,
  nPermSimple  = 10000,
  verbose      = FALSE,
  seed         = 42
)
cat("Significant Hallmark sets (FDR<0.05):", nrow(gsea_hallmark), "\n")
write.csv(as.data.frame(gsea_hallmark),
          file.path(OUTDIR_RESULTS, "gsea_hallmark_results.csv"),
          row.names = FALSE)

if (nrow(gsea_hallmark) > 0) {
  p_hm_dot <- dotplot(gsea_hallmark, showCategory = 20,
                       split = ".sign",
                       title = "Hallmark GSEA — NAFLD vs Normal (FDR<0.05)") +
    facet_grid(. ~ .sign) +
    theme_bw(base_size = 10) +
    theme(plot.title  = element_text(face = "bold"),
          axis.text.y = element_text(size = 7),
          strip.text  = element_text(face = "bold"))
  ggsave(file.path(OUTDIR_PLOTS, "gsea_hallmark_dotplot.pdf"),
         p_hm_dot, width = 14, height = 8)
  ggsave(file.path(OUTDIR_PLOTS, "gsea_hallmark_dotplot.png"),
         p_hm_dot, width = 14, height = 8, dpi = 150)

  if (nrow(gsea_hallmark) >= 5) {
    p_hm_ridge <- ridgeplot(gsea_hallmark, showCategory = 20) +
      labs(title = "Hallmark GSEA — Ridge Plot") +
      theme_bw(base_size = 10) +
      theme(plot.title  = element_text(face = "bold"),
            axis.text.y = element_text(size = 7))
    ggsave(file.path(OUTDIR_PLOTS, "gsea_hallmark_ridgeplot.pdf"),
           p_hm_ridge, width = 12, height = 8)
    ggsave(file.path(OUTDIR_PLOTS, "gsea_hallmark_ridgeplot.png"),
           p_hm_ridge, width = 12, height = 8, dpi = 150)
  }
}

## — 10b: KEGG —
cat("\nRunning KEGG GSEA...\n")
gsea_kegg <- gseKEGG(
  geneList     = ranked_entrez,
  organism     = "hsa",
  minGSSize    = 15,
  maxGSSize    = 500,
  pvalueCutoff = 0.05,
  eps          = 0,
  nPermSimple  = 10000,
  verbose      = FALSE,
  seed         = 42
)
cat("Significant KEGG pathways (FDR<0.05):", nrow(gsea_kegg), "\n")
write.csv(as.data.frame(gsea_kegg),
          file.path(OUTDIR_RESULTS, "gsea_kegg_results.csv"),
          row.names = FALSE)

if (nrow(gsea_kegg) > 0) {
  p_kegg_dot <- dotplot(gsea_kegg, showCategory = 20,
                         split = ".sign",
                         title = "KEGG GSEA — NAFLD vs Normal (FDR<0.05)") +
    facet_grid(. ~ .sign) +
    theme_bw(base_size = 10) +
    theme(plot.title  = element_text(face = "bold"),
          axis.text.y = element_text(size = 8),
          strip.text  = element_text(face = "bold"))
  ggsave(file.path(OUTDIR_PLOTS, "gsea_kegg_dotplot.pdf"),
         p_kegg_dot, width = 13, height = 8)
  ggsave(file.path(OUTDIR_PLOTS, "gsea_kegg_dotplot.png"),
         p_kegg_dot, width = 13, height = 8, dpi = 150)

  if (nrow(gsea_kegg) >= 5) {
    p_kegg_ridge <- ridgeplot(gsea_kegg, showCategory = 15) +
      labs(title = "KEGG GSEA — Ridge Plot") +
      theme_bw(base_size = 10) +
      theme(plot.title  = element_text(face = "bold"),
            axis.text.y = element_text(size = 8))
    ggsave(file.path(OUTDIR_PLOTS, "gsea_kegg_ridgeplot.pdf"),
           p_kegg_ridge, width = 11, height = 8)
    ggsave(file.path(OUTDIR_PLOTS, "gsea_kegg_ridgeplot.png"),
           p_kegg_ridge, width = 11, height = 8, dpi = 150)
  }
}

## — 10c: Reactome (C2:CP:REACTOME from MSigDB) —
cat("\nRunning Reactome GSEA...\n")
reactome_t2g <- msigdbr(species = "Homo sapiens",
                         category = "C2", subcategory = "CP:REACTOME") %>%
  dplyr::select(gs_name, gene_symbol)

gsea_reactome <- GSEA(
  geneList     = ranked_sym,
  TERM2GENE    = reactome_t2g,
  minGSSize    = 15,
  maxGSSize    = 500,
  pvalueCutoff = 0.05,
  eps          = 0,
  nPermSimple  = 10000,
  verbose      = FALSE,
  seed         = 42
)
cat("Significant Reactome pathways (FDR<0.05):", nrow(gsea_reactome), "\n")
write.csv(as.data.frame(gsea_reactome),
          file.path(OUTDIR_RESULTS, "gsea_reactome_results.csv"),
          row.names = FALSE)

if (nrow(gsea_reactome) > 0) {
  p_react_dot <- dotplot(gsea_reactome, showCategory = 20,
                          split = ".sign",
                          title = "Reactome GSEA — NAFLD vs Normal (FDR<0.05)") +
    facet_grid(. ~ .sign) +
    theme_bw(base_size = 10) +
    theme(plot.title  = element_text(face = "bold"),
          axis.text.y = element_text(size = 7),
          strip.text  = element_text(face = "bold"))
  ggsave(file.path(OUTDIR_PLOTS, "gsea_reactome_dotplot.pdf"),
         p_react_dot, width = 16, height = 9)
  ggsave(file.path(OUTDIR_PLOTS, "gsea_reactome_dotplot.png"),
         p_react_dot, width = 16, height = 9, dpi = 150)

  if (nrow(gsea_reactome) >= 5) {
    p_react_ridge <- ridgeplot(gsea_reactome, showCategory = 20) +
      labs(title = "Reactome GSEA — Ridge Plot") +
      theme_bw(base_size = 10) +
      theme(plot.title  = element_text(face = "bold"),
            axis.text.y = element_text(size = 7))
    ggsave(file.path(OUTDIR_PLOTS, "gsea_reactome_ridgeplot.pdf"),
           p_react_ridge, width = 14, height = 9)
    ggsave(file.path(OUTDIR_PLOTS, "gsea_reactome_ridgeplot.png"),
           p_react_ridge, width = 14, height = 9, dpi = 150)
  }
}

## ============================================================
## STEP 11 — Filter summary report
## ============================================================
cat("\n=== STEP 11: Filter summary ===\n")

filter_summary <- data.frame(
  Step  = c("Raw count matrix",
            "After protein-coding filter",
            "After CPM>1 in >=10% samples",
            "DEGs padj<0.05",
            "DEGs padj<0.05 & |LFC|>1 Up",
            "DEGs padj<0.05 & |LFC|>1 Down"),
  Genes = c(nrow(counts_mat),
             nrow(counts_coding),
             nrow(counts_filt),
             sum(!is.na(res_df$padj) & res_df$padj < 0.05, na.rm = TRUE),
             sum(volcano_df$sig == "Up in NAFLD"),
             sum(volcano_df$sig == "Down in NAFLD"))
)
print(filter_summary)
write.csv(filter_summary,
          file.path(OUTDIR_RESULTS, "filter_summary.csv"),
          row.names = FALSE)

cat("\n============================================================\n")
cat("All analyses complete.\n")
cat("Plots  ->", OUTDIR_PLOTS, "\n")
cat("Results->", OUTDIR_RESULTS, "\n")
