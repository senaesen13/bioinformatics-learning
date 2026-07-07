## GSE166047 — NW vs OBF subcutaneous adipose tissue RNA-seq
## Rey/Messa et al. — obesity adipose transcriptomics practice analysis
## Week 5 Day 1

suppressPackageStartupMessages({
  library(GEOquery)
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(msigdbr)
})

BASE_DIR <- "/Users/senaesen/Desktop/bioinfo-learning/week5-day1-obesity-adipose-rnaseq"
RESULTS  <- file.path(BASE_DIR, "results")
PLOTS    <- file.path(BASE_DIR, "plots")

## ── Step 1: Metadata ─────────────────────────────────────────────────────────
cat("\n=== STEP 1: Metadata ===\n")
meta_cache <- file.path(RESULTS, "metadata_cache.csv")

if (!file.exists(meta_cache)) {
  cat("Downloading GEO metadata...\n")
  gse      <- getGEO("GSE166047", GSEMatrix = TRUE, destdir = RESULTS)
  metadata <- pData(gse[[1]])
  metadata$group <- gsub("[0-9]+_S[0-9]+$", "", metadata$title)
  write.csv(metadata, meta_cache, row.names = FALSE)
  cat("Metadata cached.\n")
} else {
  cat("Loading metadata from cache.\n")
  metadata <- read.csv(meta_cache, stringsAsFactors = FALSE)
}

cat("Group table:\n"); print(table(metadata$group))

## ── Step 2: Download raw counts from RAW.tar ─────────────────────────────────
cat("\n=== STEP 2: Raw count files ===\n")
tar_url  <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE166nnn/GSE166047/suppl/GSE166047_RAW.tar"
tar_path <- file.path(RESULTS, "GSE166047_RAW.tar")

if (!file.exists(tar_path)) {
  cat("Downloading RAW.tar...\n")
  download.file(tar_url, tar_path, mode = "wb", quiet = FALSE)
} else {
  cat("RAW.tar already cached.\n")
}

target_samples <- c(
  "GSM5060703_NW1_S1.dat.gz",  "GSM5060704_NW2_S2.dat.gz",
  "GSM5060705_NW3_S3.dat.gz",  "GSM5060706_NW4_S4.dat.gz",
  "GSM5060707_NW5_S5.dat.gz",
  "GSM5060708_OBF1_S6.dat.gz", "GSM5060709_OBF2_S7.dat.gz",
  "GSM5060710_OBF3_S8.dat.gz", "GSM5060711_OBF4_S9.dat.gz",
  "GSM5060712_OBF5_S10.dat.gz"
)

extract_dir <- file.path(RESULTS, "raw_files")
dir.create(extract_dir, showWarnings = FALSE)
missing <- target_samples[!file.exists(file.path(extract_dir, target_samples))]
if (length(missing) > 0) {
  cat("Extracting", length(missing), "sample files...\n")
  untar(tar_path, files = missing, exdir = extract_dir)
}
cat("Sample files ready:", length(list.files(extract_dir)), "\n")

## ── Step 3: Build count matrix ───────────────────────────────────────────────
cat("\n=== STEP 3: Building count matrix ===\n")

read_dat <- function(filepath) {
  df <- read.table(gzfile(filepath), header = TRUE, sep = "\t",
                   stringsAsFactors = FALSE)
  df <- df[, c("gene_ID", "gene_name", "frags_locus")]
  sample_name <- sub("\\.dat\\.gz$", "", basename(filepath))
  sample_name <- sub("^GSM[0-9]+_", "", sample_name)
  colnames(df)[3] <- sample_name
  df
}

dat_list  <- lapply(file.path(extract_dir, target_samples), read_dat)
counts_df <- Reduce(function(a, b) merge(a, b, by = c("gene_ID", "gene_name")), dat_list)
rownames(counts_df) <- counts_df$gene_ID
gene_names  <- setNames(counts_df$gene_name, counts_df$gene_ID)
count_mat   <- as.matrix(counts_df[, -(1:2)])
cat("Count matrix:", dim(count_mat)[1], "genes x", dim(count_mat)[2], "samples\n")
cat("Samples:", paste(colnames(count_mat), collapse = ", "), "\n")

## ── Step 4: Protein-coding gene filter ───────────────────────────────────────
cat("\n=== STEP 4: Protein-coding gene filter ===\n")
cat("Genes before filter:", nrow(count_mat), "\n")

# mapIds is robust for large queries; multiVals="first" handles 1:many mappings
genetype_vec <- suppressMessages(
  mapIds(org.Hs.eg.db,
         keys      = rownames(count_mat),
         column    = "GENETYPE",
         keytype   = "ENSEMBL",
         multiVals = "first")
)
keep_pc      <- !is.na(genetype_vec) & genetype_vec == "protein-coding"
count_mat_pc <- count_mat[keep_pc, ]
cat("Genes after protein-coding filter:", nrow(count_mat_pc), "\n")

## ── Step 5: Mean-count filter ─────────────────────────────────────────────────
cat("\n=== STEP 5: Mean-count filter (rowMeans >= 10) ===\n")
cat("Genes before mean-count filter:", nrow(count_mat_pc), "\n")
count_mat_filt <- count_mat_pc[rowMeans(count_mat_pc) >= 10, ]
cat("Genes after mean-count filter:", nrow(count_mat_filt), "\n")

## ── Step 6: DESeq2 ───────────────────────────────────────────────────────────
cat("\n=== STEP 6: DESeq2 ===\n")

col_data <- data.frame(
  sample    = colnames(count_mat_filt),
  condition = ifelse(grepl("^NW", colnames(count_mat_filt)), "NW", "OBF"),
  row.names = colnames(count_mat_filt)
)
col_data$condition <- factor(col_data$condition, levels = c("NW", "OBF"))
cat("Sample conditions:\n"); print(col_data)

dds <- DESeqDataSetFromMatrix(
  countData = count_mat_filt,
  colData   = col_data,
  design    = ~ condition
)
dds <- DESeq(dds)
cat("Total genes tested:", nrow(dds), "\n")

res_raw    <- results(dds, contrast = c("condition", "OBF", "NW"), alpha = 0.01)
res_shrunk <- lfcShrink(dds, coef = "condition_OBF_vs_NW", type = "apeglm")

cat("\nResults summary (apeglm shrunk, alpha=0.01):\n")
summary(res_shrunk, alpha = 0.01)

res_df <- as.data.frame(res_shrunk)
res_df$gene_id   <- rownames(res_df)
res_df$gene_name <- gene_names[res_df$gene_id]
res_df$MLE_log2FC <- as.data.frame(res_raw)$log2FoldChange[match(res_df$gene_id, rownames(res_raw))]
write.csv(res_df, file.path(RESULTS, "deseq2_results_all.csv"), row.names = FALSE)

## ── Step 7: Significant genes ────────────────────────────────────────────────
cat("\n=== STEP 7: Significant genes (padj < 0.01, |MLE log2FC| > 2) ===\n")
sig <- res_df %>%
  filter(!is.na(padj), padj < 0.01, abs(MLE_log2FC) > 2) %>%
  arrange(padj)

cat("Upregulated (OBF > NW):", sum(sig$MLE_log2FC > 0), "\n")
cat("Downregulated (OBF < NW):", sum(sig$MLE_log2FC < 0), "\n")
cat("Total significant:", nrow(sig), "\n")

cat("\nTop 10 by padj:\n")
print(sig %>% select(gene_name, gene_id, MLE_log2FC, log2FoldChange, padj) %>% head(10))

write.csv(sig, file.path(RESULTS, "deseq2_significant_genes.csv"), row.names = FALSE)

## ── Step 8: Plots ────────────────────────────────────────────────────────────
cat("\n=== STEP 8: Plots ===\n")

# PCA (VST-normalized)
vst_data <- vst(dds, blind = TRUE)
pca_data <- plotPCA(vst_data, intgroup = "condition", returnData = TRUE)
pct_var  <- round(100 * attr(pca_data, "percentVar"), 1)

pca_plot <- ggplot(pca_data, aes(PC1, PC2, color = condition, label = name)) +
  geom_point(size = 4, alpha = 0.85) +
  geom_text_repel(size = 3, max.overlaps = 20) +
  scale_color_manual(values = c("NW" = "#4477AA", "OBF" = "#EE6677")) +
  labs(
    title = paste0("PCA - GSE166047 NW vs OBF Subcutaneous Adipose"),
    x     = paste0("PC1 (", pct_var[1], "%)"),
    y     = paste0("PC2 (", pct_var[2], "%)")
  ) +
  theme_bw(base_size = 13)

ggsave(file.path(PLOTS, "pca_NW_vs_OBF.pdf"), pca_plot, width = 7, height = 5)
ggsave(file.path(PLOTS, "pca_NW_vs_OBF.png"), pca_plot, width = 7, height = 5, dpi = 150)
cat("PCA plot saved.\n")

# Volcano plot
vol_df <- res_df %>%
  filter(!is.na(padj)) %>%
  mutate(
    significance = case_when(
      padj < 0.01 & MLE_log2FC >  2 ~ "Up",
      padj < 0.01 & MLE_log2FC < -2 ~ "Down",
      TRUE                           ~ "NS"
    )
  )

top_labels <- vol_df %>% filter(significance != "NS") %>% arrange(padj) %>% head(5)

volcano_plot <- ggplot(vol_df, aes(MLE_log2FC, -log10(padj), color = significance)) +
  geom_point(alpha = 0.4, size = 1.2) +
  geom_point(data = filter(vol_df, significance != "NS"), alpha = 0.8, size = 1.8) +
  geom_text_repel(data = top_labels, aes(label = gene_name),
                  size = 3.5, max.overlaps = 20, color = "black") +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c("Up" = "#CC3333", "Down" = "#3366CC", "NS" = "grey70")) +
  labs(
    title = "Volcano - GSE166047 NW vs OBF Subcutaneous Adipose",
    x     = "MLE log2 Fold Change (OBF / NW)",
    y     = "-log10(padj)"
  ) +
  theme_bw(base_size = 13)

ggsave(file.path(PLOTS, "volcano_NW_vs_OBF.pdf"), volcano_plot, width = 8, height = 6)
ggsave(file.path(PLOTS, "volcano_NW_vs_OBF.png"), volcano_plot, width = 8, height = 6, dpi = 150)
cat("Volcano plot saved.\n")

## ── Step 9: GSEA ─────────────────────────────────────────────────────────────
cat("\n=== STEP 9: GSEA ===\n")

# Ranked gene list by apeglm-shrunk log2FC, using gene symbols
ranked_df <- res_df %>%
  filter(!is.na(log2FoldChange), !is.na(gene_name), gene_name != "") %>%
  arrange(desc(log2FoldChange))
gene_list <- setNames(ranked_df$log2FoldChange, ranked_df$gene_name)
gene_list <- gene_list[!duplicated(names(gene_list))]
cat("Genes in ranked list:", length(gene_list), "\n")

run_gsea <- function(pathway_df, label) {
  cat("Running GSEA:", label, "...\n")
  set.seed(42)
  gsea_res <- GSEA(
    geneList     = gene_list,
    TERM2GENE    = pathway_df[, c("gs_name", "gene_symbol")],
    pvalueCutoff = 1,
    minGSSize    = 15,
    maxGSSize    = 500,
    eps          = 0,
    verbose      = FALSE
  )
  gsea_df <- as.data.frame(gsea_res)
  write.csv(gsea_df, file.path(RESULTS, paste0("gsea_", label, ".csv")), row.names = FALSE)

  top5_up   <- gsea_df %>% filter(NES > 0) %>% arrange(p.adjust) %>% head(5)
  top5_down <- gsea_df %>% filter(NES < 0) %>% arrange(p.adjust) %>% head(5)

  cat("\n--- Top 5 Activated (", label, ") ---\n")
  print(top5_up[, c("ID", "NES", "p.adjust")])
  cat("\n--- Top 5 Suppressed (", label, ") ---\n")
  print(top5_down[, c("ID", "NES", "p.adjust")])

  list(result = gsea_res, df = gsea_df, up = top5_up, down = top5_down)
}

# msigdbr v10+ uses collection/subcollection arguments
kegg_sets     <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_LEGACY")
hallmark_sets <- msigdbr(species = "Homo sapiens", collection = "H")

kegg_gsea     <- run_gsea(kegg_sets,     "KEGG")
hallmark_gsea <- run_gsea(hallmark_sets, "Hallmark")

cat("\n=== PIPELINE COMPLETE ===\n")
cat("Results saved to:", RESULTS, "\n")
cat("Plots saved to:", PLOTS, "\n")
