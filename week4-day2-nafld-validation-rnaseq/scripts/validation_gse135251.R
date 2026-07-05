## ============================================================
## Week 4 Day 2 — NAFLD RNA-seq Validation Pipeline
## Validation cohort: GSE135251 (Govaere et al. 2020, Sci Transl Med)
## 206 NAFLD + 10 normal liver biopsies | RNA-seq | GRCh38
## Cross-cohort comparison with GSE162694 discovery results
##
## Run from: week4-day2-nafld-validation-rnaseq/
##   e.g.  cd week4-day2-nafld-validation-rnaseq && Rscript scripts/validation_gse135251.R
##
## Requires discovery results in:
##   ../week4-day1-nafld-bulk-rnaseq/results/significant_genes.csv
##   ../week4-day1-nafld-bulk-rnaseq/results/deseq2_results.csv
## ============================================================

suppressPackageStartupMessages({
  library(GEOquery)
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(org.Hs.eg.db)
  library(dplyr)
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

disc_dir <- "../week4-day1-nafld-bulk-rnaseq"
dir.create("data/GSE135251", recursive = TRUE, showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)
dir.create("plots",   showWarnings = FALSE)

# ============================================================
# STEP 1 — Load metadata and define conditions
# ============================================================
cat("\n=== STEP 1: GSE135251 metadata ===\n")

gse_val <- getGEO("GSE135251", GSEMatrix = TRUE, destdir = "data/GSE135251/",
                  AnnotGPL = FALSE)
meta_val <- pData(gse_val[[1]])

cat("All metadata column names:\n"); print(names(meta_val))
cat("\nCharacteristics unique values:\n")
for (col in grep("characteristics", names(meta_val), value = TRUE)) {
  cat(sprintf("\n[%s]:\n", col))
  print(sort(unique(as.character(meta_val[[col]]))))
}

# Binary condition: "disease: Control" → Normal, "disease: NAFLD" → NAFLD
meta_val$condition <- ifelse(
  sub("disease: ", "", meta_val$"disease:ch1") == "Control", "Normal", "NAFLD"
)
cat("\nCondition breakdown:\n"); print(table(meta_val$condition))
cat("Group in paper breakdown:\n")
print(table(sub("group in paper: ", "", meta_val$"group in paper:ch1")))

# ============================================================
# STEP 2 — Download supplementary count files (if needed)
# ============================================================
cat("\n=== STEP 2: Supplementary count files ===\n")

supp_dir   <- "data/GSE135251"
count_files <- list.files(supp_dir, pattern = "\\.counts\\.txt\\.gz$", full.names = TRUE)

if (length(count_files) == 0) {
  tar_path <- file.path(supp_dir, "GSE135251_RAW.tar")
  if (!file.exists(tar_path)) {
    cat("Downloading GSE135251_RAW.tar (~44 MB)...\n")
    download.file(
      "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE135nnn/GSE135251/suppl/GSE135251_RAW.tar",
      destfile = tar_path, mode = "wb", quiet = FALSE
    )
  }
  cat("Extracting tar archive...\n")
  untar(tar_path, exdir = supp_dir)
  count_files <- list.files(supp_dir, pattern = "\\.counts\\.txt\\.gz$",
                              full.names = TRUE)
}
cat("Count files found:", length(count_files), "\n")

# ============================================================
# STEP 3 — Build count matrix from per-sample HTSeq files
# ============================================================
cat("\n=== STEP 3: Build count matrix ===\n")

gsm_ids <- sub("_(.*)\\.counts\\.txt\\.gz$", "", basename(count_files))
first   <- read.table(gzfile(count_files[1]), header = FALSE, sep = "\t",
                       stringsAsFactors = FALSE)
gene_ids <- first[!grepl("^__", first[, 1]), 1]
cat("Genes per file:", length(gene_ids), "\n")
cat("Reading", length(count_files), "files...\n")

count_mat <- matrix(0L, nrow = length(gene_ids), ncol = length(count_files),
                    dimnames = list(gene_ids, gsm_ids))
for (i in seq_along(count_files)) {
  d <- read.table(gzfile(count_files[i]), header = FALSE, sep = "\t",
                   stringsAsFactors = FALSE)
  d <- d[!grepl("^__", d[, 1]), ]
  count_mat[d[, 1], i] <- as.integer(d[, 2])
}
cat("Raw count matrix built:", nrow(count_mat), "x", ncol(count_mat), "\n")

# ============================================================
# STEP 4 — Align metadata to count matrix
# ============================================================
cat("\n=== STEP 4: Align samples ===\n")

rownames(meta_val) <- meta_val$geo_accession
common_val         <- intersect(colnames(count_mat), meta_val$geo_accession)
cat("Samples in common:", length(common_val), "\n")

counts_val_mat <- count_mat[, common_val]
meta_val_sub   <- meta_val[common_val, ]
storage.mode(counts_val_mat) <- "integer"
cat("Condition check:\n"); print(table(meta_val_sub$condition))

# ============================================================
# STEP 5 — PROTEIN-CODING FILTER
# Uses the same gene list derived in the discovery cohort (biomaRt)
# ============================================================
cat("\n=== STEP 5: Protein-coding filter ===\n")
cat("Genes before filter:", nrow(counts_val_mat), "\n")

pc_genes      <- readRDS("results/protein_coding_gene_list.rds")
counts_val_pc <- counts_val_mat[rownames(counts_val_mat) %in% pc_genes, ]
cat("After protein-coding filter:", nrow(counts_val_pc),
    "(removed:", nrow(counts_val_mat) - nrow(counts_val_pc), ")\n")

# ============================================================
# STEP 6 — MEAN-COUNT FILTER (rowMeans >= 10)
# ============================================================
cat("\n=== STEP 6: Mean-count filter (rowMeans >= 10) ===\n")
keep_val     <- rowMeans(counts_val_pc) >= 10
counts_val_f <- counts_val_pc[keep_val, ]
cat("After mean-count filter:", nrow(counts_val_f),
    "(removed:", sum(!keep_val), ")\n")

cat("\n--- Filtering summary GSE135251 ---\n")
cat("  Starting genes:              ", nrow(counts_val_mat), "\n")
cat("  After protein-coding filter: ", nrow(counts_val_pc),  "\n")
cat("  After mean-count filter:     ", nrow(counts_val_f),   "\n")

# ============================================================
# STEP 7 — DESeq2: Normal vs NAFLD
# ============================================================
cat("\n=== STEP 7: DESeq2 — Normal vs NAFLD ===\n")

col_val <- data.frame(
  condition = factor(meta_val_sub$condition, levels = c("Normal", "NAFLD")),
  row.names = colnames(counts_val_f)
)
cat("Samples:\n"); print(table(col_val$condition))

dds_val <- DESeqDataSetFromMatrix(countData = counts_val_f, colData = col_val,
                                   design = ~ condition)
dds_val <- DESeq(dds_val, parallel = TRUE)

# MLE (un-shrunk) for filtering; apeglm for reporting
res_val_mle <- as.data.frame(
  results(dds_val, contrast = c("condition", "NAFLD", "Normal"), alpha = 0.01)
) %>%
  tibble::rownames_to_column("ensembl_id") %>%
  rename(lfc_mle = log2FoldChange, lfcSE_mle = lfcSE,
         stat_mle = stat, pvalue_mle = pvalue, padj_mle = padj)

res_val_ape <- as.data.frame(
  lfcShrink(dds_val, coef = "condition_NAFLD_vs_Normal", type = "apeglm")
) %>%
  tibble::rownames_to_column("ensembl_id") %>%
  select(ensembl_id, lfc_apeglm = log2FoldChange, lfcSE_apeglm = lfcSE)

res_val_df <- res_val_mle %>%
  left_join(res_val_ape, by = "ensembl_id") %>%
  mutate(
    gene_symbol = mapIds(org.Hs.eg.db, keys = ensembl_id,
                         column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
  ) %>%
  arrange(padj_mle)

write.csv(res_val_df, "results/gse135251_results.csv", row.names = FALSE)
cat("Validation results saved: results/gse135251_results.csv\n")

# ============================================================
# STEP 8 — Significant genes (padj_mle < 0.01, |lfc_mle| > 2)
# ============================================================
cat("\n=== STEP 8: Significant genes (padj<0.01, |MLE LFC|>2) ===\n")

sig_val <- res_val_df %>%
  filter(!is.na(padj_mle), !is.na(lfc_mle), padj_mle < 0.01, abs(lfc_mle) > 2)

n_up_val   <- sum(sig_val$lfc_mle > 0)
n_down_val <- sum(sig_val$lfc_mle < 0)
cat("Up:   ", n_up_val,   "\n")
cat("Down: ", n_down_val, "\n")
cat("Total:", nrow(sig_val), "\n")

cat("\n=== TREM2/SPP1/GPNMB in GSE135251 ===\n")
for (g in c("TREM2", "SPP1", "GPNMB")) {
  row <- res_val_df %>% filter(gene_symbol == g)
  if (nrow(row) == 0) { cat(g, ": not found\n"); next }
  r <- row[1, ]
  in_sig <- !is.na(r$padj_mle) && r$padj_mle < 0.01 && abs(r$lfc_mle) > 2
  cat(sprintf("%s: MLE log2FC=%+.3f | apeglm=%+.3f | padj=%.2e | sig=%s\n",
              g, r$lfc_mle, r$lfc_apeglm, r$padj_mle, in_sig))
}

# ============================================================
# STEP 9 — Cross-cohort comparison
# ============================================================
cat("\n=== STEP 9: Cross-cohort comparison ===\n")

disc_sig  <- read.csv(file.path(disc_dir, "results/significant_genes.csv"))
res_disc  <- read.csv(file.path(disc_dir, "results/deseq2_results.csv"))

disc_ids <- disc_sig$ensembl_id
val_ids  <- sig_val$ensembl_id
overlap  <- intersect(disc_ids, val_ids)

cat("Discovery significant genes:  ", length(disc_ids), "\n")
cat("Validation significant genes: ", length(val_ids),  "\n")
cat("Overlap:                      ", length(overlap),  "\n")
cat("Overlap % of discovery:       ",
    round(100 * length(overlap) / length(disc_ids), 1), "%\n")

# Fisher's exact test (universe = genes tested in both cohorts)
universe_genes <- intersect(res_disc$ensembl_id, res_val_df$ensembl_id)
a  <- length(overlap)
b  <- length(setdiff(disc_ids, val_ids))
cv <- length(setdiff(val_ids, disc_ids))
d  <- length(universe_genes) - a - b - cv
ft <- fisher.test(matrix(c(a, b, cv, d), nrow = 2), alternative = "greater")
cat(sprintf("Fisher's exact: OR = %.2f | p = %.3e\n", ft$estimate, ft$p.value))

# Direction concordance
ov_disc <- disc_sig %>% filter(ensembl_id %in% overlap) %>%
  select(ensembl_id, lfc_disc = lfc_mle, gene_symbol)
ov_val  <- sig_val  %>% filter(ensembl_id %in% overlap) %>%
  select(ensembl_id, lfc_val = lfc_mle)
ov_df   <- inner_join(ov_disc, ov_val, by = "ensembl_id")
concordant <- mean(sign(ov_df$lfc_disc) == sign(ov_df$lfc_val)) * 100
cat(sprintf("Direction concordance: %.1f%%\n", concordant))

# Correlation (sig overlap)
cor_p <- cor.test(ov_df$lfc_disc, ov_df$lfc_val, method = "pearson")
cor_s <- cor.test(ov_df$lfc_disc, ov_df$lfc_val, method = "spearman")
cat(sprintf("Pearson  r (sig overlap):   %.3f (p=%.3e)\n",
            cor_p$estimate, cor_p$p.value))
cat(sprintf("Spearman rho (sig overlap): %.3f (p=%.3e)\n",
            cor_s$estimate, cor_s$p.value))

# Genome-wide correlation (all shared tested genes)
lfc_both <- inner_join(
  res_disc  %>% filter(!is.na(lfc_mle)) %>% select(ensembl_id, lfc_disc = lfc_mle, gene_symbol),
  res_val_df %>% filter(!is.na(lfc_mle)) %>% select(ensembl_id, lfc_val  = lfc_mle),
  by = "ensembl_id"
)
cor_gw <- cor.test(lfc_both$lfc_disc, lfc_both$lfc_val, method = "pearson")
cat(sprintf("Pearson r genome-wide (%d genes): %.3f\n",
            nrow(lfc_both), cor_gw$estimate))

cat("\nTREM2/SPP1/GPNMB cross-cohort:\n")
for (g in c("TREM2", "SPP1", "GPNMB")) {
  d_r   <- res_disc    %>% filter(gene_symbol == g)
  v_r   <- res_val_df  %>% filter(gene_symbol == g)
  d_sig <- nrow(d_r) > 0 && !is.na(d_r$padj_mle[1]) &&
           d_r$padj_mle[1] < 0.01 && abs(d_r$lfc_mle[1]) > 2
  v_sig <- nrow(v_r) > 0 && !is.na(v_r$padj_mle[1]) &&
           v_r$padj_mle[1] < 0.01 && abs(v_r$lfc_mle[1]) > 2
  cat(sprintf(
    "%s: disc log2FC=%s padj=%s sig=%s | val log2FC=%s padj=%s sig=%s\n", g,
    if (nrow(d_r) > 0) sprintf("%+.3f", d_r$lfc_mle[1]) else "NA",
    if (nrow(d_r) > 0) sprintf("%.2e",  d_r$padj_mle[1]) else "NA", d_sig,
    if (nrow(v_r) > 0) sprintf("%+.3f", v_r$lfc_mle[1]) else "NA",
    if (nrow(v_r) > 0) sprintf("%.2e",  v_r$padj_mle[1]) else "NA", v_sig
  ))
}

# ============================================================
# STEP 10 — Overlap plot
# ============================================================
cat("\n=== STEP 10: Overlap plot ===\n")

venn_df <- data.frame(
  Category = c("Discovery only", "Both cohorts", "Validation only"),
  Count    = c(length(disc_ids) - length(overlap), length(overlap),
               length(val_ids)  - length(overlap)),
  Group    = c("disc", "both", "val")
)
p_venn <- ggplot(venn_df, aes(x = reorder(Category, -Count), y = Count, fill = Group)) +
  geom_col(width = 0.55) +
  geom_text(aes(label = Count), vjust = -0.5, size = 5, fontface = "bold") +
  scale_fill_manual(values = c(disc = "#1565C0", both = "#6A1B9A", val = "#C62828")) +
  labs(
    title    = "Overlap: Discovery vs Validation Significant Genes",
    subtitle = sprintf("GSE162694 (n=%d) vs GSE135251 (n=%d) | padj<0.01 & |MLE LFC|>2 | %d shared",
                       length(disc_ids), length(val_ids), length(overlap)),
    x = NULL, y = "Gene count"
  ) +
  theme_bw(base_size = 13) +
  theme(plot.title = element_text(face = "bold"), legend.position = "none",
        axis.text.x = element_text(size = 11))
ggsave("plots/validation_overlap.png", p_venn, width = 7, height = 5, dpi = 150)
cat("Saved: plots/validation_overlap.png\n")

# ============================================================
# STEP 11 — LFC correlation scatter plot
# ============================================================
cat("\n=== STEP 11: LFC correlation plot ===\n")

lfc_both <- lfc_both %>%
  mutate(
    cat = case_when(
      ensembl_id %in% overlap  ~ "Both sig",
      ensembl_id %in% disc_ids ~ "Discovery only",
      ensembl_id %in% val_ids  ~ "Validation only",
      TRUE                     ~ "Neither"
    )
  )
highlight_mks <- lfc_both %>% filter(gene_symbol %in% c("TREM2", "SPP1", "GPNMB"))

p_cor <- ggplot(lfc_both %>% filter(cat == "Neither"),
                aes(x = lfc_disc, y = lfc_val)) +
  geom_point(alpha = 0.12, size = 0.4, colour = "grey60") +
  geom_point(data = filter(lfc_both, cat == "Discovery only"),
             alpha = 0.6, size = 1, colour = "#1565C0") +
  geom_point(data = filter(lfc_both, cat == "Validation only"),
             alpha = 0.6, size = 1, colour = "#C62828") +
  geom_point(data = filter(lfc_both, cat == "Both sig"),
             alpha = 0.8, size = 1.8, colour = "#6A1B9A") +
  geom_point(data = highlight_mks, size = 3.5, shape = 18, colour = "#FFD600") +
  geom_text_repel(data = highlight_mks, aes(label = gene_symbol),
                  size = 3.5, colour = "black",
                  nudge_y = 0.25, segment.size = 0.3) +
  geom_smooth(data = filter(lfc_both, cat == "Both sig"),
              method = "lm", colour = "#6A1B9A", linewidth = 0.8, se = TRUE) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "grey40") +
  labs(
    title    = "LFC Correlation: Discovery vs Validation",
    subtitle = sprintf(
      "Shared sig genes (n=%d, purple) | Pearson r=%.3f | Spearman ρ=%.3f\nGrey=neither sig | Blue=disc only | Red=val only | Yellow=TREM2/SPP1/GPNMB",
      nrow(ov_df), cor_p$estimate, cor_s$estimate),
    x = "log2FC MLE — Discovery (GSE162694)",
    y = "log2FC MLE — Validation (GSE135251)"
  ) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))
ggsave("plots/lfc_correlation.png", p_cor, width = 7.5, height = 6.5, dpi = 150)
cat("Saved: plots/lfc_correlation.png\n")

# ============================================================
# STEP 12 — Save overlap summary CSV
# ============================================================
cat("\n=== STEP 12: Save overlap summary ===\n")

summary_df <- data.frame(
  metric = c(
    "discovery_sig_genes", "validation_sig_genes", "overlap_genes",
    "overlap_pct_of_discovery", "universe_shared_tested_genes",
    "fishers_exact_OR", "fishers_exact_p",
    "direction_concordance_pct",
    "pearson_r_sig_overlap", "pearson_p_sig_overlap",
    "spearman_rho_sig_overlap", "spearman_p_sig_overlap",
    "pearson_r_genomewide",
    "TREM2_sig_discovery", "TREM2_sig_validation",
    "SPP1_sig_discovery",  "SPP1_sig_validation",
    "GPNMB_sig_discovery", "GPNMB_sig_validation"
  ),
  value = c(
    length(disc_ids), length(val_ids), length(overlap),
    round(100 * length(overlap) / length(disc_ids), 2),
    length(universe_genes),
    round(ft$estimate, 3), signif(ft$p.value, 4),
    round(concordant, 2),
    round(cor_p$estimate, 4), signif(cor_p$p.value, 4),
    round(cor_s$estimate, 4), signif(cor_s$p.value, 4),
    round(cor_gw$estimate, 4),
    any(disc_sig$gene_symbol == "TREM2"),
    any(sig_val$gene_symbol  == "TREM2"),
    any(disc_sig$gene_symbol == "SPP1"),
    any(sig_val$gene_symbol  == "SPP1"),
    any(disc_sig$gene_symbol == "GPNMB"),
    any(sig_val$gene_symbol  == "GPNMB")
  )
)
write.csv(summary_df, "results/validation_overlap_summary.csv", row.names = FALSE)
cat("Saved: results/validation_overlap_summary.csv\n")

cat("\n============================================================\n")
cat("VALIDATION COMPLETE — GSE135251\n")
cat(sprintf("  Sig genes: %d up | %d down | %d total\n",
            n_up_val, n_down_val, nrow(sig_val)))
cat(sprintf("  Overlap with discovery: %d genes\n", length(overlap)))
cat(sprintf("  Fisher's exact: OR=%.2f  p=%.3e\n", ft$estimate, ft$p.value))
cat(sprintf("  Concordance: %.1f%% | Spearman rho=%.3f\n",
            concordant, cor_s$estimate))
cat("============================================================\n")
