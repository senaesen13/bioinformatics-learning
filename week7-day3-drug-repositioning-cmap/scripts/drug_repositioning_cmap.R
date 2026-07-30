#!/usr/bin/env Rscript
## ==============================================================================
## Script: drug_repositioning_cmap.R
## Method: LINCS CMap L1000 Signature Reversal & Drug Repositioning
##
## NOTE on Tau scores: The Tau values in this script are SIMULATED PLACEHOLDERS.
## They represent the expected direction and approximate magnitude of CMap
## connectivity for these eight known NAFLD-relevant compounds, based on their
## published mechanisms of action. They are NOT from a real clue.io query.
## The .grp files in results/ are correctly formatted for submission to the
## Broad CMap portal (clue.io) for live-validated Tau scores.
## ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

# Self-contained path resolution — all paths relative to this module folder
get_module_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg  <- grep("--file=", cmd_args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(dirname(normalizePath(sub("--file=", "", file_arg[1])))))
  }
  if (exists(".rs.getScriptPath", mode = "function")) {
    sp <- .rs.getScriptPath()
    if (!is.null(sp) && nchar(sp) > 0) return(dirname(dirname(normalizePath(sp))))
  }
  return(getwd())
}

module_dir <- get_module_dir()
repo_dir   <- dirname(module_dir)
out_dir    <- file.path(module_dir, "results")
plot_dir   <- file.path(module_dir, "plots")
if (!dir.exists(out_dir))  dir.create(out_dir,  recursive = TRUE)
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

cat("============================================================\n")
cat("LINCS CMap L1000 Drug Repositioning Analysis\n")
cat("Module Directory:", module_dir, "\n")
cat("============================================================\n\n")

# 1. Load DESeq2 Results from GSE130970 (folder 03 in this repo)
res_file <- file.path(repo_dir, "03-NAFLD-second-validation-GSE130970", "results", "gse130970_results.csv")
if (!file.exists(res_file)) {
  # Optional fallback for running from a different working directory
  args <- commandArgs(trailingOnly = TRUE)
  deseq2_flag <- grep("^--deseq2=", args, value = TRUE)
  res_file <- if (length(deseq2_flag) > 0) sub("^--deseq2=", "", deseq2_flag[1]) else ""
}

if (file.exists(res_file)) {
  cat("[INFO] Loading DESeq2 results from:", res_file, "\n")
  df_raw   <- read.csv(res_file, stringsAsFactors = FALSE)
  sym_col  <- if ("gene_symbol" %in% colnames(df_raw)) "gene_symbol" else "symbol"
  lfc_col  <- if ("lfc_mle" %in% colnames(df_raw)) "lfc_mle" else if ("log2FoldChange" %in% colnames(df_raw)) "log2FoldChange" else "lfc"
  padj_col <- if ("padj_mle" %in% colnames(df_raw)) "padj_mle" else if ("padj" %in% colnames(df_raw)) "padj" else "padj"

  df_clean <- df_raw %>%
    filter(!is.na(.data[[sym_col]]), .data[[sym_col]] != "") %>%
    filter(!is.na(.data[[lfc_col]]), !is.na(.data[[padj_col]])) %>%
    mutate(gene_symbol = .data[[sym_col]], lfc = .data[[lfc_col]], padj = .data[[padj_col]])
} else {
  stop("DESeq2 results file not found. Run from the repo root or pass --deseq2=<path>.")
}

# 2. Extract Top 150 UP and Top 150 DOWN DEGs for CMap Query
up_genes <- df_clean %>%
  filter(padj < 0.05, lfc > 0.5) %>%
  arrange(desc(lfc)) %>%
  pull(gene_symbol) %>%
  unique() %>%
  head(150)

down_genes <- df_clean %>%
  filter(padj < 0.05, lfc < -0.5) %>%
  arrange(lfc) %>%
  pull(gene_symbol) %>%
  unique() %>%
  head(150)

cat(sprintf("CMap query signature -> UP genes: %d | DOWN genes: %d\n", length(up_genes), length(down_genes)))
cat(sprintf("Overlap between UP and DOWN lists: %d (should be 0)\n", length(intersect(up_genes, down_genes))))

# 3. Write Broad Institute CMap .grp query files (one gene per line, no header)
write.table(up_genes,   file.path(out_dir, "cmap_up_genes.grp"),   row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(down_genes, file.path(out_dir, "cmap_down_genes.grp"), row.names = FALSE, col.names = FALSE, quote = FALSE)
cat("[INFO] .grp query files written to results/\n")
cat("       Submit these to clue.io -> Query CMap -> L1000 for real Tau scores.\n\n")

# 4. SIMULATED Tau scores (placeholder — see note at top of script)
# These values represent the expected direction and approximate magnitude for
# each compound based on published NAFLD pharmacology. They are NOT from a real
# LINCS query. Replace with clue.io output once a real query has been run.
reference_drugs <- data.frame(
  Compound         = c("Obeticholic Acid", "Resmetirom", "Pioglitazone", "Metformin",
                       "Fenofibrate", "Vorinostat", "Pictilisib", "Rapamycin"),
  Mechanism        = c("FXR Agonist", "THR-beta Agonist", "PPAR-gamma Agonist", "AMPK Activator",
                       "PPAR-alpha Agonist", "HDAC Inhibitor", "PI3K Inhibitor", "mTOR Inhibitor"),
  Target_Class     = c("Nuclear Receptor", "Nuclear Receptor", "Nuclear Receptor", "Kinase",
                       "Nuclear Receptor", "Epigenetic", "Kinase", "Kinase"),
  Raw_Connectivity = c(-0.89, -0.85, -0.82, -0.78, -0.75, -0.71, -0.68, -0.64),
  Tau_Score        = c(-98.5, -95.2, -91.4, -86.8, -82.1, -77.5, -73.2, -68.9),
  Tau_Source       = rep("SIMULATED_PLACEHOLDER", 8),
  stringsAsFactors = FALSE
)

cat("[NOTE] Tau scores below are SIMULATED PLACEHOLDERS, not from a real LINCS query.\n")
print(reference_drugs[, c("Compound", "Mechanism", "Tau_Score", "Tau_Source")])

# 5. Export results and bar chart
write.csv(reference_drugs, file.path(out_dir, "drug_repositioning_candidates.csv"), row.names = FALSE)

p_drug <- ggplot(reference_drugs, aes(x = reorder(Compound, -Tau_Score), y = Tau_Score, fill = Mechanism)) +
  geom_bar(stat = "identity", width = 0.65, color = "black") +
  geom_text(aes(label = sprintf("Tau = %.1f", Tau_Score)), hjust = 1.1, color = "white", fontface = "bold", size = 3.8) +
  coord_flip() +
  scale_fill_brewer(palette = "Set2") +
  theme_bw(base_size = 12) +
  labs(
    title = "LINCS CMap L1000 Candidate Reversal Compounds (SIMULATED Tau)",
    subtitle = "Negative Tau = drug reverses NAFLD transcriptomic signature | Scores are placeholder values",
    x = "Drug Compound",
    y = "Connectivity Tau Score (-100 to +100)",
    fill = "Mechanism of Action"
  ) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(plot_dir, "lincs_drug_repositioning_tau_scores.png"), p_drug, width = 8, height = 5.5, dpi = 300)
cat("\n[SUCCESS] Drug repositioning pipeline completed. Results in:", out_dir, "\n")
