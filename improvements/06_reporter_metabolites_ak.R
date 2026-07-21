#!/usr/bin/env Rscript
## ==============================================================================
## Script: 06_reporter_metabolites_ak.R
## Author: Dr. Ali Kaynar (King's College London)
## Description: Patil & Nielsen Reporter Metabolite Pipeline with GEM Model Input
## ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

# Environment & Path Auto-Detection
get_root_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", cmd_args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(dirname(normalizePath(sub("--file=", "", file_arg[1])))))
  }
  if (exists(".rs.getScriptPath", mode = "function")) {
    sp <- .rs.getScriptPath()
    if (!is.null(sp) && nchar(sp) > 0) return(dirname(dirname(normalizePath(sp))))
  }
  return(getwd())
}

root_dir <- get_root_dir()
out_dir  <- file.path(root_dir, "improvements", "results_06_reporter_metabolites")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

cat("============================================================\n")
cat("AK Pipeline 06: GEM-Based Reporter Metabolite Analysis\n")
cat("Root Directory:", root_dir, "\n")
cat("Output Directory:", out_dir, "\n")
cat("============================================================\n\n")

# 1. Accept GEM Model File Input or Default to Extracted Human-GEM Network
args <- commandArgs(trailingOnly = TRUE)
gem_input_file <- if (length(args) > 0) args[1] else file.path(root_dir, "improvements", "data", "human_gem_topology_network.csv")

if (!file.exists(gem_input_file)) {
  # Fallback workspace paths
  alt_gem <- "/Users/k2254978/Desktop/Work/ToolBoxProgs/Human_GEM_2.0/model/Human-GEM.xml"
  if (file.exists(alt_gem)) {
    cat("[INFO] Running GEM XML Parser on:", alt_gem, "\n")
    system(sprintf("python3 %s/improvements/parse_gem.py", root_dir))
  }
  gem_input_file <- file.path(root_dir, "improvements", "data", "human_gem_topology_network.csv")
}

cat("[INFO] Loading GEM Model Topology Input from:", gem_input_file, "\n")

if (file.exists(gem_input_file)) {
  gem_map <- read.csv(gem_input_file, stringsAsFactors = FALSE)
} else {
  stop("GEM model input file not found.")
}

# 2. Currency Metabolites & Noise Filtering
currency_metabolites <- c("H2O", "ATP", "ADP", "AMP", "NAD+", "NADH", "NADP+", "NADPH",
                          "H+", "Pi", "PPi", "CoA", "CO2", "O2", "HCO3-", "Na+", "K+", "Cl-", "water", "oxygen", "phosphate")

gem_map_clean <- gem_map %>%
  filter(!is.na(metabolite_name), metabolite_name != "") %>%
  filter(!is.na(gene_symbol), gene_symbol != "") %>%
  filter(!tolower(metabolite_name) %in% tolower(currency_metabolites))

cat(sprintf("GEM Topology: Extracted %d metabolite-gene associations across %d unique metabolites.\n",
            nrow(gem_map_clean), length(unique(gem_map_clean$metabolite_name))))

# 3. Load DESeq2 Results (Real or Synthetic Fallback)
res_file <- file.path(root_dir, "week4-day3-nafld-gse130970-rnaseq", "results", "gse130970_results.csv")
if (!file.exists(res_file)) {
  res_file <- file.path(root_dir, "week4-day1-nafld-bulk-rnaseq", "results", "deseq2_results.csv")
}

find_col <- function(cols, candidates) {
  match <- intersect(candidates, cols)
  if (length(match) > 0) return(match[1])
  return(NULL)
}

if (file.exists(res_file)) {
  cat("[INFO] Loading differential gene expression results from:", res_file, "\n")
  df_raw <- read.csv(res_file, stringsAsFactors = FALSE)
  cols   <- colnames(df_raw)
  
  sym_col  <- find_col(cols, c("gene_symbol", "symbol", "Gene", "gene"))
  pval_col <- find_col(cols, c("pvalue_mle", "pvalue", "p.value", "pval"))
  lfc_col  <- find_col(cols, c("lfc_mle", "log2FoldChange", "lfc", "lfc_apeglm"))
  
  df_clean <- df_raw %>%
    filter(!is.na(.data[[sym_col]]), .data[[sym_col]] != "") %>%
    filter(!is.na(.data[[pval_col]])) %>%
    mutate(
      gene_symbol = .data[[sym_col]],
      pvalue = pmax(.data[[pval_col]], 1e-300),
      lfc = if (!is.null(lfc_col)) .data[[lfc_col]] else 0
    ) %>%
    group_by(gene_symbol) %>%
    slice_min(order_by = pvalue, n = 1, with_ties = FALSE) %>%
    ungroup()
} else {
  cat("[INFO] Generating synthetic gene p-value vector...\n")
  set.seed(77)
  all_genes <- unique(c(gem_map_clean$gene_symbol, paste0("GENE_", sprintf("%04d", 1:500))))
  df_clean <- data.frame(
    gene_symbol = all_genes,
    pvalue = runif(length(all_genes), 1e-5, 0.5),
    lfc = rnorm(length(all_genes), 0, 1.5),
    stringsAsFactors = FALSE
  )
}

# 4. Calculate Two-Tailed Gene Z-scores
df_clean <- df_clean %>%
  mutate(gene_z = qnorm(1 - pvalue / 2))

gene_z_vec   <- setNames(df_clean$gene_z, df_clean$gene_symbol)
gene_lfc_vec <- setNames(df_clean$lfc, df_clean$gene_symbol)
gene_universe <- names(gene_z_vec)

# 5. Patil & Nielsen Algorithm on GEM Nodes
# Filter metabolites to those with 3 to 100 neighboring genes in the dataset
metabolite_gene_counts <- gem_map_clean %>%
  filter(gene_symbol %in% gene_universe) %>%
  group_by(metabolite_name) %>%
  summarise(
    n_genes = n_distinct(gene_symbol),
    genes = list(unique(gene_symbol)),
    .groups = "drop"
  ) %>%
  filter(n_genes >= 3, n_genes <= 100)

cat(sprintf("Running Reporter Metabolites on %d qualified GEM nodes (3 <= k <= 100)...\n", nrow(metabolite_gene_counts)))

n_perm <- 1000

reporter_results <- lapply(1:nrow(metabolite_gene_counts), function(i) {
  met <- metabolite_gene_counts$metabolite_name[i]
  g_set <- metabolite_gene_counts$genes[[i]]
  k <- length(g_set)
  
  z_raw <- sum(gene_z_vec[g_set], na.rm = TRUE) / sqrt(k)
  avg_lfc <- mean(gene_lfc_vec[g_set], na.rm = TRUE)
  
  # Background Monte Carlo Permutations
  set.seed(42 + i)
  perm_z <- replicate(n_perm, {
    rand_g <- sample(gene_universe, k)
    sum(gene_z_vec[rand_g], na.rm = TRUE) / sqrt(k)
  })
  
  mu_k    <- mean(perm_z)
  sigma_k <- sd(perm_z)
  z_corr  <- if (sigma_k > 0) (z_raw - mu_k) / sigma_k else z_raw
  p_rep   <- pnorm(z_corr, lower.tail = FALSE)
  
  data.frame(
    Metabolite = met,
    Neighbor_Genes_Count = k,
    Mean_Log2FC = round(avg_lfc, 3),
    Z_Raw = round(z_raw, 3),
    Reporter_Z_Score = round(z_corr, 3),
    Pvalue = signif(p_rep, 4),
    Neighbor_Genes = paste(head(g_set, 8), collapse = "; "),
    stringsAsFactors = FALSE
  )
}) %>% bind_rows() %>%
  mutate(Padj = signif(p.adjust(Pvalue, method = "BH"), 4)) %>%
  arrange(desc(Reporter_Z_Score))

print(head(reporter_results, 15))

# 6. Export Results & Plots
write.csv(reporter_results, file.path(out_dir, "gem_reporter_metabolites_summary.csv"), row.names = FALSE)

top_plot_df <- head(reporter_results, 20)

p_rep <- ggplot(top_plot_df, aes(x = reorder(Metabolite, Reporter_Z_Score), y = Reporter_Z_Score, fill = Mean_Log2FC)) +
  geom_bar(stat = "identity", width = 0.7, color = "black") +
  coord_flip() +
  scale_fill_gradient2(low = "#1E88E5", mid = "grey90", high = "#D81B60", midpoint = 0, name = "Mean Log2FC") +
  theme_bw(base_size = 11) +
  labs(
    title = "Top 20 GEM Reporter Metabolites (Human-GEM 2.0 Input)",
    subtitle = "Patil & Nielsen Z-score correction on GEM stoichiometric topology",
    x = "Metabolite Node",
    y = "Corrected Reporter Z-score"
  ) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(out_dir, "gem_reporter_metabolites_zscores.png"), p_rep, width = 8, height = 6.5, dpi = 300)
cat("[SUCCESS] Pipeline 06 completed using GEM input. Results saved to:", out_dir, "\n")
