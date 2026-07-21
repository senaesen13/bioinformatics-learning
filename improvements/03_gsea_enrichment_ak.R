#!/usr/bin/env Rscript
## ==============================================================================
## Script: 03_gsea_enrichment_ak.R
## Author: Dr. Ali Kaynar (King's College London)
## Description: Standardized DEG Functional Enrichment, ORA, and GSEA Rank File Pipeline
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
out_dir  <- file.path(root_dir, "improvements", "results_03_enrichment_gsea")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

cat("============================================================\n")
cat("AK Pipeline 03: Standardized DEG Enrichment & GSEA Pipeline\n")
cat("Root Directory:", root_dir, "\n")
cat("Output Directory:", out_dir, "\n")
cat("============================================================\n\n")

# 1. Load Real DESeq2 Results or Create Synthetic Dataset
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
  cat("[INFO] Loading real DESeq2 results from:", res_file, "\n")
  df_raw <- read.csv(res_file, stringsAsFactors = FALSE)
  cols   <- colnames(df_raw)
  
  sym_col  <- find_col(cols, c("gene_symbol", "symbol", "Gene", "gene"))
  lfc_col  <- find_col(cols, c("lfc_mle", "log2FoldChange", "lfc", "lfc_apeglm"))
  pval_col <- find_col(cols, c("pvalue_mle", "pvalue", "p.value", "pval"))
  padj_col <- find_col(cols, c("padj_mle", "padj", "adj.P.Val", "FDR"))
  
  df_clean <- df_raw %>%
    filter(!is.na(.data[[sym_col]]), .data[[sym_col]] != "") %>%
    filter(!is.na(.data[[lfc_col]]), !is.na(.data[[pval_col]])) %>%
    mutate(
      gene_symbol = .data[[sym_col]],
      log2FC = .data[[lfc_col]],
      pvalue = pmax(.data[[pval_col]], 1e-300),
      padj = if (!is.null(padj_col)) .data[[padj_col]] else p.adjust(.data[[pval_col]], method = "BH")
    )
} else {
  cat("[INFO] Generating synthetic DESeq2 dataset for GSEA validation...\n")
  set.seed(99)
  n_genes <- 2000
  df_clean <- data.frame(
    gene_symbol = paste0("GENE_", sprintf("%04d", 1:n_genes)),
    log2FC = rnorm(n_genes, mean = 0, sd = 1.5),
    pvalue = runif(n_genes, 1e-8, 0.5),
    stringsAsFactors = FALSE
  ) %>% mutate(padj = p.adjust(pvalue, method = "BH"))
}

# 2. Standardized GSEA Ranking Metric: sign(log2FC) * -log10(pvalue)
cat("\n=== Calculating Standardized GSEA Ranking Metric ===\n")
df_ranked <- df_clean %>%
  mutate(rank_metric = sign(log2FC) * -log10(pvalue)) %>%
  group_by(gene_symbol) %>%
  slice_max(order_by = abs(rank_metric), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(desc(rank_metric))

# 3. Export Standard GSEA Rank File (.rnk format for Desktop GSEA / FGSEA)
rnk_file <- file.path(out_dir, "gsea_ranked_genes.rnk")
write.table(
  df_ranked %>% select(gene_symbol, rank_metric),
  file = rnk_file,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)
cat("[SUCCESS] Saved GSEA Desktop compatible rank file (.rnk):", rnk_file, "\n")

# 4. Over-Representation Analysis (ORA) Hypergeometric Test
up_degs   <- df_clean %>% filter(padj < 0.05, log2FC > 1) %>% pull(gene_symbol) %>% unique()
down_degs <- df_clean %>% filter(padj < 0.05, log2FC < -1) %>% pull(gene_symbol) %>% unique()
universe  <- unique(df_clean$gene_symbol)

cat(sprintf("Significance Thresholds (padj < 0.05, |log2FC| > 1) -> UP: %d | DOWN: %d | Universe: %d\n",
            length(up_degs), length(down_degs), length(universe)))

# Hallmark & Pathway Definitions
hallmark_pathways <- list(
  "HALLMARK_FATTY_ACID_METABOLISM" = c("ACADM", "ACADVL", "CPT1A", "CPT2", "FASN", "SCD", "PPARA", "HADHA", paste0("GENE_", sprintf("%04d", 1:15))),
  "HALLMARK_INFLAMMATORY_RESPONSE" = c("TNF", "IL6", "CCL2", "CD68", "TREM2", "SPP1", "NFKB1", "CXCL8", paste0("GENE_", sprintf("%04d", 16:30))),
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION" = c("COX4I1", "NDUFA9", "SDHA", "UQCRC1", "ATP5F1A", paste0("GENE_", sprintf("%04d", 31:45))),
  "HALLMARK_TGF_BETA_SIGNALING" = c("TGFB1", "SMAD2", "SMAD3", "COL1A1", "COL3A1", "ACTA2", paste0("GENE_", sprintf("%04d", 46:60)))
)

ora_runner <- function(deg_set, set_name = "UP") {
  lapply(names(hallmark_pathways), function(path_name) {
    p_genes <- hallmark_pathways[[path_name]]
    k <- length(intersect(deg_set, p_genes))
    m <- length(p_genes)
    n <- length(deg_set)
    N <- length(universe)
    p_val <- phyper(k - 1, m, N - m, n, lower.tail = FALSE)
    
    data.frame(
      Direction = set_name,
      Pathway = path_name,
      Overlap_Count = k,
      Pathway_Size = m,
      Gene_Ratio = round(k / max(1, m), 3),
      Pvalue = signif(p_val, 4),
      Padj = signif(p.adjust(p_val, method = "BH"), 4),
      Overlapping_Genes = paste(head(intersect(deg_set, p_genes), 5), collapse = "; "),
      stringsAsFactors = FALSE
    )
  }) %>% bind_rows()
}

ora_up   <- ora_runner(up_degs, "UP")
ora_down <- ora_runner(down_degs, "DOWN")
ora_all  <- rbind(ora_up, ora_down)

print(ora_all)

# 5. Export ORA & GSEA Results
write.csv(ora_all, file.path(out_dir, "ora_pathway_enrichment.csv"), row.names = FALSE)
write.csv(df_ranked, file.path(out_dir, "gsea_full_ranked_table.csv"), row.names = FALSE)

# 6. Visualization Plot
p_ora <- ggplot(ora_all %>% filter(Overlap_Count > 0), 
                aes(x = Gene_Ratio, y = Pathway, size = Overlap_Count, color = -log10(Pvalue))) +
  geom_point(alpha = 0.85) +
  facet_wrap(~ Direction, scales = "free_y") +
  scale_color_gradient(low = "#1E88E5", high = "#D81B60") +
  theme_bw(base_size = 12) +
  labs(title = "Standardized Pathway Over-Representation Analysis (ORA)",
       x = "Gene Ratio (Overlap / Pathway Size)", y = "", color = "-log10(p)", size = "Overlap Count") +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(out_dir, "ora_enrichment_dotplot.png"), p_ora, width = 8.5, height = 5.5, dpi = 300)
cat("[SUCCESS] Pipeline 03 completed. Results saved to:", out_dir, "\n")
