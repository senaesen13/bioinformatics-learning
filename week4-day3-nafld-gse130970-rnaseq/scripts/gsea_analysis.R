## ============================================================
## Week 4 Day 3 — GSEA: KEGG + Hallmark
## Dataset: GSE130970 (Entrez gene IDs — no conversion needed)
## Input:  results/gse130970_results.csv (apeglm-shrunken LFC)
## Output: results/gsea_kegg_results.csv
##         results/gsea_hallmark_results.csv
##         plots/gsea_kegg_dotplot.png  plots/gsea_kegg_ridgeplot.png
##         plots/gsea_hallmark_dotplot.png  plots/gsea_hallmark_ridgeplot.png
## Run from: week4-day3-nafld-gse130970-rnaseq/
## ============================================================

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(enrichplot)
  library(msigdbr)
  library(org.Hs.eg.db)
  library(dplyr)
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

# ============================================================
# STEP 1 — Build ranked gene list (Entrez IDs already available)
# ============================================================
cat("\n=== STEP 1: Load DESeq2 results ===\n")

res_df <- read.csv("results/gse130970_results.csv", stringsAsFactors = FALSE)
cat("Loaded", nrow(res_df), "genes\n")

res_clean <- res_df %>%
  filter(!is.na(lfc_apeglm), !is.na(entrez_id), entrez_id != "") %>%
  arrange(desc(lfc_apeglm)) %>%
  distinct(entrez_id, .keep_all = TRUE)

cat("Genes with valid apeglm LFC and Entrez ID:", nrow(res_clean), "\n")

# Build named numeric vector: Entrez ID -> apeglm LFC
gene_list_entrez <- setNames(res_clean$lfc_apeglm, res_clean$entrez_id)

cat("Final ranked list length:", length(gene_list_entrez), "\n")
cat("LFC range: [", round(min(gene_list_entrez), 3),
    ",", round(max(gene_list_entrez), 3), "]\n")

# ============================================================
# STEP 2 — KEGG GSEA
# ============================================================
cat("\n=== STEP 2: KEGG GSEA (gseKEGG) ===\n")

gsea_kegg <- gseKEGG(
  geneList      = gene_list_entrez,
  organism      = "hsa",
  minGSSize     = 15,
  maxGSSize     = 500,
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  verbose       = FALSE,
  seed          = TRUE
)
cat("KEGG significant pathways (padj<0.05):", nrow(gsea_kegg@result), "\n")

kegg_res_df <- as.data.frame(gsea_kegg)
write.csv(kegg_res_df, "results/gsea_kegg_results.csv", row.names = FALSE)
cat("Saved: results/gsea_kegg_results.csv\n")

if (nrow(kegg_res_df) > 0) {
  cat("\nTop activated KEGG pathways in NAFLD (NES > 0):\n")
  print(kegg_res_df %>% filter(NES > 0) %>% arrange(p.adjust) %>% head(10) %>%
        select(Description, NES, pvalue, p.adjust, setSize))

  cat("\nTop suppressed KEGG pathways in NAFLD (NES < 0):\n")
  print(kegg_res_df %>% filter(NES < 0) %>% arrange(p.adjust) %>% head(10) %>%
        select(Description, NES, pvalue, p.adjust, setSize))

  n_kegg_show <- min(20, nrow(kegg_res_df))

  tryCatch({
    p_kegg_dot <- dotplot(gsea_kegg, showCategory = n_kegg_show, split = ".sign") +
      facet_grid(~.sign) +
      ggtitle("KEGG GSEA — NAFLD vs Normal",
              subtitle = "GSE130970 | ranked by apeglm log2FC | padj<0.05") +
      theme_bw(base_size = 11) +
      theme(plot.title  = element_text(face = "bold", size = 13),
            axis.text.y = element_text(size = 9),
            strip.text  = element_text(face = "bold"))
    ggsave("plots/gsea_kegg_dotplot.png", p_kegg_dot,
           width = 14, height = max(6, n_kegg_show * 0.4 + 3), dpi = 150)
    cat("Saved: plots/gsea_kegg_dotplot.png\n")
  }, error = function(e) cat("KEGG dotplot failed:", conditionMessage(e), "\n"))

  tryCatch({
    p_kegg_ridge <- ridgeplot(gsea_kegg, showCategory = n_kegg_show, fill = "p.adjust") +
      labs(title    = "KEGG GSEA Ridgeplot — NAFLD vs Normal",
           subtitle = "Core-enrichment gene log2FC distribution | padj<0.05",
           x        = "log2 Fold Change (apeglm, NAFLD / Normal)") +
      theme_bw(base_size = 11) +
      theme(plot.title  = element_text(face = "bold", size = 13),
            axis.text.y = element_text(size = 9))
    ggsave("plots/gsea_kegg_ridgeplot.png", p_kegg_ridge,
           width = 12, height = max(7, n_kegg_show * 0.45 + 3), dpi = 150)
    cat("Saved: plots/gsea_kegg_ridgeplot.png\n")
  }, error = function(e) cat("KEGG ridgeplot failed:", conditionMessage(e), "\n"))

} else {
  cat("No significant KEGG pathways found.\n")
}

# ============================================================
# STEP 3 — Hallmark GSEA (MSigDB H collection)
# ============================================================
cat("\n=== STEP 3: Hallmark GSEA (MSigDB) ===\n")

hallmark_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, entrez_gene) %>%
  dplyr::mutate(entrez_gene = as.character(entrez_gene))

cat("Hallmark gene sets:", length(unique(hallmark_t2g$gs_name)), "\n")

gsea_hallmark <- GSEA(
  geneList      = gene_list_entrez,
  TERM2GENE     = hallmark_t2g,
  minGSSize     = 15,
  maxGSSize     = 500,
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  verbose       = FALSE,
  seed          = TRUE
)
cat("Hallmark significant gene sets (padj<0.05):", nrow(gsea_hallmark@result), "\n")

hallmark_res_df <- as.data.frame(gsea_hallmark)
write.csv(hallmark_res_df, "results/gsea_hallmark_results.csv", row.names = FALSE)
cat("Saved: results/gsea_hallmark_results.csv\n")

if (nrow(hallmark_res_df) > 0) {
  cat("\nTop activated Hallmark gene sets (NES > 0):\n")
  print(hallmark_res_df %>% filter(NES > 0) %>% arrange(p.adjust) %>% head(10) %>%
        select(Description, NES, pvalue, p.adjust, setSize))

  cat("\nTop suppressed Hallmark gene sets (NES < 0):\n")
  print(hallmark_res_df %>% filter(NES < 0) %>% arrange(p.adjust) %>% head(10) %>%
        select(Description, NES, pvalue, p.adjust, setSize))

  n_hall_show <- min(20, nrow(hallmark_res_df))

  gsea_hallmark_clean <- gsea_hallmark
  gsea_hallmark_clean@result$Description <- sub("^HALLMARK_", "",
    gsea_hallmark_clean@result$Description)

  tryCatch({
    p_hall_dot <- dotplot(gsea_hallmark_clean, showCategory = n_hall_show,
                          split = ".sign") +
      facet_grid(~.sign) +
      ggtitle("Hallmark GSEA — NAFLD vs Normal",
              subtitle = "MSigDB H collection | GSE130970 | padj<0.05") +
      theme_bw(base_size = 11) +
      theme(plot.title  = element_text(face = "bold", size = 13),
            axis.text.y = element_text(size = 9),
            strip.text  = element_text(face = "bold"))
    ggsave("plots/gsea_hallmark_dotplot.png", p_hall_dot,
           width = 14, height = max(6, n_hall_show * 0.4 + 3), dpi = 150)
    cat("Saved: plots/gsea_hallmark_dotplot.png\n")
  }, error = function(e) cat("Hallmark dotplot failed:", conditionMessage(e), "\n"))

  tryCatch({
    p_hall_ridge <- ridgeplot(gsea_hallmark_clean, showCategory = n_hall_show,
                               fill = "p.adjust") +
      labs(title    = "Hallmark GSEA Ridgeplot — NAFLD vs Normal",
           subtitle = "Core-enrichment gene log2FC distribution | padj<0.05",
           x        = "log2 Fold Change (apeglm, NAFLD / Normal)") +
      theme_bw(base_size = 11) +
      theme(plot.title  = element_text(face = "bold", size = 13),
            axis.text.y = element_text(size = 9))
    ggsave("plots/gsea_hallmark_ridgeplot.png", p_hall_ridge,
           width = 12, height = max(7, n_hall_show * 0.45 + 3), dpi = 150)
    cat("Saved: plots/gsea_hallmark_ridgeplot.png\n")
  }, error = function(e) cat("Hallmark ridgeplot failed:", conditionMessage(e), "\n"))

} else {
  cat("No significant Hallmark gene sets found.\n")
}

cat("\n============================================================\n")
cat("GSEA COMPLETE — GSE130970\n")
cat("============================================================\n")
cat("KEGG significant pathways:     ", nrow(kegg_res_df), "\n")
cat("Hallmark significant gene sets:", nrow(hallmark_res_df), "\n")
