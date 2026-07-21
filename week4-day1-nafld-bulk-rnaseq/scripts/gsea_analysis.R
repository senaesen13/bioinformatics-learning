## ============================================================
## Week 4 Day 1 — GSEA: KEGG + Hallmark
## Input:  results/deseq2_results.csv (apeglm-shrunken LFC)
## Output: results/gsea_kegg_results.csv
##         results/gsea_hallmark_results.csv
##         plots/gsea_kegg_dotplot.png
##         plots/gsea_kegg_ridgeplot.png
##         plots/gsea_hallmark_dotplot.png
##         plots/gsea_hallmark_ridgeplot.png
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

# Run from project root
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
# STEP 1 — Build ranked gene list from DESeq2 apeglm results
# ============================================================
cat("\n=== STEP 1: Load DESeq2 results ===\n")

res_df <- read.csv("results/deseq2_results.csv")
cat("Loaded", nrow(res_df), "genes\n")

res_clean <- res_df %>%
  filter(!is.na(lfc_apeglm), !is.na(gene_symbol), gene_symbol != "",
         !is.na(pvalue_mle), pvalue_mle > 0)

cat("Genes with valid apeglm LFC, pvalue_mle, and symbol:", nrow(res_clean), "\n")

# ============================================================
# STEP 2 — Map to Entrez IDs (required for KEGG GSEA)
# ============================================================
cat("\n=== STEP 2: Map Ensembl -> Entrez IDs ===\n")

entrez_map <- mapIds(
  org.Hs.eg.db,
  keys      = res_clean$ensembl_id,
  column    = "ENTREZID",
  keytype   = "ENSEMBL",
  multiVals = "first"
)

res_clean$entrez_id <- entrez_map[res_clean$ensembl_id]
res_with_entrez <- res_clean %>% filter(!is.na(entrez_id))
cat("Genes with Entrez ID:", nrow(res_with_entrez), "\n")

# Ranked list: Entrez ID -> sign(lfc_apeglm) * -log10(pvalue_mle), sorted descending, deduplicated
gene_list_entrez <- res_with_entrez %>%
  mutate(rank_score = sign(lfc_apeglm) * -log10(pvalue_mle)) %>%
  arrange(desc(rank_score)) %>%
  distinct(entrez_id, .keep_all = TRUE) %>%
  { setNames(.$rank_score, .$entrez_id) }

cat("Final ranked list length:", length(gene_list_entrez), "\n")
cat("Rank score range: [", round(min(gene_list_entrez), 3),
    ",", round(max(gene_list_entrez), 3), "]\n")

# ============================================================
# STEP 3 — KEGG GSEA
# ============================================================
cat("\n=== STEP 3: KEGG GSEA (gseKEGG) ===\n")

gsea_kegg <- gseKEGG(
  geneList     = gene_list_entrez,
  organism     = "hsa",
  minGSSize    = 15,
  maxGSSize    = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose      = FALSE,
  seed         = TRUE
)
cat("KEGG GSEA: significant pathways (padj<0.05) =", nrow(gsea_kegg@result), "\n")

kegg_res_df <- as.data.frame(gsea_kegg)
write.csv(kegg_res_df, "results/gsea_kegg_results.csv", row.names = FALSE)
cat("Saved: results/gsea_kegg_results.csv\n")

if (nrow(kegg_res_df) > 0) {
  cat("\nTop activated KEGG pathways in NAFLD (NES > 0):\n")
  kegg_up <- kegg_res_df %>% filter(NES > 0) %>% arrange(p.adjust) %>% head(10)
  print(kegg_up[, c("Description", "NES", "pvalue", "p.adjust", "setSize")])

  cat("\nTop suppressed KEGG pathways in NAFLD (NES < 0):\n")
  kegg_dn <- kegg_res_df %>% filter(NES < 0) %>% arrange(p.adjust) %>% head(10)
  print(kegg_dn[, c("Description", "NES", "pvalue", "p.adjust", "setSize")])

  n_kegg_show <- min(20, nrow(kegg_res_df))

  # Dotplot — split by activation direction
  tryCatch({
    p_kegg_dot <- dotplot(gsea_kegg, showCategory = n_kegg_show, split = ".sign") +
      facet_grid(~.sign) +
      ggtitle("KEGG GSEA — NAFLD vs Normal",
              subtitle = "GSE162694 | ranked by sign(lfc_apeglm) x -log10(pvalue_mle) | padj<0.05") +
      theme_bw(base_size = 11) +
      theme(
        plot.title    = element_text(face = "bold", size = 13),
        axis.text.y   = element_text(size = 9),
        strip.text    = element_text(face = "bold")
      )
    ggsave("plots/gsea_kegg_dotplot.png", p_kegg_dot,
           width = 14, height = max(6, n_kegg_show * 0.4 + 3), dpi = 150)
    cat("Saved: plots/gsea_kegg_dotplot.png\n")
  }, error = function(e) cat("KEGG dotplot failed:", conditionMessage(e), "\n"))

  # Ridgeplot
  tryCatch({
    p_kegg_ridge <- ridgeplot(gsea_kegg, showCategory = n_kegg_show,
                               fill = "p.adjust") +
      labs(
        title    = "KEGG GSEA Ridgeplot — NAFLD vs Normal",
        subtitle = "Core-enrichment gene log2FC distribution | padj<0.05",
        x        = "sign(lfc_apeglm) x -log10(pvalue_mle)"
      ) +
      theme_bw(base_size = 11) +
      theme(
        plot.title  = element_text(face = "bold", size = 13),
        axis.text.y = element_text(size = 9)
      )
    ggsave("plots/gsea_kegg_ridgeplot.png", p_kegg_ridge,
           width = 12, height = max(7, n_kegg_show * 0.45 + 3), dpi = 150)
    cat("Saved: plots/gsea_kegg_ridgeplot.png\n")
  }, error = function(e) cat("KEGG ridgeplot failed:", conditionMessage(e), "\n"))

} else {
  cat("No significant KEGG pathways found at padj<0.05. Skipping plots.\n")
}

# ============================================================
# STEP 4 — Hallmark GSEA (MSigDB H collection)
# ============================================================
cat("\n=== STEP 4: Hallmark GSEA (MSigDB) ===\n")

hallmark_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, entrez_gene) %>%
  dplyr::mutate(entrez_gene = as.character(entrez_gene))

cat("Hallmark gene sets available:", length(unique(hallmark_t2g$gs_name)), "\n")

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
cat("Hallmark GSEA: significant gene sets (padj<0.05) =", nrow(gsea_hallmark@result), "\n")

hallmark_res_df <- as.data.frame(gsea_hallmark)
write.csv(hallmark_res_df, "results/gsea_hallmark_results.csv", row.names = FALSE)
cat("Saved: results/gsea_hallmark_results.csv\n")

if (nrow(hallmark_res_df) > 0) {
  cat("\nTop activated Hallmark gene sets in NAFLD (NES > 0):\n")
  hall_up <- hallmark_res_df %>% filter(NES > 0) %>% arrange(p.adjust) %>% head(10)
  print(hall_up[, c("Description", "NES", "pvalue", "p.adjust", "setSize")])

  cat("\nTop suppressed Hallmark gene sets in NAFLD (NES < 0):\n")
  hall_dn <- hallmark_res_df %>% filter(NES < 0) %>% arrange(p.adjust) %>% head(10)
  print(hall_dn[, c("Description", "NES", "pvalue", "p.adjust", "setSize")])

  n_hall_show <- min(20, nrow(hallmark_res_df))

  # Strip "HALLMARK_" prefix from Description for cleaner labels
  gsea_hallmark_clean        <- gsea_hallmark
  gsea_hallmark_clean@result$Description <- sub(
    "^HALLMARK_", "",
    gsea_hallmark_clean@result$Description
  )

  # Dotplot
  tryCatch({
    p_hall_dot <- dotplot(gsea_hallmark_clean, showCategory = n_hall_show,
                          split = ".sign") +
      facet_grid(~.sign) +
      ggtitle("Hallmark GSEA — NAFLD vs Normal",
              subtitle = "MSigDB H collection | GSE162694 | padj<0.05") +
      theme_bw(base_size = 11) +
      theme(
        plot.title  = element_text(face = "bold", size = 13),
        axis.text.y = element_text(size = 9),
        strip.text  = element_text(face = "bold")
      )
    ggsave("plots/gsea_hallmark_dotplot.png", p_hall_dot,
           width = 14, height = max(6, n_hall_show * 0.4 + 3), dpi = 150)
    cat("Saved: plots/gsea_hallmark_dotplot.png\n")
  }, error = function(e) cat("Hallmark dotplot failed:", conditionMessage(e), "\n"))

  # Ridgeplot
  tryCatch({
    p_hall_ridge <- ridgeplot(gsea_hallmark_clean, showCategory = n_hall_show,
                               fill = "p.adjust") +
      labs(
        title    = "Hallmark GSEA Ridgeplot — NAFLD vs Normal",
        subtitle = "Core-enrichment gene log2FC distribution | padj<0.05",
        x        = "sign(lfc_apeglm) x -log10(pvalue_mle)"
      ) +
      theme_bw(base_size = 11) +
      theme(
        plot.title  = element_text(face = "bold", size = 13),
        axis.text.y = element_text(size = 9)
      )
    ggsave("plots/gsea_hallmark_ridgeplot.png", p_hall_ridge,
           width = 12, height = max(7, n_hall_show * 0.45 + 3), dpi = 150)
    cat("Saved: plots/gsea_hallmark_ridgeplot.png\n")
  }, error = function(e) cat("Hallmark ridgeplot failed:", conditionMessage(e), "\n"))

} else {
  cat("No significant Hallmark gene sets found at padj<0.05. Skipping plots.\n")
}

# ============================================================
# Summary
# ============================================================
cat("\n============================================================\n")
cat("GSEA COMPLETE — Week 4 Day 1 NAFLD\n")
cat("============================================================\n")
cat("KEGG significant pathways:    ", nrow(kegg_res_df), "\n")
cat("Hallmark significant gene sets:", nrow(hallmark_res_df), "\n")
cat("\nResults: results/gsea_kegg_results.csv\n")
cat("         results/gsea_hallmark_results.csv\n")
cat("Plots:   plots/gsea_kegg_dotplot.png\n")
cat("         plots/gsea_kegg_ridgeplot.png\n")
cat("         plots/gsea_hallmark_dotplot.png\n")
cat("         plots/gsea_hallmark_ridgeplot.png\n")
