# =============================================================================
# Week 2 — Day 2 — Gene Set Enrichment Analysis with clusterProfiler
# =============================================================================
# Builds on the Week 2 Day 1 DESeq2 results (airway dataset).
# Re-runs a minimal DESeq2 pipeline to get fresh results, then performs:
#   - ORA  (Over-Representation Analysis) for GO terms and KEGG pathways
#   - GSEA (Gene Set Enrichment Analysis) for GO terms and KEGG pathways
# See NOTES.md for a full explanation of every step.
# =============================================================================


# -----------------------------------------------------------------------------
# 1. INSTALL AND LOAD LIBRARIES
# -----------------------------------------------------------------------------

cran_pkgs <- c("ggplot2", "dplyr")
missing_cran <- cran_pkgs[!cran_pkgs %in% rownames(installed.packages())]
if (length(missing_cran) > 0) {
  install.packages(missing_cran)
}

bioc_pkgs <- c("DESeq2", "airway", "apeglm",
               "clusterProfiler", "org.Hs.eg.db", "enrichplot")
missing_bioc <- bioc_pkgs[!bioc_pkgs %in% rownames(installed.packages())]
if (length(missing_bioc) > 0) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install(missing_bioc, ask = FALSE)
}

library(DESeq2)
library(airway)
library(apeglm)
library(clusterProfiler)
library(org.Hs.eg.db)   # human gene annotation database
library(enrichplot)
library(ggplot2)
library(dplyr)

dir.create("plots",   showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)


# -----------------------------------------------------------------------------
# 2. RE-RUN MINIMAL DESeq2 TO GET RESULTS
# -----------------------------------------------------------------------------

data("airway")
dds <- DESeqDataSet(airway, design = ~ cell + dex)
dds$dex <- relevel(dds$dex, ref = "untrt")
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)

res        <- results(dds, name = "dex_trt_vs_untrt", alpha = 0.05)
res_shrunk <- lfcShrink(dds, coef = "dex_trt_vs_untrt", type = "apeglm")

# lfcShrink drops the stat column, so merge it back from the raw results
res_df <- as.data.frame(res_shrunk) %>%
  tibble::rownames_to_column("ensembl_id") %>%
  filter(!is.na(padj)) %>%
  left_join(
    as.data.frame(res) %>%
      tibble::rownames_to_column("ensembl_id") %>%
      dplyr::select(ensembl_id, wald_stat = stat),
    by = "ensembl_id"
  )

cat("DESeq2 done:", nrow(res_df), "genes tested\n")


# -----------------------------------------------------------------------------
# 3. ID CONVERSION: Ensembl → Entrez
# -----------------------------------------------------------------------------
#
# The airway dataset uses Ensembl gene IDs (e.g. ENSG00000000003).
# clusterProfiler's enrichGO() accepts Ensembl IDs directly.
# enrichKEGG() and GSEA functions require Entrez IDs (numeric, e.g. 7105).
# We convert using the org.Hs.eg.db annotation package.

id_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = res_df$ensembl_id,
  keytype = "ENSEMBL",
  columns = c("ENTREZID", "SYMBOL")
)

# Keep only 1:1 mappings (some Ensembl IDs map to multiple Entrez IDs)
id_map <- id_map[!duplicated(id_map$ENSEMBL) & !is.na(id_map$ENTREZID), ]

# Merge back onto results
res_df <- res_df %>%
  left_join(id_map, by = c("ensembl_id" = "ENSEMBL"))

cat("Genes with Entrez ID:", sum(!is.na(res_df$ENTREZID)), "\n")


# -----------------------------------------------------------------------------
# 4. DEFINE SIGNIFICANT GENE LISTS
# -----------------------------------------------------------------------------

# For ORA: we need a list of "significant" genes and a "universe" (background).
# Universe = all genes we tested. Significant = our hits.

sig_ensembl <- res_df %>%
  filter(padj < 0.05, abs(log2FoldChange) > 1) %>%
  pull(ensembl_id)

sig_entrez <- res_df %>%
  filter(padj < 0.05, abs(log2FoldChange) > 1, !is.na(ENTREZID)) %>%
  pull(ENTREZID)

universe_ensembl <- res_df$ensembl_id
universe_entrez  <- res_df %>% filter(!is.na(ENTREZID)) %>% pull(ENTREZID)

cat("Significant genes for ORA:", length(sig_ensembl), "\n")

# For GSEA: rank ALL genes by their DESeq2 Wald statistic.
# The stat column (LFC / SE) captures both direction and significance.
gsea_ranks <- res_df %>%
  filter(!is.na(ENTREZID), !is.na(wald_stat)) %>%
  arrange(desc(wald_stat)) %>%
  { setNames(.$wald_stat, .$ENTREZID) }

cat("Genes in GSEA ranked list:", length(gsea_ranks), "\n")


# =============================================================================
# PART A — OVER-REPRESENTATION ANALYSIS (ORA)
# =============================================================================
#
# ORA asks: "Are genes from a known pathway over-represented in my hit list?"
# It uses a hypergeometric test (like a one-sided Fisher's exact test).
#
# Inputs:
#   - Significant genes (the "hit list")
#   - Universe (all genes we could have detected)
#   - A gene set database (GO or KEGG)
#
# Limitation: treats all significant genes equally — ignores fold change magnitude.


# -----------------------------------------------------------------------------
# 5. ORA — Gene Ontology (GO)
# -----------------------------------------------------------------------------
#
# GO terms are organised into three ontologies:
#   BP = Biological Process (e.g. "inflammatory response")
#   MF = Molecular Function (e.g. "cytokine receptor binding")
#   CC = Cellular Component (e.g. "extracellular matrix")

ego <- enrichGO(
  gene          = sig_ensembl,
  universe      = universe_ensembl,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENSEMBL",
  ont           = "BP",          # Biological Process; use "ALL" for all three
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE           # convert IDs to gene symbols in output
)

cat("\nGO-BP enrichment results:", nrow(as.data.frame(ego)), "terms\n")
head(as.data.frame(ego)[, c("Description", "GeneRatio", "BgRatio", "p.adjust")])

# Simplify: remove redundant GO terms (parent/child terms covering same genes)
ego_simplified <- simplify(ego, cutoff = 0.7, by = "p.adjust")

# Dot plot: size = gene count, colour = adjusted p-value
dotplot(ego_simplified, showCategory = 20, title = "GO Biological Process — ORA") +
  theme(axis.text.y = element_text(size = 9))
ggsave("plots/go_ora_dotplot.png", width = 9, height = 7, dpi = 150)

# Bar plot
barplot(ego_simplified, showCategory = 20, title = "GO Biological Process — ORA")
ggsave("plots/go_ora_barplot.png", width = 9, height = 7, dpi = 150)

# Gene-concept network: links genes to the GO terms they drive
cnetplot(ego_simplified,
         showCategory = 8,
         foldChange   = setNames(res_df$log2FoldChange, res_df$ensembl_id))
ggsave("plots/go_ora_cnetplot.png", width = 10, height = 8, dpi = 150)


# -----------------------------------------------------------------------------
# 6. ORA — KEGG Pathways
# -----------------------------------------------------------------------------
#
# KEGG pathways are manually curated maps of known biological pathways.
# They require Entrez IDs and an organism code ("hsa" = Homo sapiens).

ekegg <- enrichKEGG(
  gene          = sig_entrez,
  universe      = universe_entrez,
  organism      = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2
)

cat("\nKEGG enrichment results:", nrow(as.data.frame(ekegg)), "pathways\n")
head(as.data.frame(ekegg)[, c("Description", "GeneRatio", "BgRatio", "p.adjust")])

if (nrow(as.data.frame(ekegg)) > 0) {
  dotplot(ekegg, showCategory = 20, title = "KEGG Pathways — ORA")
  ggsave("plots/kegg_ora_dotplot.png", width = 9, height = 6, dpi = 150)
}


# =============================================================================
# PART B — GENE SET ENRICHMENT ANALYSIS (GSEA)
# =============================================================================
#
# GSEA asks: "Do genes belonging to a pathway tend to cluster at the top
# (upregulated) or bottom (downregulated) of a ranked gene list?"
#
# Unlike ORA:
#   - Uses ALL genes, not just a significant cutoff
#   - Captures modest but consistent changes across a pathway
#   - Preserves directionality (up vs down)
#
# Key output: NES (Normalised Enrichment Score)
#   Positive NES → pathway genes enriched at top (upregulated in treated)
#   Negative NES → pathway genes enriched at bottom (downregulated in treated)


# -----------------------------------------------------------------------------
# 7. GSEA — Gene Ontology
# -----------------------------------------------------------------------------

set.seed(42)   # GSEA uses permutations; set seed for reproducibility

gsea_go <- gseGO(
  geneList      = gsea_ranks,     # named vector: Entrez ID → Wald stat
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  minGSSize     = 15,             # ignore pathways with fewer than 15 genes
  maxGSSize     = 500,
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  verbose       = FALSE
)

cat("\nGO-BP GSEA results:", nrow(as.data.frame(gsea_go)), "terms\n")

if (nrow(as.data.frame(gsea_go)) > 0) {
  # Dot plot coloured by NES (red = upregulated, blue = downregulated)
  dotplot(gsea_go, showCategory = 20, split = ".sign") +
    facet_grid(. ~ .sign) +
    labs(title = "GO Biological Process — GSEA")
  ggsave("plots/go_gsea_dotplot.png", width = 12, height = 7, dpi = 150)

  # Classic GSEA enrichment plot for the top enriched term
  top_term <- as.data.frame(gsea_go) %>%
    arrange(p.adjust) %>%
    slice(1) %>%
    pull(ID)

  gseaplot2(gsea_go,
            geneSetID = top_term,
            title     = paste("GSEA —", as.data.frame(gsea_go)[top_term, "Description"]))
  ggsave("plots/go_gsea_enrichmentplot.png", width = 8, height = 5, dpi = 150)
}


# -----------------------------------------------------------------------------
# 8. GSEA — KEGG Pathways
# -----------------------------------------------------------------------------

set.seed(42)

gsea_kegg <- gseKEGG(
  geneList      = gsea_ranks,
  organism      = "hsa",
  minGSSize     = 15,
  maxGSSize     = 500,
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  verbose       = FALSE
)

cat("\nKEGG GSEA results:", nrow(as.data.frame(gsea_kegg)), "pathways\n")

if (nrow(as.data.frame(gsea_kegg)) > 0) {
  dotplot(gsea_kegg, showCategory = 20, split = ".sign") +
    facet_grid(. ~ .sign) +
    labs(title = "KEGG Pathways — GSEA")
  ggsave("plots/kegg_gsea_dotplot.png", width = 12, height = 7, dpi = 150)

  # Enrichment plot for the most significant KEGG pathway
  top_kegg <- as.data.frame(gsea_kegg) %>%
    arrange(p.adjust) %>%
    slice(1) %>%
    pull(ID)

  gseaplot2(gsea_kegg,
            geneSetID = top_kegg,
            title     = paste("GSEA —", as.data.frame(gsea_kegg)[top_kegg, "Description"]))
  ggsave("plots/kegg_gsea_enrichmentplot.png", width = 8, height = 5, dpi = 150)
}


# -----------------------------------------------------------------------------
# 9. EXPORT RESULTS
# -----------------------------------------------------------------------------

write.csv(as.data.frame(ego),       "results/go_ora_results.csv",   row.names = FALSE)
write.csv(as.data.frame(ekegg),     "results/kegg_ora_results.csv", row.names = FALSE)
write.csv(as.data.frame(gsea_go),   "results/go_gsea_results.csv",  row.names = FALSE)
write.csv(as.data.frame(gsea_kegg), "results/kegg_gsea_results.csv",row.names = FALSE)

cat("\n--- Summary ---\n")
cat("GO ORA  terms:     ", nrow(as.data.frame(ego)),       "\n")
cat("KEGG ORA pathways: ", nrow(as.data.frame(ekegg)),     "\n")
cat("GO GSEA terms:     ", nrow(as.data.frame(gsea_go)),   "\n")
cat("KEGG GSEA pathways:", nrow(as.data.frame(gsea_kegg)), "\n")
cat("Done. Plots → plots/  Results → results/  (", getwd(), ")\n")
