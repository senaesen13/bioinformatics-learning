# =============================================================================
# Week 2 — Day 4 — Drug Repositioning via clusterProfiler + MSigDB C2:CGP
# =============================================================================
# Strategy: find drugs whose gene expression signatures are OPPOSITE to the
# MI disease signature. Uses the C2:CGP (Chemical & Genetic Perturbations)
# collection from MSigDB, which contains drug/compound signatures equivalent
# to the DSigDB / CMAP connectivity-map approach — available instantly via
# msigdbr with no large database downloads.
#
# Two complementary analyses:
#   1. GSEA  — ranked gene list (log2FC), finds drugs whose programmes
#              anti-correlate with MI across the full transcriptome
#   2. ORA   — split into top up / top down genes, finds drugs with
#              significant overlap with each direction independently
#
# Run from: week2-day4-drug-repositioning/
#   Rscript scripts/drug_repositioning.R
# =============================================================================


# -----------------------------------------------------------------------------
# STEP 1 — Load libraries
# All packages are pre-installed; nothing heavy is added here.
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(msigdbr)
  library(biomaRt)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(ggplot2)
  library(dplyr)
})

dir.create("plots",   showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

if (!interactive()) pdf(NULL)   # suppress Rplots.pdf when run via Rscript

cat("=== Drug Repositioning: clusterProfiler + MSigDB C2:CGP ===\n\n")


# -----------------------------------------------------------------------------
# STEP 2 — Load DESeq2 results from Week 2 Day 3
#
# We need every gene (not just significant ones) for GSEA, because GSEA ranks
# the entire transcriptome. The significant-gene subsets are used for ORA.
# -----------------------------------------------------------------------------

cat("STEP 2: Loading DESeq2 results...\n")

res_path <- "../week2-day3-kcl-mouse-deseq2/deseq2_results_MI_vs_sham.csv"
if (!file.exists(res_path))
  stop("Run week2-day3-kcl-mouse-deseq2/deseq2_mouse_mi.R first to generate results.")

res <- read.csv(res_path)
res <- res[!is.na(res$log2FoldChange) & !is.na(res$padj), ]

cat("  Genes loaded:", nrow(res), "\n")
cat("  Significant — Up:", sum(res$sig == "Up"),
    "| Down:", sum(res$sig == "Down"), "\n\n")


# -----------------------------------------------------------------------------
# STEP 3 — Convert mouse Ensembl IDs → human Entrez IDs via BioMart
#
# WHY: The C2:CGP gene sets use human Entrez IDs. Our DESeq2 data is mouse.
# We query Ensembl for high-confidence 1:1 human orthologues, then map
# human Ensembl IDs → Entrez IDs using org.Hs.eg.db (already installed).
#
# We do this once for all genes and reuse the mapping in both GSEA and ORA.
# -----------------------------------------------------------------------------

cat("STEP 3: Converting mouse genes → human Entrez IDs via BioMart...\n")

mirrors <- c("https://asia.ensembl.org", "https://useast.ensembl.org",
             "https://uswest.ensembl.org", "https://www.ensembl.org")

mart <- NULL
for (m in mirrors) {
  mart <- tryCatch(
    useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = m),
    error = function(e) NULL
  )
  if (!is.null(mart)) { cat("  Connected to:", m, "\n"); break }
}
if (is.null(mart)) stop("Could not connect to any Ensembl mirror.")

orthologs <- getBM(
  attributes = c("ensembl_gene_id",
                 "hsapiens_homolog_ensembl_gene",
                 "hsapiens_homolog_orthology_confidence"),
  filters    = "ensembl_gene_id",
  values     = unique(res$gene_id),
  mart       = mart
)

# Keep only high-confidence 1:1 orthologues
orthologs <- orthologs[
  orthologs$hsapiens_homolog_orthology_confidence == 1 &
  orthologs$hsapiens_homolog_ensembl_gene != "", ]

# Map human Ensembl → Entrez
entrez_vec <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys      = unique(orthologs$hsapiens_homolog_ensembl_gene),
  keytype   = "ENSEMBL",
  column    = "ENTREZID",
  multiVals = "first"
)
entrez_vec <- entrez_vec[!is.na(entrez_vec)]

orthologs <- merge(
  orthologs,
  data.frame(hsapiens_homolog_ensembl_gene = names(entrez_vec),
             entrez = as.character(entrez_vec)),
  by = "hsapiens_homolog_ensembl_gene"
)
orthologs <- orthologs[!is.na(orthologs$entrez), ]

cat("  Mouse genes queried:          ", nrow(res), "\n")
cat("  High-confidence orthologues:  ", nrow(orthologs), "\n")
cat("  With Entrez IDs:              ", length(unique(orthologs$entrez)), "\n\n")


# -----------------------------------------------------------------------------
# STEP 4 — Build the ranked gene list (for GSEA) and significant-gene sets (ORA)
#
# GSEA ranked list: every gene with a human ortholog, sorted descending by
# log2FoldChange. Genes at the top are most upregulated in MI; at the bottom
# most downregulated. We use Entrez IDs as names (required by clusterProfiler).
#
# ORA gene sets: top 200 most significant up-regulated and down-regulated genes.
# 200 gives good coverage of drug landmark genes without being too noisy.
# -----------------------------------------------------------------------------

cat("STEP 4: Building gene lists...\n")

# Merge DESeq2 results with ortholog table (keep one human gene per mouse gene)
res_human <- merge(res, orthologs[, c("ensembl_gene_id", "entrez")],
                   by.x = "gene_id", by.y = "ensembl_gene_id")

# If a mouse gene maps to multiple human Entrez IDs, keep the one with highest |LFC|
res_human <- res_human %>%
  group_by(entrez) %>%
  slice_max(order_by = abs(log2FoldChange), n = 1, with_ties = FALSE) %>%
  ungroup()

# --- GSEA ranked list ---
ranked <- setNames(res_human$log2FoldChange, res_human$entrez)
ranked <- sort(ranked, decreasing = TRUE)

# --- ORA gene sets ---
sig_up_entrez <- res_human %>%
  filter(sig == "Up") %>%
  arrange(padj) %>%
  slice_head(n = 200) %>%
  pull(entrez)

sig_dn_entrez <- res_human %>%
  filter(sig == "Down") %>%
  arrange(padj) %>%
  slice_head(n = 200) %>%
  pull(entrez)

cat("  Ranked gene list size (GSEA):", length(ranked), "\n")
cat("  ORA up genes:", length(sig_up_entrez),
    "| ORA down genes:", length(sig_dn_entrez), "\n\n")


# -----------------------------------------------------------------------------
# STEP 5 — Load C2:CGP and filter to verified drug/compound signatures
#
# C2:CGP (Chemical & Genetic Perturbations) contains 3555 human gene sets,
# mixing drug experiments with genetic knockouts, tissue comparisons, etc.
# We keep only gene sets whose MSigDB description contains a PubChem compound
# ID ("[PubChem=...]") — this is MSigDB's own annotation marking verified
# drug/compound perturbation experiments. This leaves ~430 gene sets covering
# ~120 distinct compounds, which is enough for meaningful analysis without
# polluting results with non-drug gene sets.
#
# Drug name extraction: we pull the compound name from the description text
# immediately before "[PubChem=...]" — more reliable than parsing gene set
# names, which follow author/study naming conventions.
# -----------------------------------------------------------------------------

cat("STEP 5: Loading C2:CGP drug gene sets (PubChem-verified) from msigdbr...\n")

cgp_raw <- msigdbr(species = "Homo sapiens", collection = "C2",
                   subcollection = "CGP")

cat("  Total C2:CGP gene sets:", n_distinct(cgp_raw$gs_name), "\n")

# Build gene set metadata table (one row per gene set)
cgp_meta <- cgp_raw %>%
  select(gs_name, gs_description) %>%
  distinct()

# Filter to drug/compound sets with a PubChem ID in their description
drug_meta <- cgp_meta %>%
  filter(grepl("[PubChem=", gs_description, fixed = TRUE)) %>%
  mutate(
    # Extract the compound name: the last non-space word immediately before "[PubChem="]
    drug_name = tolower(trimws(
      gsub("(\\S+)\\s*\\[PubChem=.*", "\\1",
           gsub(".*\\s(\\S+)\\s*\\[PubChem=.*", "\\1", gs_description))
    ))
  )

cat("  Drug gene sets (PubChem-verified):", nrow(drug_meta), "\n")
cat("  Unique compounds:", n_distinct(drug_meta$drug_name), "\n")

# Build TERM2GENE tables for clusterProfiler (gene set name → Entrez ID)
cgp_drug_genes <- cgp_raw %>%
  filter(gs_name %in% drug_meta$gs_name) %>%
  select(gs_name, entrez_gene = ncbi_gene) %>%
  mutate(entrez_gene = as.character(entrez_gene))

# Subsets for ORA connectivity analysis
cgp_dn_drug <- cgp_drug_genes %>% filter(grepl("_DN$", gs_name))
cgp_up_drug <- cgp_drug_genes %>% filter(grepl("_UP$", gs_name))

cat("  Drug _UP sets:", n_distinct(cgp_up_drug$gs_name),
    "| Drug _DN sets:", n_distinct(cgp_dn_drug$gs_name), "\n\n")


# -----------------------------------------------------------------------------
# STEP 6 — GSEA: find drug programs anti-correlated with MI
#
# We run GSEA against drug-only C2:CGP sets using the log2FC-ranked gene list.
# Interpretation:
#   NES < 0: drug's gene program is ANTI-correlated with MI → potential therapeutic
#   NES > 0: drug mimics MI → avoid
#
# pvalueCutoff = 1 captures everything; we filter to NES<0 + p.adjust<0.05 later.
# -----------------------------------------------------------------------------

cat("STEP 6: Running GSEA (log2FC-ranked genes vs drug gene sets)...\n")
cat("  Testing", n_distinct(cgp_drug_genes$gs_name), "drug gene sets —",
    "this takes ~30-60 seconds...\n")

t0 <- proc.time()

gsea_res <- GSEA(
  geneList     = ranked,
  TERM2GENE    = cgp_drug_genes,
  minGSSize    = 5,
  maxGSSize    = 500,
  pvalueCutoff = 1,
  eps          = 1e-10,
  seed         = 42,
  verbose      = FALSE
)

elapsed <- round((proc.time() - t0)[["elapsed"]])
cat("  Done in", elapsed, "seconds.\n")

gsea_df <- as.data.frame(gsea_res)
cat("  Drug gene sets tested:         ", nrow(gsea_df), "\n")
cat("  Significant (p.adjust < 0.05):", sum(gsea_df$p.adjust < 0.05, na.rm = TRUE), "\n\n")


# -----------------------------------------------------------------------------
# STEP 7 — ORA connectivity analysis
#
# Connectivity map logic:
#   MI-UP genes ↔ drug _DN sets: drugs that SUPPRESS what MI activates
#   MI-DOWN genes ↔ drug _UP sets: drugs that RESTORE what MI suppresses
#
# enricher() uses hypergeometric test; universe = all human genes in our data.
# -----------------------------------------------------------------------------

cat("STEP 7: ORA connectivity analysis...\n")

universe <- unique(res_human$entrez)

ora_up_vs_dn <- enricher(
  gene          = sig_up_entrez,
  TERM2GENE     = cgp_dn_drug,
  universe      = universe,
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  minGSSize     = 5,
  maxGSSize     = 500
)

ora_dn_vs_up <- enricher(
  gene          = sig_dn_entrez,
  TERM2GENE     = cgp_up_drug,
  universe      = universe,
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  minGSSize     = 5,
  maxGSSize     = 500
)

ora_up_df <- if (!is.null(ora_up_vs_dn)) as.data.frame(ora_up_vs_dn) else data.frame()
ora_dn_df <- if (!is.null(ora_dn_vs_up)) as.data.frame(ora_dn_vs_up) else data.frame()

cat("  ORA hits (MI-up vs drug-DN):", nrow(ora_up_df), "\n")
cat("  ORA hits (MI-dn vs drug-UP):", nrow(ora_dn_df), "\n\n")


# -----------------------------------------------------------------------------
# STEP 8 — Consolidate GSEA results with clean drug names
#
# Join GSEA scores back to the drug_meta table which has the drug_name
# extracted from the PubChem-annotated description. Aggregate by compound:
# mean NES across all tested gene sets, best p.adjust, number of sets.
# -----------------------------------------------------------------------------

cat("STEP 8: Consolidating results...\n")

gsea_anti_mi <- gsea_df %>%
  filter(NES < 0, !is.na(p.adjust)) %>%
  left_join(drug_meta[, c("gs_name", "drug_name")], by = c("ID" = "gs_name")) %>%
  arrange(p.adjust)

drug_gsea <- gsea_anti_mi %>%
  group_by(drug_name) %>%
  summarise(
    mean_NES    = mean(NES),
    best_padj   = min(p.adjust, na.rm = TRUE),
    n_sets      = n(),
    example_set = first(ID),
    .groups     = "drop"
  ) %>%
  arrange(mean_NES)

n_sig_sets  <- sum(gsea_anti_mi$p.adjust < 0.05, na.rm = TRUE)
n_sig_drugs <- n_distinct(drug_gsea$drug_name[drug_gsea$best_padj < 0.05])

cat("  Anti-MI drug sets (NES<0):      ", nrow(gsea_anti_mi), "\n")
cat("  Significant at p.adjust<0.05:   ", n_sig_sets, "sets /",
    n_sig_drugs, "compounds\n\n")


# -----------------------------------------------------------------------------
# STEP 9 — Plot 1: Top anti-MI drug candidates (lollipop)
# -----------------------------------------------------------------------------

cat("STEP 9: Generating plots...\n")

plot_drugs <- drug_gsea %>%
  filter(best_padj < 0.05) %>%
  head(25)

if (nrow(plot_drugs) == 0) plot_drugs <- head(drug_gsea, 25)

p1 <- ggplot(plot_drugs,
             aes(x = reorder(drug_name, mean_NES), y = mean_NES,
                 colour = -log10(best_padj + 1e-10), size = n_sets)) +
  geom_segment(aes(xend = reorder(drug_name, mean_NES), y = 0, yend = mean_NES),
               colour = "grey70", linewidth = 0.4) +
  geom_point() +
  coord_flip() +
  scale_colour_gradient(low = "#AED6F1", high = "#1A5276",
                        name = "-log10(adj. p)") +
  scale_size_continuous(name = "# gene sets", range = c(2, 6)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  labs(
    title    = "Drug Repositioning — Top Anti-MI Candidates (GSEA, MSigDB C2:CGP)",
    subtitle = "Negative NES = drug gene program opposes MI expression signature",
    x        = NULL,
    y        = "Mean NES (most negative = strongest anti-MI)"
  ) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 11))

ggsave("plots/gsea_top_drugs.png", p1, width = 11, height = 8, dpi = 150)
cat("  Saved: plots/gsea_top_drugs.png\n")


# -----------------------------------------------------------------------------
# STEP 10 — Plot 2: NES distribution across all drug gene sets
# -----------------------------------------------------------------------------

gsea_df_plot <- gsea_df %>%
  filter(!is.na(p.adjust)) %>%
  mutate(sig_label = case_when(
    p.adjust < 0.05 & NES < 0 ~ "Anti-MI (sig.)",
    p.adjust < 0.05 & NES > 0 ~ "Pro-MI (sig.)",
    TRUE                       ~ "Not significant"
  ))

p2 <- ggplot(gsea_df_plot, aes(x = NES, fill = sig_label)) +
  geom_histogram(bins = 40, colour = NA, alpha = 0.85) +
  scale_fill_manual(values = c("Anti-MI (sig.)"  = "#2980B9",
                               "Pro-MI (sig.)"   = "#C0392B",
                               "Not significant" = "grey75"),
                    name = NULL) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(
    title    = "NES Distribution — Drug Gene Sets vs MI Signature (GSEA)",
    subtitle = paste0(nrow(gsea_df_plot), " drug gene sets tested | ",
                      sum(gsea_df_plot$sig_label != "Not significant"),
                      " significant (adj. p < 0.05)"),
    x        = "Normalized Enrichment Score (NES)",
    y        = "Count"
  ) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 11))

ggsave("plots/gsea_nes_distribution.png", p2, width = 8, height = 5, dpi = 150)
cat("  Saved: plots/gsea_nes_distribution.png\n")


# -----------------------------------------------------------------------------
# STEP 11 — Plot 3: ORA connectivity bubble plot
# -----------------------------------------------------------------------------

if (nrow(ora_up_df) > 0 || nrow(ora_dn_df) > 0) {

  ora_combined <- bind_rows(
    if (nrow(ora_up_df) > 0)
      left_join(mutate(ora_up_df, analysis = "MI-up genes vs drug-DN sets"),
                drug_meta[, c("gs_name", "drug_name")],
                by = c("ID" = "gs_name")) else NULL,
    if (nrow(ora_dn_df) > 0)
      left_join(mutate(ora_dn_df, analysis = "MI-dn genes vs drug-UP sets"),
                drug_meta[, c("gs_name", "drug_name")],
                by = c("ID" = "gs_name")) else NULL
  ) %>%
    mutate(
      label = ifelse(!is.na(drug_name), drug_name, ID),
      GeneRatio_num = sapply(GeneRatio, function(x) {
        p <- strsplit(x, "/")[[1]]; as.numeric(p[1]) / as.numeric(p[2])
      })
    ) %>%
    arrange(p.adjust) %>%
    slice_head(n = 30)

  p3 <- ggplot(ora_combined,
               aes(x = GeneRatio_num, y = reorder(label, GeneRatio_num),
                   colour = p.adjust, size = Count)) +
    geom_point() +
    facet_wrap(~ analysis, scales = "free_y", ncol = 1) +
    scale_colour_gradient(low = "#E74C3C", high = "#BDC3C7", name = "adj. p") +
    scale_size_continuous(name = "Genes overlapping", range = c(2, 7)) +
    labs(
      title    = "ORA Connectivity: Drug Signatures Overlapping MI DEGs",
      subtitle = "Drug-DN vs MI-up (suppress activated genes) | Drug-UP vs MI-dn (restore suppressed genes)",
      x        = "Gene Ratio",
      y        = NULL
    ) +
    theme_bw(base_size = 9) +
    theme(plot.title  = element_text(face = "bold", size = 10),
          axis.text.y = element_text(size = 7))

  ggsave("plots/ora_connectivity.png", p3, width = 10, height = 10, dpi = 150)
  cat("  Saved: plots/ora_connectivity.png\n")

} else {
  cat("  No significant ORA hits at p.adjust < 0.05; skipping ORA plot.\n")
}

cat("\n")


# -----------------------------------------------------------------------------
# STEP 12 — Save results
# -----------------------------------------------------------------------------

cat("STEP 12: Saving results...\n")

write.csv(gsea_df,    "results/gsea_all_drug_sets.csv",    row.names = FALSE)
write.csv(gsea_anti_mi, "results/gsea_anti_mi_hits.csv",  row.names = FALSE)
write.csv(drug_gsea,  "results/gsea_drug_candidates.csv",  row.names = FALSE)
write.csv(drug_meta,  "results/drug_gene_set_metadata.csv", row.names = FALSE)

if (nrow(ora_up_df) > 0)
  write.csv(ora_up_df, "results/ora_MI_up_vs_drug_DN.csv", row.names = FALSE)
if (nrow(ora_dn_df) > 0)
  write.csv(ora_dn_df, "results/ora_MI_dn_vs_drug_UP.csv", row.names = FALSE)

cat("  Saved to results/\n\n")


# -----------------------------------------------------------------------------
# SUMMARY
# -----------------------------------------------------------------------------

cat("=== SUMMARY ===\n")
cat("Drug gene sets tested (GSEA):        ", nrow(gsea_df), "\n")
cat("Anti-MI sets (NES<0):                ", nrow(gsea_anti_mi), "\n")
cat("  Significant (p.adjust<0.05):      ", n_sig_sets, "sets /", n_sig_drugs, "compounds\n\n")

cat("Top repositioning candidates (GSEA, p.adjust < 0.05):\n")
top_sig <- head(drug_gsea[drug_gsea$best_padj < 0.05, ], 10)
if (nrow(top_sig) == 0) top_sig <- head(drug_gsea, 10)
print(data.frame(top_sig[, c("drug_name", "mean_NES", "best_padj", "n_sets")]),
      row.names = FALSE)

if (nrow(ora_up_df) > 0) {
  cat("\nTop ORA hits (MI-up genes vs drug-DN sets):\n")
  print(head(ora_up_df[, c("ID", "GeneRatio", "p.adjust", "Count")], 5),
        row.names = FALSE)
}

cat("\nOutput files:\n")
cat("  plots/gsea_top_drugs.png         — lollipop of top anti-MI compounds\n")
cat("  plots/gsea_nes_distribution.png  — NES histogram across all drug sets\n")
cat("  plots/ora_connectivity.png       — ORA connectivity bubble plot\n")
cat("  results/gsea_drug_candidates.csv — aggregated GSEA results by compound\n")
cat("  results/gsea_anti_mi_hits.csv    — all anti-MI gene set hits\n")
cat("  results/ora_MI_up_vs_drug_DN.csv — ORA: MI-up genes vs drug-DN sets\n")
