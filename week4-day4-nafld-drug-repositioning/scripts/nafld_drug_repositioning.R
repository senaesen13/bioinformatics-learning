## ============================================================
## Week 4 Day 4 — NAFLD Drug Repositioning
## Strategy : GSEA against MSigDB C2 CGP drug/chemical perturbation
##            gene sets using the NAFLD vs Normal DEG signature
## Logic     : Drugs whose transcriptional response OPPOSES the
##             NAFLD signature → candidate repositioned drugs
## Input     : week4-day1 DESeq2 results (NAFLD vs Normal)
## Output    : ranked drug candidates, enrichment plots, NOTES
## Run from  : week4-day4-nafld-drug-repositioning/
##   Rscript scripts/nafld_drug_repositioning.R
## ============================================================

suppressPackageStartupMessages({
  library(fgsea)
  library(msigdbr)
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(ggrepel)
  library(patchwork)
  library(pheatmap)
  library(RColorBrewer)
  library(scales)
})

set.seed(42)

dir.create("plots",   showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

cat("=== Week 4 Day 4: NAFLD Drug Repositioning ===\n\n")

## ============================================================
## STEP 1 — Load DESeq2 results from Day 1
## ============================================================
cat("STEP 1: Loading Day 1 DESeq2 results...\n")

deseq_path <- "../week4-day1-nafld-bulk-rnaseq/results/deseq2_nafld_vs_normal.csv"
if (!file.exists(deseq_path)) {
  stop("DESeq2 results not found at: ", deseq_path,
       "\nRun week4-day1 analysis first.")
}

deseq <- read.csv(deseq_path, stringsAsFactors = FALSE)
cat("  Total rows:", nrow(deseq), "\n")
cat("  Columns:", paste(colnames(deseq), collapse=", "), "\n")

# Drop rows without a gene symbol or with NA stats
deseq <- deseq %>%
  filter(!is.na(gene_symbol), gene_symbol != "", !is.na(log2FoldChange), !is.na(padj))

cat("  Rows with gene symbol and stats:", nrow(deseq), "\n")

# Significant DEGs for summary
sig_degs <- deseq %>% filter(padj < 0.05, abs(log2FoldChange) > 1)
cat("  Significant DEGs (padj<0.05, |lfc|>1):", nrow(sig_degs),
    "(", sum(sig_degs$log2FoldChange > 0), "up /",
         sum(sig_degs$log2FoldChange < 0), "down )\n")

# ── Quick DEG summary table ───────────────────────────────────
top_up <- sig_degs %>%
  filter(log2FoldChange > 0) %>%
  arrange(padj) %>%
  head(10) %>%
  select(gene_symbol, log2FoldChange, padj)

top_dn <- sig_degs %>%
  filter(log2FoldChange < 0) %>%
  arrange(padj) %>%
  head(10) %>%
  select(gene_symbol, log2FoldChange, padj)

write.csv(sig_degs, "results/significant_degs.csv", row.names = FALSE)

## ============================================================
## STEP 2 — Build ranked gene list for GSEA
## ============================================================
cat("\nSTEP 2: Building ranked gene list...\n")

# Ranking metric: signed log10 adjusted p-value
# = log2FC × -log10(padj)  — preserves direction and significance
# Handles padj = 0 by clamping to smallest observed non-zero value
min_padj <- min(deseq$padj[deseq$padj > 0], na.rm = TRUE)
deseq <- deseq %>%
  mutate(
    padj_clamped = pmax(padj, min_padj),
    rank_metric  = log2FoldChange * -log10(padj_clamped)
  )

# Keep one entry per gene symbol (use highest absolute rank metric)
ranked <- deseq %>%
  group_by(gene_symbol) %>%
  slice_max(abs(rank_metric), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(desc(rank_metric))

gene_ranks <- setNames(ranked$rank_metric, ranked$gene_symbol)

cat("  Genes in ranked list:", length(gene_ranks), "\n")
cat("  Rank range:", round(min(gene_ranks), 2), "to", round(max(gene_ranks), 2), "\n")
cat("  Top 5 upregulated:", paste(names(head(gene_ranks, 5)), collapse=", "), "\n")
cat("  Top 5 downregulated:", paste(names(tail(gene_ranks, 5)), collapse=", "), "\n")

## ============================================================
## STEP 3 — Download MSigDB C2 CGP drug signatures
## ============================================================
cat("\nSTEP 3: Loading MSigDB C2 CGP drug signatures...\n")

cgp_raw <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CGP")
cat("  Total C2 CGP gene sets:", length(unique(cgp_raw$gs_name)), "\n")

# Convert to named list for fgsea
cgp_list <- split(cgp_raw$gene_symbol, cgp_raw$gs_name)
cat("  Converted to list of", length(cgp_list), "gene sets\n")

## ============================================================
## STEP 4 — Run fgsea against all C2 CGP gene sets
## ============================================================
cat("\nSTEP 4: Running fgsea (C2 CGP, all", length(cgp_list), "gene sets)...\n")

fgsea_res <- fgsea(
  pathways   = cgp_list,
  stats      = gene_ranks,
  minSize    = 15,
  maxSize    = 500,
  nPermSimple = 10000
)

fgsea_res <- as.data.frame(fgsea_res) %>%
  arrange(pval)

cat("  Gene sets tested:", nrow(fgsea_res), "\n")
cat("  Significant (padj<0.05):", sum(fgsea_res$padj < 0.05, na.rm=TRUE), "\n")
cat("  With negative NES (padj<0.05):", sum(fgsea_res$padj < 0.05 & fgsea_res$NES < 0, na.rm=TRUE), "\n")

## ============================================================
## STEP 5 — Drug repositioning logic
## ============================================================
cat("\nSTEP 5: Identifying drug repositioning candidates...\n")

# NAFLD signature: 581 up / 16 down → heavily skewed positive.
#
# Reversal logic (CMap framework):
#   TYPE A — _DN gene set, POSITIVE NES:
#     Genes that the drug DOWNREGULATES are ENRICHED among NAFLD-upregulated genes.
#     → Drug suppresses the NAFLD upregulation. ← primary candidates
#
#   TYPE B — _UP gene set, NEGATIVE NES:
#     Genes the drug UPREGULATES are ENRICHED among NAFLD-downregulated genes.
#     → Drug restores what NAFLD depletes. (few, because only 16 genes are down)
#
#   TYPE C — undirected set, NEGATIVE NES (NES < −1.5):
#     Overall chemical/drug perturbation signature opposes NAFLD.

fgsea_res <- fgsea_res %>%
  mutate(
    drug_name = sub("_UP$|_DN$", "", pathway),
    direction = case_when(
      grepl("_UP$", pathway) ~ "UP",
      grepl("_DN$", pathway) ~ "DN",
      TRUE                   ~ "BOTH"
    ),
    is_candidate = (
      (direction == "DN"   & NES >  1.5 & padj < 0.05) |   # TYPE A (main)
      (direction == "UP"   & NES < -1.5 & padj < 0.05) |   # TYPE B
      (direction == "BOTH" & NES < -1.5 & padj < 0.05)     # TYPE C
    )
  )

candidates <- fgsea_res %>%
  filter(is_candidate) %>%
  arrange(desc(abs(NES)))

cat("  Drug repositioning candidates:", nrow(candidates), "\n")
cat("    Type A (_DN pos NES):",
    sum(candidates$direction=="DN" & candidates$NES>0), "\n")
cat("    Type B (_UP neg NES):",
    sum(candidates$direction=="UP" & candidates$NES<0), "\n")
cat("    Type C (BOTH neg NES):",
    sum(candidates$direction=="BOTH" & candidates$NES<0), "\n")

# ── Pair UP/DN evidence per drug name ─────────────────────────
drug_summary <- fgsea_res %>%
  filter(padj < 0.05) %>%
  select(drug_name, direction, NES, padj, size) %>%
  pivot_wider(
    id_cols     = drug_name,
    names_from  = direction,
    values_from = c(NES, padj, size),
    names_sep   = "_",
    values_fn   = first
  ) %>%
  mutate(
    # Reversal score: sum of type-A and type-B evidence
    reversal_score = dplyr::case_when(
      !is.na(NES_DN) & !is.na(NES_UP) ~  NES_DN - NES_UP,
      !is.na(NES_DN)                   ~  NES_DN,
      !is.na(NES_UP)                   ~ -NES_UP,
      !is.na(NES_BOTH)                 ~ -NES_BOTH,
      TRUE                             ~  0
    )
  ) %>%
  filter(reversal_score > 1.5) %>%
  arrange(desc(reversal_score))

# Top candidates to save
top_candidates <- fgsea_res %>%
  filter(is_candidate) %>%
  arrange(desc(abs(NES))) %>%
  head(60) %>%
  select(pathway, NES, pval, padj, size, drug_name, direction)

write.csv(fgsea_res %>% select(-leadingEdge),
          "results/fgsea_all_CGP_results.csv", row.names = FALSE)
write.csv(top_candidates,
          "results/top_drug_candidates.csv", row.names = FALSE)
write.csv(drug_summary,
          "results/drug_summary_paired.csv", row.names = FALSE)

cat("  Saved drug candidate tables\n")

cat("\n  Top 15 TYPE-A candidates (_DN sets, highest positive NES):\n")
top_typeA <- fgsea_res %>%
  filter(direction == "DN", padj < 0.05) %>%
  arrange(desc(NES)) %>%
  head(15) %>%
  select(pathway, NES, padj, size)
print(top_typeA, row.names = FALSE)

## ============================================================
## STEP 6 — Volcano / NES distribution plot
## ============================================================
cat("\nSTEP 6: Plots...\n")

plot_df <- fgsea_res %>%
  mutate(
    neg_log_p = -log10(pmax(pval, 1e-10)),
    sig       = padj < 0.05,
    colour_cat = case_when(
      is_candidate & direction == "DN" ~ "Type A: _DN set, NES>0\n(drug suppresses NAFLD genes)",
      is_candidate & direction == "UP" ~ "Type B: _UP set, NES<0\n(drug restores depleted genes)",
      is_candidate                     ~ "Type C: reversal candidate",
      sig & NES > 2.0                  ~ "Mimics NAFLD (strongly)",
      sig                              ~ "Significant (other)",
      TRUE                             ~ "Not significant"
    ),
    label = ifelse(is_candidate & abs(NES) > 2.2,
                   sub("_UP$|_DN$","",pathway), "")
  )

colour_map <- c(
  "Type A: _DN set, NES>0\n(drug suppresses NAFLD genes)" = "#1565C0",
  "Type B: _UP set, NES<0\n(drug restores depleted genes)" = "#6A1B9A",
  "Type C: reversal candidate"                             = "#00695C",
  "Mimics NAFLD (strongly)"                               = "#C62828",
  "Significant (other)"                                   = "#F9A825",
  "Not significant"                                       = "#BDBDBD"
)

p_nes_volcano <- ggplot(plot_df, aes(NES, neg_log_p, colour = colour_cat)) +
  geom_point(aes(size = sig), alpha = 0.55) +
  geom_text_repel(aes(label = label), size = 2.8, max.overlaps = 25,
                  segment.size = 0.3, segment.colour = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             colour = "grey40", linewidth = 0.5) +
  geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.5) +
  scale_colour_manual(values = colour_map, name = NULL) +
  scale_size_manual(values = c("TRUE" = 1.5, "FALSE" = 0.6), guide = "none") +
  labs(title    = "GSEA Drug Repositioning — C2 CGP vs NAFLD Signature",
       subtitle  = paste0(nrow(fgsea_res), " gene sets tested | blue = candidate drugs (NES<−1.5, FDR<0.05)"),
       x = "Normalised Enrichment Score (NES)",
       y = "-log10(p-value)") +
  theme_classic(base_size = 12) +
  theme(legend.position = "right",
        legend.text = element_text(size = 9))

ggsave("plots/01_nes_volcano.png", p_nes_volcano, width = 11, height = 7, dpi = 150)
cat("  Saved plots/01_nes_volcano.png\n")

## ============================================================
## STEP 7 — Bar chart: top negative-NES drug candidates
## ============================================================
# Type A: _DN gene sets with highest positive NES (primary candidates)
top_typeA_bar <- fgsea_res %>%
  filter(direction == "DN", padj < 0.05) %>%
  arrange(desc(NES)) %>%
  head(25) %>%
  mutate(label = sub("_DN$","", pathway))

p_bar_typeA <- ggplot(top_typeA_bar,
                      aes(x = reorder(label, NES), y = NES)) +
  geom_col(fill = "#1565C0", alpha = 0.85, width = 0.75) +
  geom_text(aes(label = sprintf("FDR=%.3f  n=%d", padj, size)),
            hjust = -0.05, size = 2.8, colour = "grey20") +
  coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
  labs(
    title    = "Top Drug Candidates — Type A (_DN gene sets, positive NES)",
    subtitle = "Drug downregulates genes that are UPREGULATED in NAFLD → reversal",
    x = NULL, y = "NES (positive = drug targets NAFLD-upregulated genes)"
  ) +
  theme_classic(base_size = 11)

ggsave("plots/02_top_drug_candidates_bar.png", p_bar_typeA,
       width = 13, height = 9, dpi = 150)
cat("  Saved plots/02_top_drug_candidates_bar.png\n")

## ============================================================
## STEP 8 — Top positive NES sets (drugs that mimic NAFLD)
## ============================================================
top_pos <- fgsea_res %>%
  filter(NES > 0, padj < 0.05) %>%
  arrange(desc(NES)) %>%
  head(20) %>%
  mutate(label = sub("_UP$|_DN$","", pathway))

p_bar_pos <- ggplot(top_pos, aes(x = reorder(label, NES),
                                  y = NES, fill = direction)) +
  geom_col(width = 0.75, alpha = 0.85) +
  scale_fill_manual(values = c(UP="#C62828", DN="#E65100", BOTH="#AD1457"),
                    name = "Gene set type") +
  coord_flip() +
  labs(title    = "Top 20 Gene Sets Mimicking NAFLD Signature (positive NES)",
       subtitle = "These drugs/perturbations resemble the NAFLD transcriptome",
       x = NULL, y = "NES") +
  theme_classic(base_size = 11) +
  theme(legend.position = "bottom")

ggsave("plots/03_nafld_mimicking_sets.png", p_bar_pos,
       width = 11, height = 7, dpi = 150)
cat("  Saved plots/03_nafld_mimicking_sets.png\n")

## ============================================================
## STEP 9 — GSEA enrichment plots for top 6 drug candidates
## ============================================================
cat("\nSTEP 9: GSEA enrichment plots for top candidates...\n")

top6_paths <- fgsea_res %>%
  filter(direction == "DN", padj < 0.05) %>%
  arrange(desc(NES)) %>%
  head(6) %>%
  pull(pathway)

if (length(top6_paths) > 0) {
  enrich_plots <- lapply(top6_paths, function(p) {
    nes_val <- round(fgsea_res$NES[fgsea_res$pathway == p], 2)
    fdr_val <- formatC(fgsea_res$padj[fgsea_res$pathway == p], format="e", digits=2)
    plotEnrichment(cgp_list[[p]], gene_ranks) +
      labs(title    = gsub("_", " ", sub("_DN$|_UP$","", p)),
           subtitle = paste0("NES=", nes_val, "  FDR=", fdr_val),
           x = "Gene rank", y = "Enrichment score") +
      theme_classic(base_size = 9) +
      theme(plot.title    = element_text(size = 8, face = "bold"),
            plot.subtitle = element_text(size = 7))
  })
  p_enrich <- wrap_plots(enrich_plots, ncol = 2)
  ggsave("plots/04_gsea_enrichment_top_candidates.png", p_enrich,
         width = 12, height = max(6, ceiling(length(top6_paths)/2) * 3.5), dpi = 150)
  cat("  Saved plots/04_gsea_enrichment_top_candidates.png\n")
} else {
  cat("  No Type-A candidates found for enrichment plots\n")
}

## ============================================================
## STEP 10 — Bubble plot: NES × size × significance
## ============================================================
bubble_df <- fgsea_res %>%
  filter(padj < 0.25) %>%
  mutate(
    neg_log_fdr = -log10(pmax(padj, 1e-10)),
    is_cand     = NES < -1.5 & padj < 0.05
  )

p_bubble <- ggplot(bubble_df, aes(NES, neg_log_fdr,
                                   size = size, colour = is_cand)) +
  geom_point(alpha = 0.60) +
  geom_text_repel(
    data = bubble_df %>% filter(is_cand) %>% slice_min(NES, n = 15),
    aes(label = sub("_UP$|_DN$","",pathway)),
    size = 2.6, max.overlaps = 20, segment.size = 0.3
  ) +
  scale_colour_manual(values = c("TRUE" = "#1565C0", "FALSE" = "#9E9E9E"),
                      labels = c("TRUE"="Candidate","FALSE"="Other"),
                      name = NULL) +
  scale_size_continuous(range = c(0.5, 5), name = "Gene set size") +
  geom_hline(yintercept = -log10(0.05), linetype="dashed", colour="grey50") +
  geom_vline(xintercept = 0, colour="grey50") +
  labs(title    = "Drug Repositioning Candidates — Bubble Plot",
       subtitle = "Size = gene set size | Blue = reversal candidates (NES<−1.5, FDR<0.05)",
       x = "NES", y = "-log10(FDR)") +
  theme_classic(base_size = 11)

ggsave("plots/05_bubble_plot.png", p_bubble, width = 12, height = 7, dpi = 150)
cat("  Saved plots/05_bubble_plot.png\n")

## ============================================================
## STEP 11 — Heatmap: leading-edge genes of top 10 candidates
## ============================================================
cat("\nSTEP 11: Leading-edge gene heatmap...\n")

top10_cand <- fgsea_res %>%
  filter(direction == "DN", padj < 0.05) %>%
  arrange(desc(NES)) %>%
  head(10)

if (nrow(top10_cand) >= 3) {
  # Collect leading-edge genes
  le_genes <- top10_cand$leadingEdge
  names(le_genes) <- top10_cand$pathway
  all_le <- unique(unlist(le_genes))

  # Build presence/absence matrix
  le_mat <- sapply(le_genes, function(g) as.integer(all_le %in% g))
  rownames(le_mat) <- all_le
  colnames(le_mat) <- sub("_UP$|_DN$","", names(le_genes))

  # Keep only genes appearing in ≥2 pathways
  le_mat <- le_mat[rowSums(le_mat) >= 2, , drop=FALSE]

  if (nrow(le_mat) >= 5) {
    # Annotate rows with DEG info
    le_deg <- ranked %>%
      filter(gene_symbol %in% rownames(le_mat)) %>%
      select(gene_symbol, log2FoldChange, padj) %>%
      mutate(direction = ifelse(log2FoldChange > 0, "Up in NAFLD","Down in NAFLD")) %>%
      column_to_rownames("gene_symbol")

    row_ann <- le_deg[rownames(le_mat), "direction", drop=FALSE]
    colnames(row_ann) <- "DEG direction"

    png("plots/06_leading_edge_heatmap.png",
        width = 1200, height = max(800, nrow(le_mat)*18 + 200), res = 150)
    pheatmap(
      le_mat,
      annotation_row  = row_ann,
      annotation_colors = list(
        "DEG direction" = c("Up in NAFLD"   = "#F44336",
                            "Down in NAFLD" = "#2196F3")
      ),
      color           = c("white","#1565C0"),
      cluster_rows    = TRUE,
      cluster_cols    = TRUE,
      fontsize_row    = 8,
      fontsize_col    = 9,
      main            = "Leading-Edge Genes Shared Across Top Drug Candidates",
      border_color    = "grey80",
      legend          = FALSE
    )
    dev.off()
    cat("  Saved plots/06_leading_edge_heatmap.png\n")
  }
}

## ============================================================
## STEP 12 — DEG volcano coloured by drug leading-edge membership
## ============================================================
# Highlight genes in leading edge of top 3 drug candidates
top3_paths <- top10_cand$pathway[1:min(3, nrow(top10_cand))]
le_top3     <- unique(unlist(top10_cand$leadingEdge[top10_cand$pathway %in% top3_paths]))

vol_df <- ranked %>%
  mutate(
    neg_log_p = -log10(pmax(padj, 1e-300)),
    in_le     = gene_symbol %in% le_top3,
    sig       = padj < 0.05 & abs(log2FoldChange) > 1,
    dot_cat   = case_when(
      in_le & log2FoldChange < 0 ~ "Leading-edge (down in NAFLD)\n→ restored by drug",
      in_le & log2FoldChange > 0 ~ "Leading-edge (up in NAFLD)\n→ suppressed by drug",
      sig   ~ "Other significant DEG",
      TRUE  ~ "Not significant"
    ),
    label = ifelse(in_le & abs(log2FoldChange) > 2, gene_symbol, "")
  )

dot_cols <- c(
  "Leading-edge (down in NAFLD)\n→ restored by drug"  = "#1565C0",
  "Leading-edge (up in NAFLD)\n→ suppressed by drug"  = "#AD1457",
  "Other significant DEG"                              = "#F9A825",
  "Not significant"                                    = "#BDBDBD"
)

p_vol_le <- ggplot(vol_df, aes(log2FoldChange, neg_log_p, colour = dot_cat)) +
  geom_point(aes(size = sig), alpha = 0.5) +
  geom_text_repel(aes(label = label), size = 2.8, max.overlaps = 30,
                  segment.size = 0.3) +
  geom_hline(yintercept = -log10(0.05), linetype="dashed", colour="grey50") +
  geom_vline(xintercept = c(-1,1), linetype="dashed", colour="grey50") +
  scale_colour_manual(values = dot_cols, name = NULL) +
  scale_size_manual(values = c("TRUE"=1.4,"FALSE"=0.5), guide="none") +
  labs(title    = "NAFLD DEG Volcano — Drug Leading-Edge Genes Highlighted",
       subtitle = paste0("Leading-edge genes from top 3 drug candidates: ",
                         paste(sub("_UP$|_DN$","",top3_paths), collapse=", ")),
       x = "log2 Fold-Change (NAFLD vs Normal)",
       y = "-log10(padj)") +
  theme_classic(base_size = 11) +
  theme(legend.position = "right",
        legend.text = element_text(size = 8))

ggsave("plots/07_deg_volcano_leading_edge.png", p_vol_le,
       width = 12, height = 7, dpi = 150)
cat("  Saved plots/07_deg_volcano_leading_edge.png\n")

## ============================================================
## FINAL SUMMARY
## ============================================================
cat("\n========================================================\n")
cat("  ANALYSIS COMPLETE — NAFLD Drug Repositioning\n")
cat("========================================================\n")
cat("  DEG input:          ", nrow(deseq), "genes\n")
cat("  Sig DEGs (in/out):  ", sum(sig_degs$log2FoldChange > 0), "up /",
                              sum(sig_degs$log2FoldChange < 0), "down\n")
cat("  Gene sets tested:   ", nrow(fgsea_res), "\n")
cat("  Sig (FDR<0.05):     ", sum(fgsea_res$padj < 0.05, na.rm=TRUE), "\n")
cat("  Type-A candidates (_DN, FDR<0.05):",
    sum(fgsea_res$direction=="DN" & !is.na(fgsea_res$padj) & fgsea_res$padj < 0.05), "\n")
cat("\n  Top 10 Type-A drug candidates (_DN sets, highest NES):\n")
fgsea_res %>%
  filter(direction == "DN", padj < 0.05) %>%
  arrange(desc(NES)) %>%
  head(10) %>%
  mutate(NES=round(NES,3), padj=formatC(padj,format="e",digits=2)) %>%
  select(pathway, NES, padj, size) %>%
  { print(as.data.frame(.), row.names=FALSE) }
cat("========================================================\n")
cat("  Plots  → plots/\n")
cat("  Tables → results/\n")
cat("========================================================\n")
