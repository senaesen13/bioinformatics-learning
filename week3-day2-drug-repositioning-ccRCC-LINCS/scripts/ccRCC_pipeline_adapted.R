################################################################################
# Week 3 Day 2 — ccRCC Drug Repositioning Pipeline (Adapted for Mac)
# Based on: Li et al. (2022) EBioMedicine 78:103963
# Original workshop code by Elif's group (KCL)
#
# What this script does:
#   Part 1 — Kaplan-Meier survival analysis on 500 TCGA KIRC genes
#   Part 2 — GO term enrichment on prognostic oncogenes
#   Part 3 — Co-expression network (random walk modules)
#   Part 4 — Drug repositioning via MSigDB C2:CGP GSEA
#             (replaces CMap .gctx approach; same connectivity map logic)
#
# Run from: week3-day2-drug-repositioning-ccRCC-LINCS/
#   Rscript scripts/ccRCC_pipeline_adapted.R
################################################################################

library(survival)
library(ggplot2)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)
library(igraph)
library(matrixStats)
library(data.table)
library(reshape2)
library(ggrepel)
library(pheatmap)

# ── Paths ────────────────────────────────────────────────────────────────────
data_dir    <- "/Users/senaesen/Desktop/Drug_Repositioning/ccRCC_drug_repositioning/Example_data/"
out_plots   <- "/Users/senaesen/Desktop/bioinfo-learning/week3-day2-drug-repositioning-ccRCC-LINCS/plots/"
out_results <- "/Users/senaesen/Desktop/bioinfo-learning/week3-day2-drug-repositioning-ccRCC-LINCS/results/"
scripts_dir <- "/Users/senaesen/Desktop/bioinfo-learning/week3-day2-drug-repositioning-ccRCC-LINCS/scripts/"

source(file.path(scripts_dir, "networkson_change.R"))
source(file.path(scripts_dir, "deg_GoTerm_clusterProfiler.R"))

# Patch makeCorNet for igraph >= 2.0 (add.colnames/add.rownames were removed)
makeCorNet <- function(corMat) {
  corNet <- igraph::graph_from_adjacency_matrix(corMat, mode = "undirected",
                                                 weighted = TRUE, diag = FALSE)
  V(corNet)$name <- rownames(corMat)
  corNet <- igraph::simplify(corNet, remove.multiple = TRUE, remove.loops = TRUE)
  return(corNet)
}

cancerType <- "TCGA_KIRC"

################################################################################
# PART 1: Kaplan-Meier survival analysis
################################################################################
cat("\n=== PART 1: Kaplan-Meier survival analysis ===\n")
cat("Goal: find genes where high expression predicts worse survival in ccRCC\n\n")

# Load TCGA KIRC data: 528 patients x (7 clinical + 500 gene) columns
exp_raw <- read.table(file.path(data_dir, paste0(cancerType, "_trans_exp_TPM_1.txt")),
                      header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cat("Data loaded:", nrow(exp_raw), "patients x", ncol(exp_raw) - 7, "genes\n")
cat("Dead patients:", sum(exp_raw$Status == "dead"), "/", nrow(exp_raw), "\n")

clinical <- exp_raw[, 1:7]
expr_mat  <- exp_raw[, 8:ncol(exp_raw)]   # patients x genes
gene_list <- colnames(expr_mat)            # Ensembl IDs

# For each gene: find optimal cutoff, run log-rank test, get Cox hazard ratio
run_km_gene <- function(j) {
  EXP        <- as.numeric(expr_mat[, j])
  LivingDays <- as.numeric(clinical$LivingDays)
  DeadInd    <- clinical$Status == "dead"
  keep       <- !is.na(EXP) & !is.na(LivingDays)
  EXP <- EXP[keep]; LivingDays <- LivingDays[keep]; DeadInd <- DeadInd[keep]

  # Test each cutoff between 20th-80th percentile; pick the most significant
  cutoffs <- sort(unique(EXP))
  cutoffs <- cutoffs[cutoffs > quantile(EXP, 0.2) & cutoffs <= quantile(EXP, 0.8)]
  if (length(cutoffs) < 2) return(NULL)

  best_p <- 1; best_cutoff <- median(EXP)
  for (co in cutoffs) {
    group <- ifelse(EXP >= co, 1, 0)
    if (length(unique(group)) < 2) next
    res <- tryCatch(
      survdiff(Surv(LivingDays, DeadInd) ~ group),
      error = function(e) NULL)
    if (is.null(res)) next
    p <- pchisq(res$chisq, 1, lower.tail = FALSE)
    if (p < best_p) { best_p <- p; best_cutoff <- co }
  }

  # Cox proportional hazards — gives hazard ratio direction
  group  <- ifelse(EXP >= best_cutoff, 1, 0)
  cox    <- tryCatch(coxph(Surv(LivingDays, DeadInd) ~ group), error = function(e) NULL)
  coef   <- if (!is.null(cox)) coef(cox)[1] else NA

  data.frame(gene = gene_list[j], p_logrank = best_p,
             cutoff = best_cutoff, coef = coef, stringsAsFactors = FALSE)
}

cat("Running KM for", length(gene_list), "genes...\n")
results_list <- vector("list", length(gene_list))
for (i in seq_along(gene_list)) {
  if (i %% 100 == 0) cat("  Progress:", i, "/", length(gene_list), "\n")
  results_list[[i]] <- run_km_gene(i)
}

km_results          <- do.call(rbind, Filter(Negate(is.null), results_list))
km_results$p_adj    <- p.adjust(km_results$p_logrank, method = "BH")
km_results$symbol   <- mapIds(org.Hs.eg.db, keys = km_results$gene,
                               keytype = "ENSEMBL", column = "SYMBOL",
                               multiVals = "first")

n_sig <- sum(km_results$p_adj < 0.05, na.rm = TRUE)
cat("Significant (BH p.adj < 0.05):", n_sig, "genes\n")

write.csv(km_results, file.path(out_results, "km_survival_all_genes.csv"), row.names = FALSE)
cat("Full KM table saved: results/km_survival_all_genes.csv\n")

# ── Plot 1: Volcano — hazard ratio vs significance ──────────────────────────
km_results$sig <- !is.na(km_results$p_adj) & km_results$p_adj < 0.05 &
                  !is.na(km_results$coef) & km_results$coef > 0

label_df <- km_results %>%
  filter(sig, !is.na(symbol)) %>%
  arrange(p_logrank) %>%
  head(20)

p1 <- ggplot(km_results, aes(x = coef, y = -log10(p_logrank + 1e-10), color = sig)) +
  geom_point(alpha = 0.6, size = 1.8) +
  scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "red3"),
                     labels = c("Not significant", "Prognostic oncogene (p.adj<0.05)")) +
  geom_text_repel(data = label_df, aes(label = symbol), size = 3,
                  max.overlaps = 25, color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.4) +
  labs(title = "Prognostic genes — TCGA KIRC (ccRCC)",
       subtitle = "Red: high expression → worse survival. Labelled: top 20 by p-value.",
       x = "Cox hazard ratio", y = "-log10(log-rank p-value)", color = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(out_plots, "01_km_volcano.png"), p1, width = 8, height = 6, dpi = 150)
cat("Plot 1 saved: 01_km_volcano.png\n")

# ── Plot 2: KM curves for top prognostic genes ──────────────────────────────
top_genes <- km_results %>%
  filter(sig) %>%
  arrange(p_logrank) %>%
  head(9)

if (nrow(top_genes) > 0) {
  km_frames <- lapply(seq_len(nrow(top_genes)), function(i) {
    gene   <- top_genes$gene[i]
    sym    <- ifelse(is.na(top_genes$symbol[i]), gene, top_genes$symbol[i])
    EXP    <- as.numeric(expr_mat[, which(gene_list == gene)])
    co     <- top_genes$cutoff[i]
    df_km  <- data.frame(
      time   = as.numeric(clinical$LivingDays) / 365,
      status = as.integer(clinical$Status == "dead"),
      group  = ifelse(EXP >= co, "High", "Low")
    )
    df_km <- df_km[!is.na(df_km$time), ]
    sf    <- survfit(Surv(time, status) ~ group, data = df_km)
    strata_labels <- gsub("group=", "", names(sf$strata))
    data.frame(
      time   = sf$time,
      surv   = sf$surv,
      group  = rep(strata_labels, sf$strata),
      panel  = paste0(sym, "\np=",
                      formatC(top_genes$p_logrank[i], format = "e", digits = 1))
    )
  })

  km_all <- do.call(rbind, km_frames)
  p2 <- ggplot(km_all, aes(x = time, y = surv, color = group)) +
    geom_step(linewidth = 0.9) +
    facet_wrap(~panel, ncol = 3) +
    scale_color_manual(values = c("High" = "red3", "Low" = "steelblue")) +
    labs(title = "Kaplan-Meier curves — top prognostic genes (ccRCC)",
         x = "Time (years)", y = "Survival probability", color = "Expression") +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom", strip.text = element_text(size = 8))

  ggsave(file.path(out_plots, "02_km_top_genes.png"), p2, width = 10, height = 10, dpi = 150)
  cat("Plot 2 saved: 02_km_top_genes.png\n")
} else {
  cat("No significant genes for KM grid plot.\n")
}

################################################################################
# PART 2: GO term enrichment on prognostic oncogenes
################################################################################
cat("\n=== PART 2: GO term enrichment ===\n")
cat("Goal: understand what biological processes the ccRCC oncogenes are involved in\n\n")

prog_ensembl <- km_results %>% filter(sig) %>% pull(gene)
cat("Prognostic oncogenes:", length(prog_ensembl), "\n")

bg_symbols   <- na.omit(mapIds(org.Hs.eg.db, keys = gene_list, keytype = "ENSEMBL",
                                column = "SYMBOL", multiVals = "first"))
prog_symbols <- na.omit(mapIds(org.Hs.eg.db, keys = prog_ensembl, keytype = "ENSEMBL",
                                column = "SYMBOL", multiVals = "first"))

cat("Gene symbols converted:", length(prog_symbols), "input /",
    length(bg_symbols), "background\n")

go_result <- enrichGO(gene          = prog_symbols,
                      OrgDb         = "org.Hs.eg.db",
                      keyType       = "SYMBOL",
                      ont           = "BP",
                      universe      = bg_symbols,
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.2,
                      qvalueCutoff  = 0.2,
                      minGSSize     = 3,
                      maxGSSize     = 500)

if (!is.null(go_result) && nrow(go_result@result) > 0) {
  go_table <- deg_GoTerm_clusterProfiler(go_result)
  write.csv(go_table, file.path(out_results, "go_enrichment_prognostic_genes.csv"),
            row.names = FALSE)
  cat("GO terms found:", nrow(go_table), "— saved to results/\n")

  p3 <- barplot(go_result, showCategory = min(20, nrow(go_result@result)),
                title = "GO Biological Process — ccRCC prognostic oncogenes") +
    theme_bw(base_size = 11)
  ggsave(file.path(out_plots, "03_go_barplot.png"), p3, width = 10, height = 7, dpi = 150)
  cat("Plot 3 saved: 03_go_barplot.png\n")
} else {
  cat("No significant GO terms (too few genes or too narrow background).\n")
  cat("This is expected with only 500 background genes.\n")
}

################################################################################
# PART 3: Co-expression network
################################################################################
cat("\n=== PART 3: Co-expression network ===\n")
cat("Goal: find clusters of genes that are co-regulated across ccRCC tumours\n\n")

# Build gene x patient matrix; keep genes with mean TPM > 1
expr_t    <- t(as.matrix(sapply(expr_mat, as.numeric)))  # genes x patients
rownames(expr_t) <- gene_list
expr_filt <- expr_t[rowMeans(expr_t) > 1, ]
cat("Genes after mean TPM > 1 filter:", nrow(expr_filt), "of", nrow(expr_t), "\n")

# Spearman co-expression: keep top 1% edges
cat("Computing Spearman co-expression (this takes ~30s)...\n")
corMatrix <- makeCorTable(expr_filt, cutoff = 0.99, mode = "spearman",
                          self = FALSE, debug = FALSE)

# Build igraph network
corNet <- makeCorNet(corMatrix)
cat("Network:", vcount(corNet), "nodes,", ecount(corNet), "edges\n")

# Random walk community detection
cat("Detecting modules (random walk clustering)...\n")
moduleList <- makeModuleList(corNet, debug = TRUE)

# Filter modules: ≥5 genes, annotate by cluster transitivity (connectivity score)
cytomat      <- annotateModulesByCC(corNet, moduleList,
                                    cutCluster = 5, cutCC = 0.4, debug = FALSE)
module_table <- cytomat$nodeTable
write.csv(module_table, file.path(out_results, "coexpression_modules.csv"), row.names = FALSE)
cat("Modules with ≥5 genes:", nrow(module_table), "— saved to results/\n")

# Add gene symbols to network nodes
node_symbols <- mapIds(org.Hs.eg.db, keys = V(corNet)$name,
                       keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first")
V(corNet)$symbol <- ifelse(is.na(node_symbols), V(corNet)$name, node_symbols)
V(corNet)$module <- membership(cluster_walktrap(corNet))

# ── Plot 4: Network coloured by module ──────────────────────────────────────
set.seed(42)
n_modules <- max(V(corNet)$module)
pal       <- rainbow(n_modules, alpha = 0.8)

png(file.path(out_plots, "04_coexpression_network.png"), width = 1400, height = 1100, res = 130)
plot(corNet,
     vertex.color = pal[V(corNet)$module],
     vertex.size  = 5,
     vertex.label = NA,
     edge.width   = 0.3,
     edge.color   = "grey85",
     main = "Co-expression network — TCGA KIRC (top 1% Spearman, random walk modules)")
dev.off()
cat("Plot 4 saved: 04_coexpression_network.png\n")

# ── Plot 5: Module size bar chart ────────────────────────────────────────────
mod_sizes <- sort(table(V(corNet)$module), decreasing = TRUE)
top_n     <- min(20, length(mod_sizes))
png(file.path(out_plots, "05_module_sizes.png"), width = 900, height = 500, res = 130)
barplot(mod_sizes[1:top_n],
        main  = paste("Top", top_n, "co-expression module sizes"),
        xlab  = "Module ID", ylab = "Number of genes",
        col   = pal[as.integer(names(mod_sizes[1:top_n]))],
        las   = 2, cex.names = 0.8)
dev.off()
cat("Plot 5 saved: 05_module_sizes.png\n")

################################################################################
# PART 4: Drug repositioning via MSigDB C2:CGP GSEA
################################################################################
cat("\n=== PART 4: Drug repositioning (MSigDB C2:CGP GSEA) ===\n")
cat("Note: Original pipeline uses CMap .gctx files (HA1E_sh.gctx, HA1E_cp.gctx)\n")
cat("      These are GB-scale and not in the example data folder.\n")
cat("      Using MSigDB C2:CGP drug gene sets as equivalent proxy:\n")
cat("      same connectivity map logic — find drugs that reverse the ccRCC signature.\n\n")

# Rank all genes by their prognostic score:
#   score = -log10(p) × sign(coef)
#   Positive = high expression → worse survival (oncogene-like)
#   Negative = low expression → worse survival (tumour-suppressor-like)
# Drugs with NES < 0 have their "downregulated" set enriched at the top
# → the drug suppresses ccRCC oncogenes → repositioning candidate
km_ranked <- km_results %>%
  filter(!is.na(coef), !is.na(p_logrank)) %>%
  mutate(score = -log10(p_logrank + 1e-10) * sign(coef)) %>%
  arrange(desc(score))

km_ranked$entrez <- mapIds(org.Hs.eg.db, keys = km_ranked$gene,
                           keytype = "ENSEMBL", column = "ENTREZID",
                           multiVals = "first")
km_ranked <- km_ranked %>% filter(!is.na(entrez))

ranked_vec <- setNames(km_ranked$score, km_ranked$entrez)
ranked_vec <- sort(ranked_vec, decreasing = TRUE)
ranked_vec <- ranked_vec[!duplicated(names(ranked_vec))]
cat("Ranked gene list:", length(ranked_vec), "genes\n")

# Load MSigDB C2:CGP; filter to PubChem-verified drug sets
cgp_all   <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CGP")
cgp_drugs <- cgp_all %>%
  filter(grepl("[PubChem=", gs_description, fixed = TRUE)) %>%
  dplyr::select(gs_name, entrez_gene = ncbi_gene)
cgp_drugs$entrez_gene <- as.character(cgp_drugs$entrez_gene)
cat("PubChem-verified drug gene sets:", length(unique(cgp_drugs$gs_name)), "\n")

# Run GSEA
set.seed(42)
gsea_res <- GSEA(geneList     = ranked_vec,
                 TERM2GENE    = cgp_drugs,
                 minGSSize    = 5,
                 maxGSSize    = 500,
                 pvalueCutoff = 1,
                 eps          = 1e-10,
                 seed         = 42,
                 verbose      = FALSE)

gsea_df <- as.data.frame(gsea_res)

# Extract drug name: last word before [PubChem=...]
extract_drug <- function(desc) {
  tolower(trimws(gsub(".*\\s(\\S+)\\s*\\[PubChem=.*", "\\1", desc)))
}
gsea_df$drug_name <- extract_drug(gsea_df$Description)

# Drug candidates: NES < 0 = drug anti-correlates with ccRCC oncogenic signature
candidates <- gsea_df %>%
  filter(p.adjust < 0.05, NES < 0) %>%
  arrange(NES) %>%
  dplyr::select(drug_name, ID, NES, pvalue, p.adjust, setSize)

write.csv(gsea_df,    file.path(out_results, "gsea_all_drugs.csv"),       row.names = FALSE)
write.csv(candidates, file.path(out_results, "gsea_drug_candidates.csv"), row.names = FALSE)

cat("\nDrug candidates (p.adj < 0.05, NES < 0):", nrow(candidates), "\n")
if (nrow(candidates) > 0) {
  print(candidates[, c("drug_name", "NES", "p.adjust")])
} else {
  cat("No hits at p.adj < 0.05. This is expected with only 500 ranked genes —\n")
  cat("GSEA power improves greatly with genome-wide ranked lists (~20,000 genes).\n")
  cat("Top 10 drug sets by NES (most negative = best candidates):\n")
  top10 <- gsea_df %>% filter(!is.na(NES)) %>% arrange(NES) %>%
    head(10) %>% dplyr::select(drug_name, NES, pvalue, setSize)
  print(top10)
  write.csv(top10, file.path(out_results, "gsea_top10_candidates_nominal.csv"), row.names=FALSE)
}

# ── Plot 6: Drug candidate lollipop ─────────────────────────────────────────
n_show  <- 12
top_neg <- gsea_df %>% filter(!is.na(NES)) %>% arrange(NES) %>% head(n_show)
top_pos <- gsea_df %>% filter(!is.na(NES)) %>% arrange(desc(NES)) %>% head(n_show)
plot_df <- bind_rows(top_neg, top_pos) %>%
  mutate(drug_label = factor(drug_name, levels = drug_name[order(NES)]),
         direction  = ifelse(NES < 0,
                             "Anti-correlated (candidate)",
                             "Co-correlated (contra-indicated)"))

p6 <- ggplot(plot_df, aes(x = NES, y = drug_label, color = direction)) +
  geom_segment(aes(xend = 0, yend = drug_label), linewidth = 0.5) +
  geom_point(size = 3) +
  scale_color_manual(values = c("Anti-correlated (candidate)"    = "steelblue",
                                "Co-correlated (contra-indicated)" = "tomato")) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  labs(title    = "Drug repositioning candidates — ccRCC",
       subtitle = "NES < 0: drug reverses ccRCC oncogenic signature → potential treatment",
       x = "Normalised Enrichment Score (NES)", y = "Drug", color = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

ggsave(file.path(out_plots, "06_drug_candidates_lollipop.png"), p6,
       width = 9, height = 7, dpi = 150)
cat("Plot 6 saved: 06_drug_candidates_lollipop.png\n")

################################################################################
# BONUS: LINCS metadata check — which of our target genes have shRNA profiles?
################################################################################
cat("\n=== BONUS: LINCS metadata ===\n")
cat("Checking siginfo_beta_yuan.Rdata for shRNA profiles of our top genes...\n")

load(file.path(data_dir, "siginfo_beta_yuan.Rdata"))  # loads sig_info

top_symbols <- km_results %>%
  filter(sig, !is.na(symbol)) %>%
  arrange(p_logrank) %>%
  head(15) %>%
  pull(symbol)

lincs_check <- as.data.frame(sig_info) %>%
  filter(pert_type %in% c("trt_sh", "trt_sh.cgs"),
         cmap_name  %in% top_symbols) %>%
  group_by(cmap_name, cell_iname) %>%
  summarise(n_signatures = dplyr::n(), .groups = "drop") %>%
  arrange(cmap_name)

if (nrow(lincs_check) > 0) {
  write.csv(lincs_check, file.path(out_results, "lincs_shrna_availability.csv"),
            row.names = FALSE)
  cat("LINCS shRNA signatures found for our target genes:\n")
  print(lincs_check)
  cat("\nThese genes have CMap profiles — full LINCS pipeline possible with .gctx files\n")
} else {
  cat("No LINCS shRNA signatures in metadata for our top prognostic genes.\n")
  cat("Target genes:", paste(top_symbols, collapse=", "), "\n")
  cat("Full pipeline would require downloading shRNA .gctx files from https://clue.io\n")
}

cat("\n==============================\n")
cat("Pipeline complete.\n")
cat("Plots:   plots/01_km_volcano.png\n")
cat("         plots/02_km_top_genes.png\n")
cat("         plots/03_go_barplot.png\n")
cat("         plots/04_coexpression_network.png\n")
cat("         plots/05_module_sizes.png\n")
cat("         plots/06_drug_candidates_lollipop.png\n")
cat("Results: results/km_survival_all_genes.csv\n")
cat("         results/go_enrichment_prognostic_genes.csv\n")
cat("         results/coexpression_modules.csv\n")
cat("         results/gsea_drug_candidates.csv\n")
cat("==============================\n")
