## ============================================================
## Week 4 Day 4 — Pairwise Cross-Cohort Comparison
## GSE130970 (Day 3) vs GSE162694 (Day 1, discovery)
## GSE130970 (Day 3) vs GSE135251 (Day 2, validation)
## Matching by gene symbol (datasets use different gene ID systems)
## ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
})

if (interactive()) {
  tryCatch({
    proj_dir <- dirname(dirname(rstudioapi::getSourceEditorContext()$path))
    setwd(proj_dir)
  }, error = function(e) invisible(NULL))
}
cat("Working directory:", getwd(), "\n")
dir.create("results", showWarnings = FALSE)

# Paths
d1_dir <- "../week4-day1-nafld-bulk-rnaseq"
d2_dir <- "../week4-day2-nafld-validation-rnaseq"
d3_dir <- "../week4-day3-nafld-gse130970-rnaseq"

# ============================================================
# Load all three cohort results
# ============================================================
cat("\n=== Loading cohort results ===\n")

# Day 1 (GSE162694) — already filtered significant_genes.csv at new thresholds
res_d1  <- read.csv(file.path(d1_dir, "results/deseq2_results.csv"),
                    stringsAsFactors = FALSE)
sig_d1  <- read.csv(file.path(d1_dir, "results/significant_genes.csv"),
                    stringsAsFactors = FALSE)

# Day 2 (GSE135251) — apply thresholds directly to full results
res_d2_full <- read.csv(file.path(d2_dir, "results/gse135251_results.csv"),
                        stringsAsFactors = FALSE)
sig_d2 <- res_d2_full %>%
  filter(!is.na(padj_mle), !is.na(lfc_mle)) %>%
  filter(padj_mle < 0.05, abs(lfc_mle) > 1)

# Day 3 (GSE130970)
res_d3  <- read.csv(file.path(d3_dir, "results/gse130970_results.csv"),
                    stringsAsFactors = FALSE)
sig_d3  <- read.csv(file.path(d3_dir, "results/significant_genes.csv"),
                    stringsAsFactors = FALSE)

cat("Day 1 sig genes:", nrow(sig_d1), "| total tested:", nrow(res_d1), "\n")
cat("Day 2 sig genes:", nrow(sig_d2), "| total tested:", nrow(res_d2_full), "\n")
cat("Day 3 sig genes:", nrow(sig_d3), "| total tested:", nrow(res_d3), "\n")

# Deduplicate by gene_symbol (keep lowest padj_mle per symbol)
dedup <- function(df, sym_col = "gene_symbol", padj_col = "padj_mle") {
  df %>%
    filter(!is.na(.data[[sym_col]]), .data[[sym_col]] != "") %>%
    group_by(.data[[sym_col]]) %>%
    slice_min(order_by = .data[[padj_col]], n = 1, with_ties = FALSE) %>%
    ungroup()
}

res_d1_u  <- dedup(res_d1)
res_d2_u  <- dedup(res_d2_full)
res_d3_u  <- dedup(res_d3)
sig_d1_u  <- dedup(sig_d1)
sig_d2_u  <- dedup(sig_d2)
sig_d3_u  <- dedup(sig_d3)

cat("\nAfter dedup by gene_symbol:\n")
cat("  Day 1: sig", nrow(sig_d1_u), "/ total", nrow(res_d1_u), "\n")
cat("  Day 2: sig", nrow(sig_d2_u), "/ total", nrow(res_d2_u), "\n")
cat("  Day 3: sig", nrow(sig_d3_u), "/ total", nrow(res_d3_u), "\n")

# ============================================================
# Helper function: one pairwise comparison
# ============================================================
pairwise_compare <- function(res_a, sig_a, label_a,
                             res_b, sig_b, label_b) {
  cat("\n============================================================\n")
  cat(sprintf("COMPARISON: %s vs %s\n", label_a, label_b))
  cat("============================================================\n")

  # Significant gene sets (by symbol)
  sym_a_sig <- sig_a$gene_symbol
  sym_b_sig <- sig_b$gene_symbol
  overlap    <- intersect(sym_a_sig, sym_b_sig)

  cat(sprintf("  Sig genes %s: %d\n", label_a, length(sym_a_sig)))
  cat(sprintf("  Sig genes %s: %d\n", label_b, length(sym_b_sig)))
  cat(sprintf("  Overlap:       %d\n", length(overlap)))
  cat(sprintf("  Overlap %% of %s: %.1f%%\n",
              label_a, 100 * length(overlap) / max(length(sym_a_sig), 1)))

  # Fisher's exact test
  # Universe = genes with valid gene_symbol tested in both cohorts
  universe_syms <- intersect(res_a$gene_symbol, res_b$gene_symbol)
  n_u <- length(universe_syms)
  a   <- length(overlap)
  b   <- length(setdiff(sym_a_sig, sym_b_sig))
  cv  <- length(setdiff(sym_b_sig, sym_a_sig))
  d   <- n_u - a - b - cv
  ft  <- fisher.test(matrix(c(a, b, cv, d), nrow = 2), alternative = "greater")
  cat(sprintf("  Universe (shared tested genes): %d\n", n_u))
  cat(sprintf("  Fisher's exact: OR = %.2f | p = %.3e\n", ft$estimate, ft$p.value))

  # Direction concordance on overlap genes
  ov_a <- sig_a %>% filter(gene_symbol %in% overlap) %>%
    dplyr::select(gene_symbol, lfc_a = lfc_mle) %>%
    group_by(gene_symbol) %>% slice(1) %>% ungroup()
  ov_b <- sig_b %>% filter(gene_symbol %in% overlap) %>%
    dplyr::select(gene_symbol, lfc_b = lfc_mle) %>%
    group_by(gene_symbol) %>% slice(1) %>% ungroup()
  ov_df <- inner_join(ov_a, ov_b, by = "gene_symbol")

  n_conc <- sum(sign(ov_df$lfc_a) == sign(ov_df$lfc_b))
  concordant <- 100 * n_conc / max(nrow(ov_df), 1)
  binom_res <- binom.test(n_conc, nrow(ov_df), p = 0.5, alternative = "greater")
  cat(sprintf("  Direction concordance: %.1f%% (%d/%d) | binomial p = %.3e\n",
              concordant, n_conc, nrow(ov_df), binom_res$p.value))

  # Shared-tested gene set (needed for plots)
  lfc_all <- inner_join(
    res_a %>% filter(!is.na(lfc_mle), !is.na(gene_symbol)) %>%
      group_by(gene_symbol) %>% slice(1) %>% ungroup() %>%
      dplyr::select(gene_symbol, lfc_a = lfc_mle),
    res_b %>% filter(!is.na(lfc_mle), !is.na(gene_symbol)) %>%
      group_by(gene_symbol) %>% slice(1) %>% ungroup() %>%
      dplyr::select(gene_symbol, lfc_b = lfc_mle),
    by = "gene_symbol"
  )
  cat(sprintf("  Shared tested genes: %d\n", nrow(lfc_all)))

  # TREM2 / SPP1 / GPNMB status
  cat("\n  TREM2 / SPP1 / GPNMB:\n")
  for (g in c("TREM2", "SPP1", "GPNMB")) {
    a_r  <- res_a %>% filter(gene_symbol == g)
    b_r  <- res_b %>% filter(gene_symbol == g)
    a_sig <- nrow(a_r) > 0 && !is.na(a_r$padj_mle[1]) &&
             a_r$padj_mle[1] < 0.05 && abs(a_r$lfc_mle[1]) > 1
    b_sig <- nrow(b_r) > 0 && !is.na(b_r$padj_mle[1]) &&
             b_r$padj_mle[1] < 0.05 && abs(b_r$lfc_mle[1]) > 1
    cat(sprintf("    %s: %s log2FC=%s padj=%s sig=%s | %s log2FC=%s padj=%s sig=%s\n",
                g,
                label_a,
                if(nrow(a_r)>0) sprintf("%+.2f",a_r$lfc_mle[1]) else "NA",
                if(nrow(a_r)>0) sprintf("%.2e",a_r$padj_mle[1]) else "NA",
                a_sig,
                label_b,
                if(nrow(b_r)>0) sprintf("%+.2f",b_r$lfc_mle[1]) else "NA",
                if(nrow(b_r)>0) sprintf("%.2e",b_r$padj_mle[1]) else "NA",
                b_sig))
  }

  # Return a summary list for writing the markdown
  list(
    label_a = label_a, label_b = label_b,
    n_sig_a = length(sym_a_sig), n_sig_b = length(sym_b_sig),
    overlap = length(overlap),
    overlap_pct = round(100 * length(overlap) / max(length(sym_a_sig), 1), 1),
    n_universe = n_u,
    fisher_or = round(ft$estimate, 3),
    fisher_p  = signif(ft$p.value, 4),
    concordance = round(concordant, 1),
    n_concordant = n_conc,
    n_concordant_pairs = nrow(ov_df),
    binom_p_concordance = signif(binom_res$p.value, 4),
    ov_df = ov_df,
    lfc_all = lfc_all %>%
      mutate(cat = case_when(
        gene_symbol %in% overlap  ~ "Both sig",
        gene_symbol %in% sym_a_sig ~ sprintf("%s only", label_a),
        gene_symbol %in% sym_b_sig ~ sprintf("%s only", label_b),
        TRUE ~ "Neither"
      ))
  )
}

# ============================================================
# Run comparisons
# ============================================================
res_A <- pairwise_compare(
  res_d3, sig_d3, "GSE130970",
  res_d1, sig_d1, "GSE162694"
)

res_B <- pairwise_compare(
  res_d3, sig_d3, "GSE130970",
  res_d2_u, sig_d2_u, "GSE135251"
)

# ============================================================
# Save LFC correlation plots
# ============================================================
plot_lfc_corr <- function(comp, title_str, out_path) {
  lfc_dat <- comp$lfc_all
  mks <- lfc_dat %>% filter(gene_symbol %in% c("TREM2","SPP1","GPNMB"))
  both_sig <- lfc_dat %>% filter(cat == "Both sig")

  lab_a <- comp$label_a; lab_b <- comp$label_b

  p <- ggplot(lfc_dat %>% filter(cat == "Neither"),
              aes(x = lfc_a, y = lfc_b)) +
    geom_point(alpha = 0.12, size = 0.4, colour = "grey60") +
    geom_point(data = filter(lfc_dat, cat != "Neither" & cat != "Both sig"),
               alpha = 0.5, size = 0.8, colour = "#1565C0") +
    geom_point(data = both_sig, alpha = 0.8, size = 1.6, colour = "#6A1B9A") +
    geom_point(data = mks, size = 3.5, shape = 18, colour = "#FFD600") +
    geom_text_repel(data = mks, aes(label = gene_symbol),
                    size = 3.5, colour = "black",
                    nudge_y = 0.25, segment.size = 0.3) +
    { if (nrow(both_sig) >= 3)
        geom_smooth(data = both_sig, method = "lm",
                    colour = "#6A1B9A", linewidth = 0.8, se = TRUE)
      else list() } +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "grey40") +
    labs(
      title    = title_str,
      subtitle = sprintf(
        "Shared sig genes (n=%d, purple) | Direction concordance=%.1f%% (binom.p=%.2e)\nGrey=neither | Blue=one-cohort only | Yellow=TREM2/SPP1/GPNMB",
        nrow(both_sig), comp$concordance, comp$binom_p_concordance),
      x = sprintf("log2FC MLE — %s", lab_a),
      y = sprintf("log2FC MLE — %s", lab_b)
    ) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))

  ggsave(out_path, p, width = 7.5, height = 6.5, dpi = 150)
  cat("Saved:", out_path, "\n")
}

plot_lfc_corr(res_A,
  "LFC Correlation: GSE130970 vs GSE162694 (Day 1 Discovery)",
  "results/lfc_corr_d3_vs_d1.png")

plot_lfc_corr(res_B,
  "LFC Correlation: GSE130970 vs GSE135251 (Day 2 Validation)",
  "results/lfc_corr_d3_vs_d2.png")

# ============================================================
# Save summary CSV
# ============================================================
cat("\n=== Saving summary ===\n")

summary_rows <- lapply(list(res_A, res_B), function(r) {
  data.frame(
    comparison              = sprintf("%s vs %s", r$label_a, r$label_b),
    sig_genes_a             = r$n_sig_a,
    sig_genes_b             = r$n_sig_b,
    overlap                 = r$overlap,
    overlap_pct_of_a        = r$overlap_pct,
    universe_genes          = r$n_universe,
    fisher_or               = r$fisher_or,
    fisher_p                = r$fisher_p,
    concordance_pct         = r$concordance,
    n_concordant            = r$n_concordant,
    n_overlap_pairs         = r$n_concordant_pairs,
    binom_p_concordance     = r$binom_p_concordance,
    stringsAsFactors = FALSE
  )
})

summary_df <- do.call(rbind, summary_rows)
write.csv(summary_df, "results/pairwise_summary.csv", row.names = FALSE)
cat("Saved: results/pairwise_summary.csv\n")
print(t(summary_df))

cat("\n============================================================\n")
cat("PAIRWISE COMPARISON COMPLETE\n")
cat("============================================================\n")
