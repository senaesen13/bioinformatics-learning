## ============================================================
## Week 4 Day 2 — Re-filter with new thresholds + cross-cohort validation
## New thresholds: padj < 0.05, |MLE LFC| > 1
## Reads existing CSVs — no DESeq2 re-run needed
## ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
})

if (interactive()) {
  tryCatch({
    proj_dir <- dirname(dirname(rstudioapi::getSourceEditorContext()$path))
    setwd(proj_dir)
  }, error = function(e) invisible(NULL))
}
cat("Working directory:", getwd(), "\n")

disc_dir <- "../week4-day1-nafld-bulk-rnaseq"

# ============================================================
# Load validation results
# ============================================================
res_val_df <- read.csv("results/gse135251_results.csv", stringsAsFactors = FALSE)
cat("Loaded", nrow(res_val_df), "genes from gse135251_results.csv\n")

# ============================================================
# Re-filter: padj_mle < 0.05 AND |lfc_mle| > 1
# ============================================================
sig_val <- res_val_df %>%
  filter(!is.na(padj_mle), !is.na(lfc_mle)) %>%
  filter(padj_mle < 0.05, abs(lfc_mle) > 1) %>%
  arrange(padj_mle)

n_up_val   <- sum(sig_val$lfc_mle > 0)
n_down_val <- sum(sig_val$lfc_mle < 0)
cat("\n=== Validation sig genes (padj<0.05, |MLE LFC|>1) ===\n")
cat("  Up:   ", n_up_val, "\n")
cat("  Down: ", n_down_val, "\n")
cat("  Total:", nrow(sig_val), "\n")

# ============================================================
# TREM2 / SPP1 / GPNMB
# ============================================================
cat("\n=== TREM2 / SPP1 / GPNMB (new thresholds) ===\n")
for (g in c("TREM2", "SPP1", "GPNMB")) {
  row <- res_val_df %>% filter(gene_symbol == g)
  if (nrow(row) == 0) { cat(g, ": not found\n"); next }
  r <- row[1, ]
  in_sig <- !is.na(r$padj_mle) && r$padj_mle < 0.05 && abs(r$lfc_mle) > 1
  cat(sprintf("%s: MLE log2FC=%+.3f | apeglm=%+.3f | padj=%.2e | sig=%s\n",
              g, r$lfc_mle, r$lfc_apeglm, r$padj_mle, in_sig))
}

# ============================================================
# Cross-cohort comparison
# ============================================================
cat("\n=== Cross-cohort comparison ===\n")

disc_sig <- read.csv(file.path(disc_dir, "results/significant_genes.csv"),
                     stringsAsFactors = FALSE)
res_disc <- read.csv(file.path(disc_dir, "results/deseq2_results.csv"),
                     stringsAsFactors = FALSE)

disc_ids <- disc_sig$ensembl_id
val_ids  <- sig_val$ensembl_id
overlap  <- intersect(disc_ids, val_ids)

cat("Discovery significant genes:  ", length(disc_ids), "\n")
cat("Validation significant genes: ", length(val_ids), "\n")
cat("Overlap:                      ", length(overlap), "\n")
cat("Overlap % of discovery:       ",
    round(100 * length(overlap) / length(disc_ids), 1), "%\n")

# Fisher's exact test
universe_genes <- intersect(res_disc$ensembl_id, res_val_df$ensembl_id)
a  <- length(overlap)
b  <- length(setdiff(disc_ids, val_ids))
cv <- length(setdiff(val_ids, disc_ids))
d  <- length(universe_genes) - a - b - cv
ft <- fisher.test(matrix(c(a, b, cv, d), nrow = 2), alternative = "greater")
cat(sprintf("Fisher's exact: OR = %.2f | p = %.3e\n", ft$estimate, ft$p.value))

# Direction concordance
ov_disc <- disc_sig %>%
  filter(ensembl_id %in% overlap) %>%
  select(ensembl_id, lfc_disc = lfc_mle, gene_symbol)
ov_val <- sig_val %>%
  filter(ensembl_id %in% overlap) %>%
  select(ensembl_id, lfc_val = lfc_mle)
ov_df <- inner_join(ov_disc, ov_val, by = "ensembl_id")
concordant <- mean(sign(ov_df$lfc_disc) == sign(ov_df$lfc_val)) * 100
cat(sprintf("Direction concordance: %.1f%%\n", concordant))

# Correlations (sig overlap)
cor_p <- cor.test(ov_df$lfc_disc, ov_df$lfc_val, method = "pearson")
cor_s <- cor.test(ov_df$lfc_disc, ov_df$lfc_val, method = "spearman")
cat(sprintf("Pearson  r (sig overlap):   %.4f (p=%.3e)\n",
            cor_p$estimate, cor_p$p.value))
cat(sprintf("Spearman rho (sig overlap): %.4f (p=%.3e)\n",
            cor_s$estimate, cor_s$p.value))

# Genome-wide correlation
lfc_both <- inner_join(
  res_disc   %>% filter(!is.na(lfc_mle)) %>% select(ensembl_id, lfc_disc = lfc_mle, gene_symbol),
  res_val_df %>% filter(!is.na(lfc_mle)) %>% select(ensembl_id, lfc_val  = lfc_mle),
  by = "ensembl_id"
)
cor_gw <- cor.test(lfc_both$lfc_disc, lfc_both$lfc_val, method = "pearson")
cat(sprintf("Pearson r genome-wide (%d genes): %.4f\n",
            nrow(lfc_both), cor_gw$estimate))

# Overlap gene table
cat("\n--- Overlap genes ---\n")
ov_df_full <- ov_df %>%
  mutate(concordant_dir = sign(lfc_disc) == sign(lfc_val))
print(ov_df_full %>% select(gene_symbol, lfc_disc, lfc_val, concordant_dir))

# Cross-cohort marker check
cat("\nTREM2/SPP1/GPNMB cross-cohort:\n")
for (g in c("TREM2", "SPP1", "GPNMB")) {
  d_r   <- res_disc   %>% filter(gene_symbol == g)
  v_r   <- res_val_df %>% filter(gene_symbol == g)
  d_sig <- nrow(d_r) > 0 && !is.na(d_r$padj_mle[1]) &&
           d_r$padj_mle[1] < 0.05 && abs(d_r$lfc_mle[1]) > 1
  v_sig <- nrow(v_r) > 0 && !is.na(v_r$padj_mle[1]) &&
           v_r$padj_mle[1] < 0.05 && abs(v_r$lfc_mle[1]) > 1
  cat(sprintf(
    "%s: disc log2FC=%s padj=%s sig=%s | val log2FC=%s padj=%s sig=%s\n", g,
    if (nrow(d_r) > 0) sprintf("%+.3f", d_r$lfc_mle[1]) else "NA",
    if (nrow(d_r) > 0) sprintf("%.2e",  d_r$padj_mle[1]) else "NA", d_sig,
    if (nrow(v_r) > 0) sprintf("%+.3f", v_r$lfc_mle[1]) else "NA",
    if (nrow(v_r) > 0) sprintf("%.2e",  v_r$padj_mle[1]) else "NA", v_sig
  ))
}

# ============================================================
# Overlap bar chart
# ============================================================
cat("\n=== Overlap plot ===\n")

venn_df <- data.frame(
  Category = c("Discovery only", "Both cohorts", "Validation only"),
  Count    = c(length(disc_ids) - length(overlap),
               length(overlap),
               length(val_ids)  - length(overlap)),
  Group    = c("disc", "both", "val")
)
p_venn <- ggplot(venn_df, aes(x = reorder(Category, -Count), y = Count, fill = Group)) +
  geom_col(width = 0.55) +
  geom_text(aes(label = Count), vjust = -0.5, size = 5, fontface = "bold") +
  scale_fill_manual(values = c(disc = "#1565C0", both = "#6A1B9A", val = "#C62828")) +
  labs(
    title    = "Overlap: Discovery vs Validation Significant Genes",
    subtitle = sprintf(
      "GSE162694 (n=%d) vs GSE135251 (n=%d) | padj<0.05 & |MLE LFC|>1 | %d shared",
      length(disc_ids), length(val_ids), length(overlap)),
    x = NULL, y = "Gene count"
  ) +
  theme_bw(base_size = 13) +
  theme(plot.title = element_text(face = "bold"), legend.position = "none",
        axis.text.x = element_text(size = 11))
ggsave("plots/validation_overlap.png", p_venn, width = 7, height = 5, dpi = 150)
cat("Saved: plots/validation_overlap.png\n")

# ============================================================
# LFC correlation scatter
# ============================================================
cat("\n=== LFC correlation plot ===\n")

lfc_both <- lfc_both %>%
  mutate(cat = case_when(
    ensembl_id %in% overlap  ~ "Both sig",
    ensembl_id %in% disc_ids ~ "Discovery only",
    ensembl_id %in% val_ids  ~ "Validation only",
    TRUE                     ~ "Neither"
  ))
highlight_mks <- lfc_both %>% filter(gene_symbol %in% c("TREM2", "SPP1", "GPNMB"))

p_cor <- ggplot(lfc_both %>% filter(cat == "Neither"),
                aes(x = lfc_disc, y = lfc_val)) +
  geom_point(alpha = 0.12, size = 0.4, colour = "grey60") +
  geom_point(data = filter(lfc_both, cat == "Discovery only"),
             alpha = 0.6, size = 1, colour = "#1565C0") +
  geom_point(data = filter(lfc_both, cat == "Validation only"),
             alpha = 0.6, size = 1, colour = "#C62828") +
  geom_point(data = filter(lfc_both, cat == "Both sig"),
             alpha = 0.8, size = 1.8, colour = "#6A1B9A") +
  geom_point(data = highlight_mks, size = 3.5, shape = 18, colour = "#FFD600") +
  geom_text_repel(data = highlight_mks, aes(label = gene_symbol),
                  size = 3.5, colour = "black",
                  nudge_y = 0.25, segment.size = 0.3) +
  geom_smooth(data = filter(lfc_both, cat == "Both sig"),
              method = "lm", colour = "#6A1B9A", linewidth = 0.8, se = TRUE) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "grey40") +
  labs(
    title    = "LFC Correlation: Discovery vs Validation",
    subtitle = sprintf(
      "Shared sig genes (n=%d, purple) | Pearson r=%.3f | Spearman ρ=%.3f\nGrey=neither sig | Blue=disc only | Red=val only | Yellow=TREM2/SPP1/GPNMB",
      nrow(ov_df), cor_p$estimate, cor_s$estimate),
    x = "log2FC MLE — Discovery (GSE162694)",
    y = "log2FC MLE — Validation (GSE135251)"
  ) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))
ggsave("plots/lfc_correlation.png", p_cor, width = 7.5, height = 6.5, dpi = 150)
cat("Saved: plots/lfc_correlation.png\n")

# ============================================================
# Save updated overlap summary CSV
# ============================================================
cat("\n=== Saving overlap summary CSV ===\n")

summary_df <- data.frame(
  metric = c(
    "discovery_sig_genes", "validation_sig_genes", "overlap_genes",
    "overlap_pct_of_discovery", "universe_shared_tested_genes",
    "fishers_exact_OR", "fishers_exact_p",
    "direction_concordance_pct",
    "pearson_r_sig_overlap", "pearson_p_sig_overlap",
    "spearman_rho_sig_overlap", "spearman_p_sig_overlap",
    "pearson_r_genomewide",
    "TREM2_sig_discovery", "TREM2_sig_validation",
    "SPP1_sig_discovery",  "SPP1_sig_validation",
    "GPNMB_sig_discovery", "GPNMB_sig_validation"
  ),
  value = c(
    length(disc_ids), length(val_ids), length(overlap),
    round(100 * length(overlap) / length(disc_ids), 2),
    length(universe_genes),
    round(ft$estimate, 3), signif(ft$p.value, 4),
    round(concordant, 2),
    round(cor_p$estimate, 4), signif(cor_p$p.value, 4),
    round(cor_s$estimate, 4), signif(cor_s$p.value, 4),
    round(cor_gw$estimate, 4),
    any(disc_sig$gene_symbol == "TREM2"),
    any(sig_val$gene_symbol  == "TREM2"),
    any(disc_sig$gene_symbol == "SPP1"),
    any(sig_val$gene_symbol  == "SPP1"),
    any(disc_sig$gene_symbol == "GPNMB"),
    any(sig_val$gene_symbol  == "GPNMB")
  )
)
write.csv(summary_df, "results/validation_overlap_summary.csv", row.names = FALSE)
cat("Saved: results/validation_overlap_summary.csv\n")

cat("\n============================================================\n")
cat("DAY 2 RE-THRESHOLD COMPLETE\n")
cat(sprintf("  Validation: %d up | %d down | %d total\n",
            n_up_val, n_down_val, nrow(sig_val)))
cat(sprintf("  Overlap: %d genes\n", length(overlap)))
cat(sprintf("  Fisher's: OR=%.2f  p=%.3e\n", ft$estimate, ft$p.value))
cat(sprintf("  Concordance: %.1f%% | Spearman rho=%.4f\n",
            concordant, cor_s$estimate))
cat("============================================================\n")
