## ============================================================
## Week 4 Day 1 — Re-filter with new thresholds
## New thresholds: padj < 0.05, |MLE LFC| > 1
## Reads existing deseq2_results.csv — no DESeq2 re-run needed
## ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
})

# Run from week4-day1-nafld-bulk-rnaseq/
if (interactive()) {
  tryCatch({
    proj_dir <- dirname(dirname(rstudioapi::getSourceEditorContext()$path))
    setwd(proj_dir)
  }, error = function(e) invisible(NULL))
}
cat("Working directory:", getwd(), "\n")

# ============================================================
# Load full DESeq2 results (already saved)
# ============================================================
res_df <- read.csv("results/deseq2_results.csv", stringsAsFactors = FALSE)
cat("Loaded", nrow(res_df), "genes from deseq2_results.csv\n")

# ============================================================
# Re-filter: padj_mle < 0.05 AND |lfc_mle| > 1
# ============================================================
sig_genes <- res_df %>%
  filter(!is.na(padj_mle), !is.na(lfc_mle)) %>%
  filter(padj_mle < 0.05, abs(lfc_mle) > 1) %>%
  arrange(padj_mle)

n_up   <- sum(sig_genes$lfc_mle > 0)
n_down <- sum(sig_genes$lfc_mle < 0)
n_tot  <- nrow(sig_genes)

cat("\n=== New significant genes (padj<0.05, |MLE LFC|>1) ===\n")
cat("  Up in NAFLD:  ", n_up, "\n")
cat("  Down in NAFLD:", n_down, "\n")
cat("  Total:        ", n_tot, "\n")

write.csv(sig_genes, "results/significant_genes.csv", row.names = FALSE)
cat("Saved: results/significant_genes.csv\n")

# ============================================================
# TREM2 / SPP1 / GPNMB under new thresholds
# ============================================================
cat("\n=== TREM2 / SPP1 / GPNMB (new thresholds) ===\n")
for (g in c("TREM2", "SPP1", "GPNMB")) {
  row <- res_df %>% filter(gene_symbol == g)
  if (nrow(row) == 0) { cat(g, ": not found\n"); next }
  r <- row[1, ]
  in_sig <- !is.na(r$padj_mle) && r$padj_mle < 0.05 && abs(r$lfc_mle) > 1
  cat(sprintf("%s: MLE log2FC=%+.3f | apeglm=%+.3f | padj=%.2e | sig=%s\n",
              g, r$lfc_mle, r$lfc_apeglm, r$padj_mle, in_sig))
}

# ============================================================
# Volcano plot — new thresholds
# ============================================================
cat("\n=== Volcano plot ===\n")

highlight_genes <- c("TREM2", "SPP1", "GPNMB")

volcano_df <- res_df %>%
  filter(!is.na(padj_mle), !is.na(lfc_mle)) %>%
  mutate(
    label          = ifelse(!is.na(gene_symbol), gene_symbol, ensembl_id),
    neg_log10_padj = -log10(padj_mle + 1e-300),
    sig = case_when(
      padj_mle < 0.05 & lfc_mle >  1 ~ "Up in NAFLD",
      padj_mle < 0.05 & lfc_mle < -1 ~ "Down in NAFLD",
      TRUE ~ "NS"
    )
  )

top_sig_labels <- volcano_df %>%
  filter(sig != "NS") %>%
  slice_max(order_by = neg_log10_padj, n = 20)

highlight_df <- volcano_df %>% filter(label %in% highlight_genes)
label_df     <- bind_rows(top_sig_labels, highlight_df) %>%
  distinct(ensembl_id, .keep_all = TRUE)

p_volcano <- ggplot(volcano_df, aes(x = lfc_mle, y = neg_log10_padj, colour = sig)) +
  geom_point(alpha = 0.4, size = 0.6) +
  geom_point(data = highlight_df, size = 2.5, shape = 18, colour = "#FFD600") +
  geom_text_repel(
    data         = label_df,
    aes(label    = label),
    size         = 2.5,
    max.overlaps = 25,
    segment.size = 0.3,
    segment.alpha = 0.6,
    colour       = "black"
  ) +
  scale_colour_manual(values = c(
    "Up in NAFLD"   = "#C62828",
    "Down in NAFLD" = "#1565C0",
    "NS"            = "grey70"
  )) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed",
             colour = "black", linewidth = 0.4) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             colour = "black", linewidth = 0.4) +
  labs(
    title    = "Volcano Plot: NAFLD vs Normal Liver",
    subtitle = sprintf(
      "GSE162694 | padj<0.05 & |MLE LFC|>1 | %d up, %d down | TREM2/SPP1/GPNMB = yellow diamonds",
      n_up, n_down),
    x        = "MLE log2 Fold Change (NAFLD / Normal)",
    y        = expression(-log[10](p[adj])),
    colour   = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), legend.position = "top")

ggsave("plots/volcano.png", p_volcano, width = 9, height = 7, dpi = 150)
cat("Saved: plots/volcano.png\n")

cat("\n============================================================\n")
cat("DAY 1 RE-THRESHOLD COMPLETE\n")
cat(sprintf("  Up: %d | Down: %d | Total: %d\n", n_up, n_down, n_tot))
cat("============================================================\n")
