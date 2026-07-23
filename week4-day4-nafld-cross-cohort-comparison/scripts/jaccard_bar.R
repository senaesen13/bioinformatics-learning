## ============================================================
## Jaccard Index bar chart — three pairwise NAFLD cohort comparisons
## Jaccard = |A ∩ B| / (|A| + |B| - |A ∩ B|)
## Overlap counts sourced from pairwise_comparison.R output
## ============================================================

suppressPackageStartupMessages(library(ggplot2))

if (interactive()) {
  tryCatch({
    proj_dir <- dirname(dirname(rstudioapi::getSourceEditorContext()$path))
    setwd(proj_dir)
  }, error = function(e) invisible(NULL))
}
cat("Working directory:", getwd(), "\n")
dir.create("results", showWarnings = FALSE)

# Overlap counts from pairwise_comparison analysis
# sig threshold: padj < 0.05 AND |lfc_mle| > 1, matched by gene_symbol
jaccard_df <- data.frame(
  comparison = c("A: D3 vs D1", "B: D3 vs D2", "C: D1 vs D2"),
  sig_a      = c(180,  180,  485),
  sig_b      = c(485, 1058, 1058),
  overlap    = c( 38,   52,  139),
  stringsAsFactors = FALSE
)
jaccard_df$union   <- jaccard_df$sig_a + jaccard_df$sig_b - jaccard_df$overlap
jaccard_df$jaccard <- jaccard_df$overlap / jaccard_df$union

cat("\nJaccard Index values:\n")
print(jaccard_df[, c("comparison", "sig_a", "sig_b", "overlap", "union", "jaccard")])

jaccard_df$comparison <- factor(jaccard_df$comparison,
  levels = c("A: D3 vs D1", "B: D3 vs D2", "C: D1 vs D2"))

p <- ggplot(jaccard_df, aes(x = comparison, y = jaccard, fill = comparison)) +
  geom_col(width = 0.55, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.3f", jaccard)),
            vjust = -0.45, size = 4, fontface = "bold") +
  scale_fill_manual(values = c(
    "A: D3 vs D1" = "#1565C0",
    "B: D3 vs D2" = "#00838F",
    "C: D1 vs D2" = "#6A1B9A"
  )) +
  scale_y_continuous(limits = c(0, 0.125), expand = c(0, 0)) +
  labs(
    title    = "Jaccard Index: Pairwise NAFLD Cohort Overlap",
    subtitle = "sig genes: padj < 0.05 & |log2FC| > 1  |  Jaccard = overlap / (|A| + |B| − overlap)",
    x = NULL,
    y = "Jaccard Index"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title  = element_text(face = "bold"),
    axis.text.x = element_text(size = 11)
  )

ggsave("results/jaccard_bar.png", p, width = 6.5, height = 5, dpi = 150)
cat("Saved: results/jaccard_bar.png\n")
