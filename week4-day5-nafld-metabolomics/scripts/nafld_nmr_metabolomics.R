library(dplyr)
library(limma)
library(ggplot2)
library(pheatmap)
library(ggrepel)

# ── Paths ──────────────────────────────────────────────────────────────────────
maf_file  <- "data/m_MTBLS174_hna_fld_metabolite_profiling_NMR_spectroscopy_v2_maf.tsv"
meta_file <- "data/s_MTBLS174.txt"
dir.create("plots",   showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

# ── Load metabolite abundance matrix ──────────────────────────────────────────
maf <- read.delim(maf_file, check.names = FALSE, stringsAsFactors = FALSE)

# Columns 1–18 are metadata; 19+ are sample abundances
meta_cols <- 18
sample_ids <- colnames(maf)[(meta_cols + 1):ncol(maf)]

# Extract metabolite names and abundance matrix
metabolite_names <- maf$metabolite_identification
abund <- maf[, sample_ids, drop = FALSE]
rownames(abund) <- metabolite_names
abund <- as.matrix(abund)
mode(abund) <- "numeric"

# Log2 transform (add 1 to avoid log(0))
abund_log <- log2(abund + 1)

# ── Load and clean sample metadata ────────────────────────────────────────────
meta_raw <- read.delim(meta_file, check.names = FALSE, stringsAsFactors = FALSE)

meta <- meta_raw %>%
  select(
    sample_name  = `Sample Name`,
    gender       = `Factor Value[Gender]`,
    bmi          = `Factor Value[BMI]`,
    age          = `Factor Value[Age]`,
    steatosis    = `Factor Value[Steatosis]`
  ) %>%
  filter(sample_name %in% sample_ids)

# Parse steatosis to numeric (take midpoint for ranges, strip < / > symbols)
parse_steatosis <- function(x) {
  x <- gsub("<", "0", x)           # <5 → treat as 0
  x <- gsub(">", "", x)
  x <- gsub("%", "", x)
  sapply(x, function(v) {
    parts <- as.numeric(strsplit(v, "-")[[1]])
    mean(parts, na.rm = TRUE)
  })
}

meta$steatosis_num <- parse_steatosis(meta$steatosis)

# Categorise steatosis: Low (<20%), Medium (20–40%), High (>40%)
meta$steatosis_grade <- cut(
  meta$steatosis_num,
  breaks = c(-Inf, 19.9, 40, Inf),
  labels = c("Low", "Medium", "High")
)

# Align sample order to abundance matrix columns
meta <- meta[match(sample_ids, meta$sample_name), ]
rownames(meta) <- meta$sample_name

cat("Samples:", nrow(meta), "\n")
cat("Metabolites:", nrow(abund_log), "\n")
cat("Steatosis grades:\n"); print(table(meta$steatosis_grade))

# ── 1. PCA plot ───────────────────────────────────────────────────────────────
pca <- prcomp(t(abund_log), scale. = TRUE)
pca_df <- as.data.frame(pca$x[, 1:2])
pca_df$sample <- rownames(pca_df)
pca_df <- left_join(pca_df, meta, by = c("sample" = "sample_name"))

pct_var <- round(100 * summary(pca)$importance[2, 1:2], 1)

p_pca <- ggplot(pca_df, aes(PC1, PC2, colour = steatosis_grade, label = sample)) +
  geom_point(size = 3) +
  geom_text_repel(size = 2.5, max.overlaps = 15) +
  scale_colour_manual(
    values = c(Low = "#2196F3", Medium = "#FF9800", High = "#F44336"),
    na.value = "grey60"
  ) +
  labs(
    title  = "PCA of NMR Metabolomics – NAFLD Serum Samples",
    x      = paste0("PC1 (", pct_var[1], "%)"),
    y      = paste0("PC2 (", pct_var[2], "%)"),
    colour = "Steatosis Grade"
  ) +
  theme_bw()

ggsave("plots/pca_steatosis.png", p_pca, width = 7, height = 5, dpi = 150)
cat("Saved: plots/pca_steatosis.png\n")

# ── 2. Differential metabolite analysis: High vs Low steatosis ────────────────
keep <- meta$steatosis_grade %in% c("Low", "High")
abund_sub <- abund_log[, keep]
grade_sub  <- factor(meta$steatosis_grade[keep], levels = c("Low", "High"))

design <- model.matrix(~ grade_sub)
fit    <- lmFit(abund_sub, design)
fit    <- eBayes(fit)

top <- topTable(fit, coef = 2, number = Inf, sort.by = "P")
top$metabolite <- rownames(top)
top$significant <- top$adj.P.Val < 0.05 & abs(top$logFC) > 0.5

write.csv(top, "results/differential_metabolites_high_vs_low.csv", row.names = FALSE)
cat("Saved: results/differential_metabolites_high_vs_low.csv\n")
cat("Significant metabolites (FDR<0.05, |logFC|>0.5):", sum(top$significant, na.rm = TRUE), "\n")

# ── 3. Heatmap of top 15 metabolites ─────────────────────────────────────────
top_mets <- head(top$metabolite[order(top$adj.P.Val)], 15)
heat_mat  <- abund_log[top_mets, ]
heat_mat  <- t(scale(t(heat_mat)))   # z-score across samples

anno_col <- data.frame(
  Steatosis = meta$steatosis_grade,
  row.names  = meta$sample_name
)
anno_colors <- list(
  Steatosis = c(Low = "#2196F3", Medium = "#FF9800", High = "#F44336")
)

png("plots/heatmap_top15.png", width = 900, height = 700, res = 120)
pheatmap(
  heat_mat,
  annotation_col  = anno_col,
  annotation_colors = anno_colors,
  color            = colorRampPalette(c("#3F51B5", "white", "#E91E63"))(100),
  main             = "Top 15 Differential Metabolites (High vs Low Steatosis)",
  fontsize_row     = 9,
  fontsize_col     = 8,
  show_colnames    = TRUE
)
dev.off()
cat("Saved: plots/heatmap_top15.png\n")

# ── 4. Volcano plot ───────────────────────────────────────────────────────────
top$neg_log10_p <- -log10(top$P.Value)
top$label <- ifelse(top$significant, top$metabolite, NA)

p_volcano <- ggplot(top, aes(logFC, neg_log10_p, colour = significant, label = label)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_text_repel(size = 3, na.rm = TRUE, max.overlaps = 20, colour = "black") +
  scale_colour_manual(values = c("FALSE" = "grey60", "TRUE" = "#E53935"),
                      labels = c("Not significant", "Significant")) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", colour = "grey40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey40") +
  labs(
    title   = "Volcano Plot: Metabolites High vs Low Steatosis",
    x       = "log2 Fold Change",
    y       = "-log10(P-value)",
    colour  = NULL
  ) +
  theme_bw()

ggsave("plots/volcano_plot.png", p_volcano, width = 7, height = 5, dpi = 150)
cat("Saved: plots/volcano_plot.png\n")

# ── 5. Pathway enrichment (MetaboAnalyst-style KEGG mapping via enrichment) ──
# Use significant metabolites; perform simple over-representation using
# manually curated KEGG pathway groupings for the 16 NMR metabolites detected.
pathway_map <- list(
  "Amino acid metabolism"   = c("L-alanine", "Glutamine", "Glycine", "Valine",
                                "Leucine / Isoleucine", "Threonine", "Tyrosine",
                                "Phenylalanine", "Histidine", "Lysine"),
  "Energy / TCA cycle"      = c("Citrate", "Succinate", "Fumarate", "Lactate",
                                "Pyruvate", "Formate"),
  "Lipid / choline metabolism" = c("Choline", "Phosphocholine", "GPC",
                                   "VLDL / LDL lipids", "HDL lipids",
                                   "3-Hydroxybutyrate"),
  "One-carbon metabolism"   = c("Betaine", "Creatine", "Creatinine",
                                "Dimethylglycine"),
  "Gut microbiome / SCFA"   = c("Acetate", "Propionate", "Butyrate",
                                "Trimethylamine N-oxide")
)

sig_mets <- top$metabolite[top$significant %in% TRUE]
all_mets  <- top$metabolite

enrich <- lapply(names(pathway_map), function(pw) {
  pw_mets <- pathway_map[[pw]]
  sig_in  <- sum(sig_mets %in% pw_mets)
  all_in  <- sum(all_mets %in% pw_mets)
  total   <- length(all_mets)
  sig_tot <- length(sig_mets)
  # Hypergeometric p-value
  p <- phyper(sig_in - 1, all_in, total - all_in, sig_tot, lower.tail = FALSE)
  data.frame(pathway = pw, sig_in = sig_in, all_in = all_in,
             p_value = p, stringsAsFactors = FALSE)
}) %>% bind_rows() %>% arrange(p_value)

enrich$neg_log10_p <- -log10(enrich$p_value + 1e-10)
write.csv(enrich, "results/pathway_enrichment.csv", row.names = FALSE)
cat("Saved: results/pathway_enrichment.csv\n")

p_path <- ggplot(enrich, aes(reorder(pathway, neg_log10_p), neg_log10_p,
                              fill = neg_log10_p)) +
  geom_col() +
  coord_flip() +
  scale_fill_gradient(low = "#90CAF9", high = "#1565C0") +
  labs(
    title = "Pathway Enrichment of Significant Metabolites",
    x     = NULL,
    y     = "-log10(p-value)",
    fill  = "-log10(p)"
  ) +
  theme_bw() +
  theme(legend.position = "right")

ggsave("plots/pathway_enrichment.png", p_path, width = 7, height = 4, dpi = 150)
cat("Saved: plots/pathway_enrichment.png\n")

cat("\nAll done!\n")
