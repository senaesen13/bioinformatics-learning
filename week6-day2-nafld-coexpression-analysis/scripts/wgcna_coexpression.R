## ============================================================
## WGCNA-style coexpression: 139 real NAFLD genes, GSE162694 VST
## Method: Pearson correlation, soft-thresholding power selection,
##         adjacency matrix, Topological Overlap Matrix (TOM),
##         hierarchical clustering, module eigengenes (first PC).
## ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(cluster)   # for silhouette()
})

if (interactive()) {
  tryCatch({
    proj_dir <- dirname(dirname(rstudioapi::getSourceEditorContext()$path))
    setwd(proj_dir)
  }, error = function(e) invisible(NULL))
}
cat("Working directory:", getwd(), "\n")
dir.create("results", showWarnings = FALSE)
dir.create("plots",   showWarnings = FALSE)

set.seed(42)

# ============================================================
# 1. Load real VST expression matrix (genes × samples)
# ============================================================
cat("\n=== 1. Loading VST expression matrix ===\n")
vst      <- read.csv("results/vst_matrix_139genes.csv", row.names = 1, check.names = FALSE)
expr_mat <- as.matrix(vst)
cat(sprintf("Matrix: %d genes × %d samples\n", nrow(expr_mat), ncol(expr_mat)))

# ============================================================
# 2. Pearson correlation matrix (gene × gene, across 143 samples)
# ============================================================
cat("\n=== 2. Computing correlation matrix ===\n")
cor_mat <- cor(t(expr_mat), method = "pearson")
cat(sprintf("Correlation range: [%.3f, %.3f]\n",
            min(cor_mat[lower.tri(cor_mat)]),
            max(cor_mat[lower.tri(cor_mat)])))

# ============================================================
# 3. Soft-thresholding power selection (scale-free topology fit)
# ============================================================
cat("\n=== 3. Soft-thresholding power selection ===\n")
powers  <- 1:20
sft_df  <- do.call(rbind, lapply(powers, function(beta) {
  adj <- abs(cor_mat)^beta
  diag(adj) <- 0
  k   <- rowSums(adj)

  # Connectivity distribution → scale-free R²
  h    <- hist(k, breaks = 10, plot = FALSE)
  pk   <- h$counts / sum(h$counts)
  kv   <- h$mids
  keep <- pk > 0 & kv > 0
  if (sum(keep) < 3) {
    return(data.frame(power = beta, r_sq = NA_real_, mean_k = round(mean(k), 2)))
  }
  fit  <- lm(log10(pk[keep]) ~ log10(kv[keep]))
  r_sq <- round(summary(fit)$r.squared, 4)
  data.frame(power = beta, r_sq = r_sq, mean_k = round(mean(k), 2))
}))

cat("Power selection table:\n")
print(sft_df)
write.csv(sft_df, "results/wgcna_soft_threshold.csv", row.names = FALSE)

# Select: first power where R² >= 0.80, else best available
good <- sft_df$power[!is.na(sft_df$r_sq) & sft_df$r_sq >= 0.80]
if (length(good) > 0) {
  soft_power <- min(good)
  cat(sprintf("\nSelected β = %d (first power achieving R² >= 0.80)\n", soft_power))
} else {
  soft_power <- sft_df$power[which.max(sft_df$r_sq)]
  best_r2    <- max(sft_df$r_sq, na.rm = TRUE)
  cat(sprintf(
    "\nNo power reaches R² >= 0.80 (max = %.3f); n = 139 genes limits scale-free sensitivity.\n",
    best_r2))
  cat(sprintf("Selected β = %d (best available R²).\n", soft_power))
}

# SFT plot
p_sft <- ggplot(sft_df[!is.na(sft_df$r_sq), ], aes(x = power, y = r_sq)) +
  geom_line(colour = "#1565C0", linewidth = 0.8) +
  geom_point(size = 2.5) +
  geom_vline(xintercept = soft_power, linetype = "dashed",
             colour = "red", linewidth = 0.7) +
  geom_hline(yintercept = 0.80, linetype = "dotted", colour = "grey40") +
  annotate("text", x = soft_power + 0.7, y = 0.05,
           label = sprintf("β = %d", soft_power),
           colour = "red", hjust = 0, size = 3.5) +
  scale_x_continuous(breaks = powers) +
  labs(
    title    = "Soft-Thresholding Power Selection (Scale-Free Topology Fit)",
    subtitle = sprintf("139 genes × 143 samples  |  Selected β = %d", soft_power),
    x = "Soft-thresholding power (β)",
    y = "Scale-free topology R²"
  ) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))
ggsave("plots/wgcna_soft_threshold.png", p_sft, width = 7, height = 4.5, dpi = 150)
cat("Saved: plots/wgcna_soft_threshold.png\n")

# ============================================================
# 4. Soft-thresholded adjacency matrix
# ============================================================
cat(sprintf("\n=== 4. Adjacency matrix (β = %d) ===\n", soft_power))
adj_mat <- abs(cor_mat)^soft_power
diag(adj_mat) <- 0
cat(sprintf("Mean connectivity: %.4f\n", mean(rowSums(adj_mat))))

# ============================================================
# 5. Topological Overlap Matrix (TOM)
#    TOM(i,j) = (Σ_u a_iu * a_uj + a_ij) / (min(k_i, k_j) + 1 - a_ij)
# ============================================================
cat("\n=== 5. Computing TOM ===\n")
k_conn  <- rowSums(adj_mat)
l_mat   <- adj_mat %*% adj_mat          # Σ_u a_iu * a_uj
min_k   <- outer(k_conn, k_conn, FUN = pmin)
tom     <- (l_mat + adj_mat) / (min_k + 1 - adj_mat)
diag(tom) <- 1
dist_tom  <- 1 - tom
cat(sprintf("TOM computed. dist_tom range: [%.4f, %.4f]\n",
            min(dist_tom[upper.tri(dist_tom)]),
            max(dist_tom[upper.tri(dist_tom)])))

# ============================================================
# 6. Hierarchical clustering + silhouette module selection
# ============================================================
cat("\n=== 6. Hierarchical clustering (average linkage on 1 - TOM) ===\n")
gene_tree <- hclust(as.dist(dist_tom), method = "average")

sil_k <- sapply(2:8, function(kk) {
  ct  <- cutree(gene_tree, k = kk)
  sil <- silhouette(ct, dist = as.dist(dist_tom))
  mean(sil[, 3])
})
names(sil_k) <- 2:8

cat("Silhouette by k:\n")
print(round(sil_k, 4))

k_best <- as.integer(names(which.max(sil_k)))
cat(sprintf("Selected k = %d (silhouette = %.4f)\n", k_best, max(sil_k)))

mod_cut <- cutree(gene_tree, k = k_best)

# WGCNA color convention: largest module → turquoise, then blue, brown, yellow …
wgcna_palette <- c("turquoise", "blue", "brown", "yellow", "green", "red", "black", "pink")
mod_by_size   <- order(table(mod_cut), decreasing = TRUE)
color_map     <- setNames(wgcna_palette[seq_along(mod_by_size)],
                          as.character(mod_by_size))
module_colors <- color_map[as.character(mod_cut)]
names(module_colors) <- rownames(vst)

gene_modules <- data.frame(
  gene         = rownames(vst),
  wgcna_module = paste0("ME_", module_colors),
  module_num   = mod_cut,
  kWithin      = rowSums(adj_mat),
  stringsAsFactors = FALSE
)

cat("\nModule sizes (WGCNA color names):\n")
print(sort(table(gene_modules$wgcna_module), decreasing = TRUE))

# ============================================================
# 7. Module eigengenes (first principal component per module)
# ============================================================
cat("\n=== 7. Computing module eigengenes ===\n")

mod_names_sorted <- names(sort(table(gene_modules$wgcna_module), decreasing = TRUE))

eigen_rows <- lapply(mod_names_sorted, function(mod) {
  mod_genes <- gene_modules$gene[gene_modules$wgcna_module == mod]

  if (length(mod_genes) == 1) {
    # Single-gene module: use its scaled expression as the eigengene
    ev <- as.numeric(scale(expr_mat[mod_genes, ]))
  } else {
    # Multi-gene module: first PC across samples
    x_scaled <- scale(t(expr_mat[mod_genes, , drop = FALSE]))  # samples × genes
    pca_res  <- prcomp(x_scaled, center = FALSE, scale. = FALSE)
    ev       <- pca_res$x[, 1]
    # Sign convention: positive correlation with majority of module genes
    if (mean(cor(ev, x_scaled) > 0) < 0.5) ev <- -ev
  }
  c(module = mod, setNames(as.numeric(ev), colnames(expr_mat)))
})

eigengene_df <- as.data.frame(do.call(rbind, eigen_rows), stringsAsFactors = FALSE)
# Ensure numeric columns are numeric
for (col in colnames(eigengene_df)[-1]) {
  eigengene_df[[col]] <- as.numeric(eigengene_df[[col]])
}
write.csv(eigengene_df, "results/wgcna_module_eigengenes.csv", row.names = FALSE)
cat(sprintf("Module eigengenes saved: %d modules × %d samples\n",
            nrow(eigengene_df), ncol(eigengene_df) - 1))

# ============================================================
# 8. Key gene module assignments
# ============================================================
cat("\n=== 8. Key gene module assignments ===\n")
key_genes <- c("TREM2", "SPP1", "GPNMB", "FASN", "FABP4", "CD68", "COL1A1")
key_results <- lapply(key_genes, function(g) {
  row <- gene_modules[gene_modules$gene == g, ]
  if (nrow(row) > 0) {
    cat(sprintf("  %-10s → %s  (kWithin = %.2f)\n", g, row$wgcna_module, row$kWithin))
    data.frame(gene = g, wgcna_module = row$wgcna_module, in_overlap = TRUE,
               kWithin = round(row$kWithin, 2), stringsAsFactors = FALSE)
  } else {
    cat(sprintf("  %-10s → NOT in 139-gene overlap\n", g))
    data.frame(gene = g, wgcna_module = NA, in_overlap = FALSE,
               kWithin = NA_real_, stringsAsFactors = FALSE)
  }
})
key_df <- do.call(rbind, key_results)

# ============================================================
# 9. Dendrogram + module color bar
# ============================================================
cat("\n=== 9. Generating dendrogram plot ===\n")
n_genes          <- nrow(vst)
colors_by_order  <- module_colors[gene_tree$order]

png("plots/wgcna_dendrogram_modules.png", width = 1400, height = 700, res = 150)
layout(matrix(c(1, 2), nrow = 2), heights = c(3, 1))

par(mar = c(0, 4, 3, 1))
plot(gene_tree, labels = FALSE, hang = 0.04,
     main = sprintf("WGCNA Gene Dendrogram  (1 − TOM, β = %d, k = %d)",
                    soft_power, k_best),
     xlab = "", sub = "", ylab = "1 − TOM distance")

par(mar = c(2, 4, 0.5, 1))
plot.new()
plot.window(xlim = c(0, n_genes), ylim = c(0, 1))
for (i in seq_along(colors_by_order)) {
  rect(i - 1, 0, i, 1, col = colors_by_order[i], border = NA)
}
mtext("Module", side = 2, las = 1, line = 0.3, cex = 0.8)

dev.off()
cat("Saved: plots/wgcna_dendrogram_modules.png\n")

# ============================================================
# 10. Compare with original correlation-based modules (M1, M2)
# ============================================================
cat("\n=== 10. Comparison with original hclust/cutree modules ===\n")
old_mods <- read.csv("results/module_assignments.csv", stringsAsFactors = FALSE)
# column names: gene_symbol, module
colnames(old_mods) <- c("gene", "old_module")

comp_df <- merge(gene_modules[, c("gene", "wgcna_module", "kWithin")],
                 old_mods, by = "gene", all.x = TRUE)

cat("\nCross-tabulation: WGCNA module vs original module\n")
print(table(WGCNA = comp_df$wgcna_module, Original = comp_df$old_module))

# Concordance within M1 and M2
for (m in c("M1", "M2")) {
  genes_in_old  <- old_mods$gene[old_mods$old_module == m]
  wgcna_of_old  <- gene_modules$wgcna_module[gene_modules$gene %in% genes_in_old]
  cat(sprintf("\nOriginal %s genes (%d) distribute into WGCNA modules:\n",
              m, length(wgcna_of_old)))
  print(sort(table(wgcna_of_old), decreasing = TRUE))
}

write.csv(gene_modules, "results/wgcna_gene_modules.csv",        row.names = FALSE)
write.csv(comp_df,      "results/wgcna_vs_old_comparison.csv",   row.names = FALSE)
write.csv(key_df,       "results/wgcna_key_gene_assignments.csv", row.names = FALSE)

cat("\n============================================================\n")
cat("WGCNA coexpression analysis complete.\n")
cat(sprintf("  β = %d  |  k = %d modules  |  139 genes × 143 samples\n",
            soft_power, k_best))
cat("============================================================\n")
