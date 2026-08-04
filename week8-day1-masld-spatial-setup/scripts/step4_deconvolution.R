#!/usr/bin/env Rscript
# Step 4: Cell-type deconvolution for MASLD Visium dataset
# Loads saved combined Seurat object; applies AddModuleScore with GSE136103 reference

suppressMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

PROJ_ROOT <- "/Users/senaesen/Desktop/bioinfo-learning/week8-day1-masld-spatial-setup"
REF_RDS   <- "/Users/senaesen/Desktop/bioinfo-learning/week6-day1-nafld-scrna-gse136103/results/nafld_seurat_annotated.rds"
RESULTS   <- file.path(PROJ_ROOT, "results")
DECON_DIR <- file.path(RESULTS, "deconvolution")
dir.create(DECON_DIR, recursive=TRUE, showWarnings=FALSE)

cat("=== STEP 4: Cell-type Deconvolution (AddModuleScore) ===\n")
cat("Loading saved combined spatial object...\n")
combined <- readRDS(file.path(RESULTS, "masld_spatial_combined.rds"))
cat(sprintf("  %d spots × %d genes, %d clusters\n",
    ncol(combined), nrow(combined), length(unique(combined$seurat_clusters))))

cat("Loading GSE136103 scRNA-seq reference...\n")
ref_obj <- readRDS(REF_RDS)
cat(sprintf("  Reference: %d genes × %d cells, cell types: %s\n",
    nrow(ref_obj), ncol(ref_obj),
    paste(sort(unique(ref_obj$cell_type)), collapse=", ")))

# Build marker gene sets per cell type (top 50 avg-expressing genes per type)
cat("Building cell-type marker gene lists...\n")
DefaultAssay(ref_obj) <- "RNA"
ref_obj <- NormalizeData(ref_obj, verbose=FALSE)
ref_norm <- GetAssayData(ref_obj, assay="RNA", layer="data")

cell_types <- sort(unique(ref_obj$cell_type))
marker_genes <- setNames(vector("list", length(cell_types)), cell_types)

for (ct in cell_types) {
  ct_cells <- which(ref_obj$cell_type == ct)
  if (length(ct_cells) < 5) next
  avg_exp <- rowMeans(ref_norm[, ct_cells, drop=FALSE])
  top50   <- names(sort(avg_exp, decreasing=TRUE))[1:50]
  in_spatial <- intersect(top50, rownames(combined))
  if (length(in_spatial) >= 5) {
    marker_genes[[ct]] <- in_spatial
    cat(sprintf("  %-35s : %d/%d marker genes found in spatial data\n",
                ct, length(in_spatial), 50))
  }
}
marker_genes <- Filter(Negate(is.null), marker_genes)
cat(sprintf("  %d cell types with sufficient markers\n\n", length(marker_genes)))

# Score each cell type in combined spatial object
cat("Scoring cell types with AddModuleScore...\n")
DefaultAssay(combined) <- "SCT"

ct_module_names <- character(length(marker_genes))
for (i in seq_along(marker_genes)) {
  ct <- names(marker_genes)[i]
  safe_name <- paste0("CT", i, "_", gsub("[^A-Za-z0-9]", "_", substr(ct, 1, 20)))
  combined <- AddModuleScore(
    combined,
    features = list(marker_genes[[ct]]),
    name     = safe_name,
    ctrl     = 100,
    verbose  = FALSE
  )
  # AddModuleScore appends "1" to the name
  ct_module_names[i] <- paste0(safe_name, "1")
  cat(sprintf("  Scored: %-35s → column %s\n", ct, ct_module_names[i]))
}

# Rename for readability
names(ct_module_names) <- names(marker_genes)

# Extract spot-level scores
score_df <- combined@meta.data[, c("array_id", "seurat_clusters", ct_module_names)]
colnames(score_df)[colnames(score_df) %in% ct_module_names] <- names(ct_module_names)
score_df$barcode <- rownames(score_df)

write.csv(score_df, file.path(DECON_DIR, "spot_celltype_scores.csv"), row.names=FALSE)
cat("\nSpot-level scores saved.\n")

# Mean score per cluster (for cluster biological annotation)
score_cols <- names(marker_genes)
mean_by_cluster <- score_df %>%
  group_by(seurat_clusters) %>%
  summarise(across(all_of(score_cols), mean), .groups="drop")
write.csv(mean_by_cluster, file.path(DECON_DIR, "mean_celltype_scores_per_cluster.csv"), row.names=FALSE)

# Mean score per array
mean_by_array <- score_df %>%
  group_by(array_id) %>%
  summarise(across(all_of(score_cols), mean), .groups="drop")
write.csv(mean_by_array, file.path(DECON_DIR, "mean_celltype_scores_per_array.csv"), row.names=FALSE)

# ── Heatmap: cell-type scores by cluster ──────────────────────────────────────
mat_clust <- as.matrix(mean_by_cluster[, score_cols])
rownames(mat_clust) <- paste0("C", mean_by_cluster$seurat_clusters)

# Scale scores across clusters per cell type (z-score) for better color contrast
mat_z <- apply(mat_clust, 2, function(x) (x - mean(x)) / (sd(x) + 1e-8))

df_long <- as.data.frame(mat_z) %>%
  mutate(Cluster = rownames(mat_z)) %>%
  pivot_longer(-Cluster, names_to="CellType", values_to="ZScore")

# Order clusters and cell types
df_long$Cluster  <- factor(df_long$Cluster, levels=paste0("C", 0:23))
df_long$CellType <- factor(df_long$CellType,
  levels = score_cols[order(apply(mat_z, 2, which.max))])

p_heat_clust <- ggplot(df_long, aes(CellType, Cluster, fill=ZScore)) +
  geom_tile(color="white", linewidth=0.2) +
  scale_fill_gradient2(
    low="steelblue", mid="white", high="firebrick",
    midpoint=0, name="Z-score\n(per cell type)"
  ) +
  theme_bw(base_size=9) +
  theme(
    axis.text.x = element_text(angle=50, hjust=1, size=7),
    axis.text.y = element_text(size=7),
    panel.grid  = element_blank()
  ) +
  labs(
    title = "Cell-Type Activity Scores per Spatial Cluster",
    subtitle = "Z-scored AddModuleScore means; reference: GSE136103 human liver (Ramachandran 2019)",
    x = NULL, y = "Spatial Cluster"
  )

ggsave(file.path(DECON_DIR, "celltype_scores_heatmap_by_cluster.png"),
       p_heat_clust, width=14, height=8, dpi=150)

# ── Heatmap: cell-type scores by array ───────────────────────────────────────
mat_arr <- as.matrix(mean_by_array[, score_cols])
rownames(mat_arr) <- mean_by_array$array_id
mat_arr_z <- apply(mat_arr, 2, function(x) (x - mean(x)) / (sd(x) + 1e-8))

df_arr_long <- as.data.frame(mat_arr_z) %>%
  mutate(Array = rownames(mat_arr_z)) %>%
  pivot_longer(-Array, names_to="CellType", values_to="ZScore")
df_arr_long$CellType <- factor(df_arr_long$CellType, levels=levels(df_long$CellType))

p_heat_arr <- ggplot(df_arr_long, aes(CellType, Array, fill=ZScore)) +
  geom_tile(color="white", linewidth=0.3) +
  scale_fill_gradient2(low="steelblue", mid="white", high="firebrick",
                       midpoint=0, name="Z-score") +
  theme_bw(base_size=9) +
  theme(axis.text.x=element_text(angle=50, hjust=1, size=7)) +
  labs(title="Mean Cell-Type Activity per Visium Array",
       x=NULL, y=NULL)

ggsave(file.path(DECON_DIR, "celltype_scores_heatmap_by_array.png"),
       p_heat_arr, width=14, height=5, dpi=150)

# ── Identify dominant cell type per cluster ───────────────────────────────────
dominant <- apply(mat_clust, 1, function(row) names(which.max(row)))
cluster_annotation <- data.frame(
  cluster         = rownames(mat_clust),
  n_spots         = table(combined$seurat_clusters),
  dominant_celltype = dominant
)
write.csv(cluster_annotation, file.path(DECON_DIR, "cluster_dominant_celltype.csv"), row.names=FALSE)

cat("\nDominant cell type per cluster:\n")
print(cluster_annotation[, c("cluster","dominant_celltype")])

# ── Save updated combined object with scores ──────────────────────────────────
saveRDS(combined, file.path(RESULTS, "masld_spatial_combined.rds"))

cat("\n=== Step 4 complete ===\n")
cat("Outputs:\n")
cat("  deconvolution/spot_celltype_scores.csv\n")
cat("  deconvolution/mean_celltype_scores_per_cluster.csv\n")
cat("  deconvolution/mean_celltype_scores_per_array.csv\n")
cat("  deconvolution/cluster_dominant_celltype.csv\n")
cat("  deconvolution/celltype_scores_heatmap_by_cluster.png\n")
cat("  deconvolution/celltype_scores_heatmap_by_array.png\n")
cat("  results/masld_spatial_combined.rds  (updated with scores)\n")
