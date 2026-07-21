################################################################################
# Week 3 Day 3 — Spatial Transcriptomics: Human Heart (10x Visium)
# Dataset: V1_Human_Heart — 10x Genomics Visium spatial gene expression
# Run from: week3-day3-spatial-transcriptomics-heart/
#   Rscript scripts/spatial_heart.R
################################################################################

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

options(future.globals.maxSize = 2000 * 1024^2)  # 2 GB for SCTransform

# ── Paths ─────────────────────────────────────────────────────────────────────
data_dir    <- "data/"
out_plots   <- "plots/"
out_results <- "results/"

################################################################################
# Step 1 — Load data
################################################################################
cat("\n=== Step 1: Load 10x Visium data ===\n")

# Seurat 5 Load10X_Spatial defaults to .h5; use Read10X + Read10X_Image for MTX format
counts <- Read10X(data.dir = file.path(data_dir, "filtered_feature_bc_matrix"))
heart  <- CreateSeuratObject(counts, assay = "Spatial")

img <- Read10X_Image(
  image.dir  = file.path(data_dir, "spatial"),
  image.name = "tissue_lowres_image.png",
  filter.matrix = TRUE
)
img <- img[Cells(heart)]
DefaultAssay(object = img) <- "Spatial"
heart[["heart"]] <- img

cat("Spots:", ncol(heart), "\n")
cat("Genes:", nrow(heart), "\n")

################################################################################
# Step 2 — QC
################################################################################
cat("\n=== Step 2: Quality control ===\n")

heart[["percent.mt"]] <- PercentageFeatureSet(heart, pattern = "^MT-")

# Plot QC metrics on tissue
p_qc1 <- SpatialFeaturePlot(heart, features = "nCount_Spatial") +
  labs(title = "UMI counts per spot") +
  theme(legend.position = "right")

p_qc2 <- SpatialFeaturePlot(heart, features = "nFeature_Spatial") +
  labs(title = "Genes detected per spot") +
  theme(legend.position = "right")

p_qc3 <- SpatialFeaturePlot(heart, features = "percent.mt") +
  labs(title = "Mitochondrial % per spot") +
  theme(legend.position = "right")

ggsave(file.path(out_plots, "01_qc_spatial.png"),
       p_qc1 | p_qc2 | p_qc3,
       width = 15, height = 5, dpi = 150)
cat("Plot 01 saved: QC spatial maps\n")

# Violin plots of QC metrics
p_vln <- VlnPlot(heart,
                 features = c("nCount_Spatial", "nFeature_Spatial", "percent.mt"),
                 pt.size = 0.1, ncol = 3) &
  theme(axis.text.x = element_blank())
ggsave(file.path(out_plots, "02_qc_violin.png"), p_vln,
       width = 12, height = 4, dpi = 150)
cat("Plot 02 saved: QC violin\n")

n_spots_before_qc <- ncol(heart)   # FIX: capture BEFORE subsetting (see qc_summary below)
cat("Before QC — spots:", n_spots_before_qc, "\n")
cat("nCount_Spatial range:", min(heart$nCount_Spatial), "-", max(heart$nCount_Spatial), "\n")
cat("percent.mt range:", round(min(heart$percent.mt), 2), "-", round(max(heart$percent.mt), 2), "\n")

# Heart tissue has naturally high mitochondrial content (cardiac cells are
# ~30-50% MT by design). Use a permissive MT threshold; filter only empty
# droplets / off-tissue spots by UMI count.
heart <- subset(heart,
                nCount_Spatial > 200 &
                nFeature_Spatial > 100 &
                percent.mt < 60)
cat("After QC — spots:", ncol(heart), "\n")

################################################################################
# Step 3 — Normalisation with SCTransform
################################################################################
cat("\n=== Step 3: SCTransform normalisation ===\n")
cat("SCTransform accounts for technical variation in UMI capture across spots.\n")

heart <- SCTransform(heart, assay = "Spatial", verbose = FALSE)
cat("SCTransform complete.\n")

################################################################################
# Step 4 — PCA and UMAP
################################################################################
cat("\n=== Step 4: Dimensionality reduction (PCA + UMAP) ===\n")

heart <- RunPCA(heart, assay = "SCT", verbose = FALSE)

# Elbow plot
p_elbow <- ElbowPlot(heart, ndims = 30) +
  labs(title = "PCA elbow plot — Visium heart") +
  theme_bw()
ggsave(file.path(out_plots, "03_pca_elbow.png"), p_elbow,
       width = 6, height = 4, dpi = 150)
cat("Plot 03 saved: PCA elbow\n")

heart <- RunUMAP(heart, dims = 1:30, verbose = FALSE)
cat("PCA + UMAP complete.\n")

################################################################################
# Step 5 — Clustering
################################################################################
cat("\n=== Step 5: Clustering ===\n")

heart <- FindNeighbors(heart, dims = 1:30, verbose = FALSE)
heart <- FindClusters(heart, resolution = 0.5, verbose = FALSE)

n_clusters <- length(unique(heart$seurat_clusters))
cat("Clusters found:", n_clusters, "\n")
cat("Cluster sizes:\n")
print(table(heart$seurat_clusters))

# UMAP coloured by cluster
p_umap <- DimPlot(heart, reduction = "umap", label = TRUE, pt.size = 1.5) +
  labs(title = "UMAP — Visium heart clusters") +
  theme_bw()
ggsave(file.path(out_plots, "04_umap_clusters.png"), p_umap,
       width = 7, height = 6, dpi = 150)
cat("Plot 04 saved: UMAP clusters\n")

# Spatial plot coloured by cluster
p_spatial_cluster <- SpatialDimPlot(heart, label = TRUE, label.size = 3) +
  labs(title = "Spatial clusters — Human Heart")
ggsave(file.path(out_plots, "05_spatial_clusters.png"), p_spatial_cluster,
       width = 7, height = 7, dpi = 150)
cat("Plot 05 saved: spatial cluster map\n")

################################################################################
# Step 6 — Spatially variable genes
################################################################################
cat("\n=== Step 6: Spatially variable genes ===\n")
cat("Finding genes whose expression varies across the tissue section...\n")

heart <- FindSpatiallyVariableFeatures(
  heart,
  assay    = "SCT",
  features = VariableFeatures(heart)[1:1000],
  selection.method = "moransi"
)

top_svf <- head(SpatiallyVariableFeatures(heart, method = "moransi"), 12)
cat("Top spatially variable genes:\n")
print(top_svf)

write.csv(
  data.frame(gene = SpatiallyVariableFeatures(heart, method = "moransi")),
  file.path(out_results, "spatially_variable_genes.csv"),
  row.names = FALSE
)
cat("Spatially variable genes saved to results/\n")

# Spatial feature plots for top 6 SVGs
p_svg <- SpatialFeaturePlot(
  heart,
  features = top_svf[1:6],
  ncol = 3,
  alpha = c(0.1, 1)
)
ggsave(file.path(out_plots, "06_spatially_variable_genes.png"), p_svg,
       width = 14, height = 9, dpi = 150)
cat("Plot 06 saved: top spatially variable genes\n")

################################################################################
# Step 7 — Cluster marker genes
################################################################################
cat("\n=== Step 7: Cluster marker genes ===\n")

DefaultAssay(heart) <- "SCT"
markers <- FindAllMarkers(
  heart,
  only.pos         = TRUE,
  min.pct          = 0.25,
  logfc.threshold  = 0.25,
  verbose          = FALSE
)

top5 <- markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 5) %>%
  ungroup()

write.csv(markers, file.path(out_results, "cluster_markers_all.csv"),  row.names = FALSE)
write.csv(top5,    file.path(out_results, "cluster_markers_top5.csv"), row.names = FALSE)
cat("Marker genes saved. Total markers:", nrow(markers), "\n")

# Dot plot of top 2 markers per cluster
top2 <- markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 2) %>%
  ungroup() %>%
  pull(gene) %>%
  unique()

p_dot <- DotPlot(heart, features = top2, assay = "SCT") +
  RotatedAxis() +
  labs(title = "Top 2 marker genes per cluster — Human Heart") +
  theme_bw(base_size = 10)
ggsave(file.path(out_plots, "07_marker_dotplot.png"), p_dot,
       width = 12, height = 5, dpi = 150)
cat("Plot 07 saved: marker dot plot\n")

# Spatial feature plots for top marker of each cluster
top1_per_cluster <- markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 1) %>%
  ungroup() %>%
  pull(gene) %>%
  unique()

n_plot <- min(9, length(top1_per_cluster))
p_markers_spatial <- SpatialFeaturePlot(
  heart,
  features = top1_per_cluster[1:n_plot],
  ncol     = 3,
  alpha    = c(0.1, 1)
)
ggsave(file.path(out_plots, "08_top_marker_spatial.png"), p_markers_spatial,
       width = 14, height = ceiling(n_plot / 3) * 5, dpi = 150)
cat("Plot 08 saved: spatial maps of top cluster markers\n")

################################################################################
# Step 8 — Known cardiac marker genes
################################################################################
cat("\n=== Step 8: Known cardiac marker genes ===\n")

cardiac_markers <- c(
  "TNNT2",   # cardiomyocytes
  "MYH7",    # ventricular myosin heavy chain
  "MYH6",    # atrial myosin heavy chain
  "VIM",     # fibroblasts / mesenchyme
  "DCN",     # fibroblasts
  "PECAM1",  # endothelial cells (CD31)
  "ACTA2",   # smooth muscle
  "CD68"     # macrophages
)

# Only plot markers present in the dataset
cardiac_present <- cardiac_markers[cardiac_markers %in% rownames(heart)]
cat("Cardiac markers found in dataset:", paste(cardiac_present, collapse=", "), "\n")

if (length(cardiac_present) >= 2) {
  p_cardiac <- SpatialFeaturePlot(
    heart,
    features = cardiac_present,
    ncol     = 4,
    alpha    = c(0.1, 1)
  )
  ggsave(file.path(out_plots, "09_cardiac_markers_spatial.png"), p_cardiac,
         width = 18, height = ceiling(length(cardiac_present) / 4) * 5, dpi = 150)
  cat("Plot 09 saved: cardiac marker spatial maps\n")
}

# Violin plots of cardiac markers across clusters
if (length(cardiac_present) >= 2) {
  p_vln_cardiac <- VlnPlot(heart, features = cardiac_present[1:min(4, length(cardiac_present))],
                            ncol = 4, pt.size = 0) &
    theme(axis.text.x = element_text(angle=45, hjust=1, size=8))
  ggsave(file.path(out_plots, "10_cardiac_markers_violin.png"), p_vln_cardiac,
         width = 14, height = 4, dpi = 150)
  cat("Plot 10 saved: cardiac marker violin plots\n")
}

################################################################################
# Step 9 — Save results
################################################################################
cat("\n=== Step 9: Saving results ===\n")

# Cluster composition table
cluster_table <- as.data.frame(table(Cluster = heart$seurat_clusters))
write.csv(cluster_table, file.path(out_results, "cluster_sizes.csv"), row.names = FALSE)

# QC summary
qc_summary <- data.frame(
  metric  = c("spots_before_qc", "spots_after_qc", "n_genes", "n_clusters"),
  # FIX: use the pre-subset count captured before filtering. The old expression
  # `ncol(heart) + sum(heart$nCount_Spatial <= 500)` ran AFTER subsetting to
  # nCount_Spatial > 200, so it double-counted surviving spots and never counted
  # the ones actually removed — the "before" number was wrong.
  value   = c(n_spots_before_qc,
              ncol(heart), nrow(heart), n_clusters)
)
write.csv(qc_summary, file.path(out_results, "qc_summary.csv"), row.names = FALSE)

# Save Seurat object (excluded from git via .gitignore)
saveRDS(heart, file.path(out_results, "heart_seurat.rds"))
cat("Seurat object saved: results/heart_seurat.rds\n")

cat("\n==============================\n")
cat("Pipeline complete.\n")
cat("Plots:   10 files in plots/\n")
cat("Results: cluster_sizes, cluster_markers, spatially_variable_genes in results/\n")
cat("==============================\n")
