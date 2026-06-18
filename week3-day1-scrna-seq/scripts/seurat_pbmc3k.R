# =============================================================================
# Week 3 — Day 1 — Seurat PBMC 3k scRNA-seq Pipeline
# =============================================================================
# Classic 10x Genomics tutorial dataset: 2,700 peripheral blood mononuclear
# cells (PBMCs) from a healthy donor, sequenced on the 10x Chromium platform.
# This is the standard learning dataset for scRNA-seq — every paper and package
# tutorial uses it.
#
# Pipeline:
#   1. Load data
#   2. Quality control (QC)
#   3. Normalisation
#   4. Highly variable genes
#   5. Scaling
#   6. PCA
#   7. Clustering (KNN graph + Leiden/Louvain)
#   8. UMAP
#   9. Marker genes
#  10. Cell type annotation
#
# Run from: week3-day1-scrna-seq/
#   Rscript scripts/seurat_pbmc3k.R
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

if (!interactive()) pdf(NULL)

cat("=== Seurat PBMC 3k Pipeline ===\n\n")


# =============================================================================
# STEP 1 — Load the count matrix
# =============================================================================
# Read10X() reads the three files that make up a 10x count matrix:
#   barcodes.tsv — one cell barcode per line
#   genes.tsv    — one gene per line (Ensembl ID + gene symbol)
#   matrix.mtx   — sparse matrix of UMI counts (genes × cells)
#
# CreateSeuratObject() wraps this into a Seurat object — the core data
# structure used throughout the analysis. It stores the count matrix,
# metadata per cell, and all downstream results (PCA, UMAP, clusters).
#
# min.cells = 3:  keep only genes detected in at least 3 cells
# min.features = 200: keep only cells with at least 200 detected genes
# =============================================================================

cat("STEP 1: Loading count matrix...\n")

pbmc.data <- Read10X(data.dir = "filtered_gene_bc_matrices/hg19/")
pbmc <- CreateSeuratObject(
  counts      = pbmc.data,
  project     = "PBMC3k",
  min.cells   = 3,
  min.features = 200
)

cat("  Cells loaded:", ncol(pbmc), "\n")
cat("  Genes loaded:", nrow(pbmc), "\n\n")


# =============================================================================
# STEP 2 — Quality Control (QC)
# =============================================================================
# Three metrics identify low-quality cells:
#
#   nFeature_RNA  — number of distinct genes detected per cell
#                   too low  → empty droplet or dead cell
#                   too high → probable doublet (two cells in one droplet)
#
#   nCount_RNA    — total UMI count per cell
#                   mirrors nFeature_RNA; low = poor quality, high = doublet
#
#   percent.mt    — % of UMIs from mitochondrial genes
#                   when a cell's outer membrane is damaged, cytoplasmic RNA
#                   leaks out. Mitochondria (membrane-enclosed) stay intact.
#                   Result: dying/dead cells have disproportionately high MT %.
#                   Human mitochondrial genes all start with "MT-".
#
# Thresholds used:
#   nFeature_RNA > 200  (remove empties)
#   nFeature_RNA < 2500 (remove doublets)
#   percent.mt   < 5    (remove dead/dying cells)
# =============================================================================

cat("STEP 2: Quality control...\n")

# Calculate mitochondrial percentage
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Plot QC metrics — violin plots per cell
p_qc <- VlnPlot(
  pbmc,
  features  = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol      = 3,
  pt.size   = 0.1
) & theme(plot.title = element_text(size = 10))

ggsave("plots/01_qc_violin.png", p_qc, width = 12, height = 5, dpi = 150)

# Scatter plots: nCount vs percent.mt, nCount vs nFeature
p_scatter1 <- FeatureScatter(pbmc, "nCount_RNA", "percent.mt")
p_scatter2 <- FeatureScatter(pbmc, "nCount_RNA", "nFeature_RNA")
ggsave("plots/02_qc_scatter.png",
       p_scatter1 + p_scatter2, width = 10, height = 5, dpi = 150)

cat("  Before QC: ", ncol(pbmc), "cells\n")

# Apply QC filters
pbmc <- subset(pbmc,
               subset = nFeature_RNA > 200 &
                        nFeature_RNA < 2500 &
                        percent.mt < 5)

cat("  After QC:  ", ncol(pbmc), "cells\n\n")


# =============================================================================
# STEP 3 — Normalisation
# =============================================================================
# Different cells capture different amounts of RNA due to technical variation,
# not biology. We normalise to remove this before comparing cells.
#
# NormalizeData() applies:
#   1. Divide each cell's counts by its total UMIs × 10,000
#      (puts all cells on the same "per 10k molecules" scale — like CPM)
#   2. log1p: log(value + 1)
#      Compresses the dynamic range; makes expression more normally distributed.
#      The +1 prevents log(0) = -infinity for undetected genes.
#
# The normalised values are stored in pbmc[["RNA"]]$data
# =============================================================================

cat("STEP 3: Normalising...\n")

pbmc <- NormalizeData(pbmc,
                      normalization.method = "LogNormalize",
                      scale.factor         = 10000,
                      verbose              = FALSE)

cat("  Done. Normalised counts stored in RNA$data\n\n")


# =============================================================================
# STEP 4 — Highly Variable Genes (HVGs)
# =============================================================================
# We have ~13,000 genes but not all are informative. Most are "housekeeping"
# genes that are expressed at similar levels in every cell — they add noise,
# not signal.
#
# FindVariableFeatures() identifies genes that vary MORE than expected for
# their average expression level (using a variance-stabilising method).
# These are the genes that distinguish different cell types and states.
#
# nfeatures = 2000 is the Seurat default — a good starting point for PBMCs.
# Only these 2000 genes are used for PCA and clustering.
# =============================================================================

cat("STEP 4: Finding highly variable genes...\n")

pbmc <- FindVariableFeatures(pbmc,
                             selection.method = "vst",
                             nfeatures        = 2000,
                             verbose          = FALSE)

# Plot: top 10 most variable genes labelled
top10 <- head(VariableFeatures(pbmc), 10)
p_hvg <- VariableFeaturePlot(pbmc)
p_hvg <- LabelPoints(plot = p_hvg, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)
ggsave("plots/03_variable_features.png", p_hvg, width = 9, height = 5, dpi = 150)

cat("  Top 10 HVGs:", paste(top10, collapse = ", "), "\n\n")


# =============================================================================
# STEP 5 — Scaling
# =============================================================================
# Before PCA, all genes need to be on the same scale. Without scaling, genes
# with high absolute expression would dominate the PCA simply because of their
# magnitude — not because they're more biologically informative.
#
# ScaleData() z-scores each gene across cells:
#   scaled_value = (value - mean) / standard_deviation
#
# After scaling, every gene has mean = 0 and variance = 1.
# This is applied only to the 2000 HVGs (default) to save time.
# Scaled values are stored in pbmc[["RNA"]]$scale.data
# =============================================================================

cat("STEP 5: Scaling data...\n")
pbmc <- ScaleData(pbmc, verbose = FALSE)
cat("  Done. Scaled matrix stored in RNA$scale.data\n\n")


# =============================================================================
# STEP 6 — PCA (Principal Component Analysis)
# =============================================================================
# Even 2000 genes is too many dimensions to cluster in directly. PCA compresses
# these into a smaller set of components (PCs) that capture most of the variance.
#
# Each PC is a linear combination of genes that explains a dimension of
# variation across cells:
#   PC1 might separate T cells from B cells
#   PC2 might separate monocytes from lymphocytes
#   PC3 might reflect cell cycle state
#
# npcs = 50: compute the top 50 PCs. We'll use only the top 20 or so for
# clustering (determined by the elbow plot). Using too many PCs adds noise
# from components that capture technical variation rather than biology.
# =============================================================================

cat("STEP 6: Running PCA...\n")

pbmc <- RunPCA(pbmc, features = VariableFeatures(pbmc),
               npcs = 50, verbose = FALSE)

# Elbow plot: variance explained by each PC
# The "elbow" — where the curve flattens — shows how many PCs capture real signal
p_elbow <- ElbowPlot(pbmc, ndims = 50) +
  geom_vline(xintercept = 10, linetype = "dashed", colour = "red") +
  labs(subtitle = "Red dashed line = 10 PCs used for clustering")
ggsave("plots/04_elbow_plot.png", p_elbow, width = 7, height = 5, dpi = 150)

# Heatmap of top genes driving PC1 and PC2
png("plots/05_pca_heatmap.png", width = 1200, height = 800, res = 120)
DimHeatmap(pbmc, dims = 1:6, cells = 500, balanced = TRUE, fast = FALSE)
dev.off()

cat("  Top genes in PC1:", paste(head(pbmc@reductions$pca@feature.loadings[
  order(abs(pbmc@reductions$pca@feature.loadings[,1]), decreasing=TRUE), 1], 5)
  |> names(), collapse=", "), "\n\n")


# =============================================================================
# STEP 7 — Clustering
# =============================================================================
# Two sub-steps:
#
# FindNeighbors(): builds a K-nearest-neighbour (KNN) graph in PCA space.
#   For each cell, find its k=20 most similar cells (by Euclidean distance
#   in the top 10 PCs). Draw edges between them. Result: a graph where
#   connected cells are transcriptionally similar.
#   dims = 1:10 — use top 10 PCs (as suggested by elbow plot)
#
# FindClusters(): applies the Louvain community detection algorithm to the
#   KNN graph. Finds groups of cells that are more connected to each other
#   than to the rest of the graph.
#   resolution = 0.5: controls granularity.
#     Higher → more, smaller clusters. Lower → fewer, larger clusters.
#     0.5 is a good starting point for ~2,700 PBMCs.
# =============================================================================

cat("STEP 7: Clustering...\n")

pbmc <- FindNeighbors(pbmc, dims = 1:10, verbose = FALSE)
pbmc <- FindClusters(pbmc, resolution = 0.5, verbose = FALSE)

n_clusters <- length(levels(pbmc$seurat_clusters))
cat("  Clusters found:", n_clusters, "\n")
cat("  Cells per cluster:\n")
print(table(pbmc$seurat_clusters))
cat("\n")


# =============================================================================
# STEP 8 — UMAP
# =============================================================================
# UMAP (Uniform Manifold Approximation and Projection) projects the high-
# dimensional PCA embedding into 2D for visualisation.
#
# It tries to preserve the local structure: cells that are neighbours in
# 10-dimensional PCA space end up close together on the 2D plot.
#
# IMPORTANT: UMAP is ONLY for visualisation. The 2D coordinates are not
# used for clustering or statistics. Distances between clusters on UMAP
# are not meaningful — two clusters far apart on UMAP may not actually be
# more different than two clusters sitting close together.
#
# We use dims = 1:10 (same PCs as clustering) for consistency.
# =============================================================================

cat("STEP 8: Running UMAP...\n")

pbmc <- RunUMAP(pbmc, dims = 1:10, verbose = FALSE, seed.use = 42)

# UMAP coloured by cluster
p_umap_clusters <- DimPlot(pbmc, reduction = "umap", label = TRUE,
                            label.size = 5, repel = TRUE) +
  labs(title = "PBMC 3k — Seurat clusters (res = 0.5)",
       subtitle = "UMAP coloured by cluster identity") +
  theme_bw(base_size = 11)

ggsave("plots/06_umap_clusters.png", p_umap_clusters,
       width = 8, height = 7, dpi = 150)

cat("  UMAP done. Saved to plots/06_umap_clusters.png\n\n")


# =============================================================================
# STEP 9 — Find Marker Genes per Cluster
# =============================================================================
# FindAllMarkers() compares each cluster against all other clusters to find
# genes that are specifically expressed in that cluster (Wilcoxon rank-sum test).
#
# Parameters:
#   only.pos = TRUE  — only report upregulated markers (easier to interpret)
#   min.pct  = 0.25  — gene must be expressed in ≥25% of cells in the cluster
#   logfc.threshold = 0.25 — minimum log fold change (speeds up computation)
#
# Output: table with avg_log2FC (fold change), p_val_adj (FDR), and pct.1/pct.2
# (fraction of cells expressing the gene in the cluster vs. all other clusters)
# =============================================================================

cat("STEP 9: Finding cluster marker genes...\n")

markers <- FindAllMarkers(
  pbmc,
  only.pos        = TRUE,
  min.pct         = 0.25,
  logfc.threshold = 0.25,
  verbose         = FALSE
)

# Top 5 markers per cluster
top5 <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5) %>%
  ungroup()

cat("  Top marker per cluster:\n")
top1 <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 1)
print(data.frame(cluster = top1$cluster, gene = top1$gene,
                 log2FC = round(top1$avg_log2FC, 2)), row.names = FALSE)
cat("\n")

# Heatmap of top 5 markers per cluster
pbmc <- ScaleData(pbmc, features = unique(top5$gene), verbose = FALSE)
png("plots/07_marker_heatmap.png", width = 1400, height = 900, res = 120)
DoHeatmap(pbmc, features = unique(top5$gene), size = 3) +
  theme(axis.text.y = element_text(size = 6))
dev.off()
cat("  Saved: plots/07_marker_heatmap.png\n")

# Dot plot of canonical PBMC markers across clusters
canonical_markers <- c(
  "MS4A1", "CD79A",           # B cells
  "CD3D", "CD3E",             # T cells (all)
  "CD8A", "CD8B",             # CD8+ T cells
  "IL7R", "CCR7",             # CD4+ naive T
  "S100A4",                   # Memory T
  "GNLY", "NKG7",             # NK cells
  "CD14", "LYZ",              # CD14+ monocytes
  "FCGR3A", "MS4A7",          # CD16+ monocytes
  "FCER1A", "CST3",           # Dendritic cells
  "PPBP"                      # Platelets
)

p_dot <- DotPlot(pbmc, features = canonical_markers) +
  RotatedAxis() +
  labs(title = "Canonical PBMC marker expression across clusters") +
  theme(axis.text.x = element_text(size = 8))

ggsave("plots/08_marker_dotplot.png", p_dot, width = 12, height = 6, dpi = 150)
cat("  Saved: plots/08_marker_dotplot.png\n\n")


# =============================================================================
# STEP 10 — Cell Type Annotation
# =============================================================================
# Based on canonical marker genes, we assign a cell type identity to each
# cluster. This combines:
#   - The marker gene output from Step 9
#   - Known biology of peripheral blood cell types
#   - The dot plot showing which markers are expressed in which clusters
#
# Mapping based on well-established PBMC markers:
#   MS4A1 / CD79A → B cells
#   CD3D + IL7R + CCR7 → Naive CD4+ T cells
#   CD3D + IL7R + S100A4 → Memory CD4+ T cells
#   CD3D + CD8A → CD8+ T cells
#   GNLY + NKG7 → NK cells
#   CD14 + LYZ → CD14+ Monocytes
#   FCGR3A + MS4A7 → CD16+ Monocytes
#   FCER1A + CST3 → Dendritic cells
#   PPBP → Platelets
#
# These assignments are specific to this dataset at resolution = 0.5.
# A different resolution or dataset requires re-evaluation.
# =============================================================================

cat("STEP 10: Annotating cell types...\n")

new_labels <- c(
  "0" = "CD4+ Naive T",      # IL7R+, CCR7+
  "1" = "CD14+ Monocytes",   # CD14+, LYZ+, CST3+
  "2" = "CD4+ Memory T",     # IL7R+, S100A4+, lower CCR7
  "3" = "B cells",           # MS4A1+, CD79A+
  "4" = "CD8+ T",            # CD8A+, NKG7+
  "5" = "CD16+ Monocytes",   # FCGR3A+, MS4A7+
  "6" = "NK cells",          # GNLY++, NKG7++
  "7" = "Dendritic cells",   # FCER1A+, CST3+, LYZ+
  "8" = "Platelets"          # PPBP++
)

pbmc <- RenameIdents(pbmc, new_labels)

# Add annotation to metadata
pbmc$cell_type <- as.character(Idents(pbmc))

# UMAP coloured by cell type
cell_type_colours <- c(
  "CD4+ Naive T"     = "#4DBBD5",
  "CD14+ Monocytes"  = "#E64B35",
  "CD4+ Memory T"    = "#00A087",
  "B cells"          = "#F39B7F",
  "CD8+ T"           = "#3C5488",
  "CD16+ Monocytes"  = "#91D1C2",
  "NK cells"         = "#8491B4",
  "Dendritic cells"  = "#DC0000",
  "Platelets"        = "#7E6148"
)

p_umap_annotated <- DimPlot(
  pbmc, reduction = "umap", label = TRUE,
  label.size = 4, repel = TRUE,
  cols = cell_type_colours
) +
  labs(title   = "PBMC 3k — Cell Type Annotation",
       subtitle = "Labels based on canonical marker genes") +
  theme_bw(base_size = 11)

ggsave("plots/09_umap_annotated.png", p_umap_annotated,
       width = 9, height = 7, dpi = 150)
cat("  Saved: plots/09_umap_annotated.png\n")

# Side-by-side: clusters vs cell types
p_side <- p_umap_clusters + p_umap_annotated +
  plot_annotation(title = "PBMC 3k: clusters (left) vs cell type annotation (right)")
ggsave("plots/10_umap_comparison.png", p_side, width = 16, height = 7, dpi = 150)
cat("  Saved: plots/10_umap_comparison.png\n\n")


# =============================================================================
# STEP 11 — Save Results
# =============================================================================

cat("STEP 11: Saving results...\n")

# Cell counts per annotated cell type
cell_counts <- as.data.frame(table(pbmc$cell_type))
colnames(cell_counts) <- c("cell_type", "n_cells")
cell_counts <- cell_counts[order(-cell_counts$n_cells), ]
write.csv(cell_counts, "results/cell_type_counts.csv", row.names = FALSE)

# All markers with p-values
write.csv(markers, "results/cluster_markers_all.csv", row.names = FALSE)

# Top 5 markers per cluster
write.csv(top5, "results/cluster_markers_top5.csv", row.names = FALSE)

# Save Seurat object for future use
saveRDS(pbmc, "results/pbmc3k_seurat.rds")

cat("  Saved to results/\n\n")


# =============================================================================
# SUMMARY
# =============================================================================

cat("=== SUMMARY ===\n")
cat("Cells after QC:   ", ncol(pbmc), "\n")
cat("Genes in object:  ", nrow(pbmc), "\n")
cat("Clusters found:   ", n_clusters, "\n")
cat("Cell types annotated:\n")
print(cell_counts, row.names = FALSE)
cat("\nOutput files:\n")
cat("  plots/01_qc_violin.png\n")
cat("  plots/02_qc_scatter.png\n")
cat("  plots/03_variable_features.png\n")
cat("  plots/04_elbow_plot.png\n")
cat("  plots/05_pca_heatmap.png\n")
cat("  plots/06_umap_clusters.png\n")
cat("  plots/07_marker_heatmap.png\n")
cat("  plots/08_marker_dotplot.png\n")
cat("  plots/09_umap_annotated.png\n")
cat("  plots/10_umap_comparison.png\n")
cat("  results/cell_type_counts.csv\n")
cat("  results/cluster_markers_all.csv\n")
cat("  results/cluster_markers_top5.csv\n")
cat("  results/pbmc3k_seurat.rds\n")
