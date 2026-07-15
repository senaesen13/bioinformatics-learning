## Week 6 Day 1 — GSE136103: Seurat pipeline (NAFLD subset)
## 5 healthy liver donors + Cirrhotic1 (NAFLD) + Cirrhotic4 (NAFLD)
## Stops at clustering + marker genes; no cell-type annotation yet.

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

# ── Parameters ────────────────────────────────────────────────────────────────
MIN_FEATURES <- 200
MAX_FEATURES <- 2500
MAX_MT_PCT   <- 5
N_VARIABLE   <- 2000
N_PCS        <- 20    # chosen after reviewing elbow plot
RESOLUTION   <- 0.5
SEED         <- 42

# ── Load NAFLD subset metadata ────────────────────────────────────────────────
meta     <- read.csv("../results/nafld_subset_metadata.csv", stringsAsFactors = FALSE)
data_dir <- "../data"

# ── Helper: read one GSM directory (non-standard filenames) ──────────────────
read_gsm <- function(gsm_id, meta_row) {
  dir_path <- file.path(data_dir, gsm_id)
  files    <- list.files(dir_path, full.names = TRUE)

  mtx_f  <- files[grepl("matrix\\.mtx\\.gz$",  files)]
  bar_f  <- files[grepl("barcodes\\.tsv\\.gz$", files)]
  gene_f <- files[grepl("genes\\.tsv\\.gz$",    files)]

  counts <- ReadMtx(
    mtx             = mtx_f,
    cells           = bar_f,
    features        = gene_f,
    feature.column  = 2,   # use gene symbol, not Ensembl ID
    unique.features = TRUE
  )

  obj <- CreateSeuratObject(
    counts       = counts,
    project      = meta_row$donor,
    min.cells    = 3,
    min.features = MIN_FEATURES
  )
  obj$gsm_id         <- gsm_id
  obj$donor          <- meta_row$donor
  obj$group          <- meta_row$group
  obj$etiology       <- meta_row$etiology
  obj$fraction       <- meta_row$fraction
  obj$disease_status <- ifelse(meta_row$group == "Healthy liver",
                               "healthy", "NAFLD_cirrhosis")
  obj
}

# ── 1. Load all 15 libraries ─────────────────────────────────────────────────
message("\n=== Step 1: Loading ", nrow(meta), " libraries ===\n")
seurat_list <- vector("list", nrow(meta))
for (i in seq_len(nrow(meta))) {
  row <- meta[i, ]
  message("  Loading ", row$geo_accession,
          " [", row$donor, " | ", row$fraction, "]")
  seurat_list[[i]] <- read_gsm(row$geo_accession, row)
  message("    -> ", ncol(seurat_list[[i]]), " cells")
}
names(seurat_list) <- meta$geo_accession

# ── 2. Merge into one object ─────────────────────────────────────────────────
message("\n=== Step 2: Merging ===\n")
merged <- merge(
  x          = seurat_list[[1]],
  y          = seurat_list[-1],
  add.cell.ids = meta$geo_accession
)
cells_before_qc <- ncol(merged)
message("Merged: ", cells_before_qc, " cells x ", nrow(merged), " features")

# Seurat v5: merge keeps per-sample layers separate; join before any DE testing
merged <- JoinLayers(merged)

# ── 3. QC ─────────────────────────────────────────────────────────────────────
message("\n=== Step 3: QC ===\n")
merged[["percent.mt"]] <- PercentageFeatureSet(merged, pattern = "^MT-")

message("QC summary BEFORE filtering:")
cat("nFeature_RNA: "); print(summary(merged$nFeature_RNA))
cat("percent.mt  : "); print(summary(merged$percent.mt))

p_pre <- VlnPlot(
  merged,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  group.by = "donor", ncol = 3, pt.size = 0
) + plot_annotation(title = "QC before filtering — by donor")
ggsave("../plots/qc_violin_before.png", p_pre, width = 16, height = 6, dpi = 150)

merged <- subset(
  merged,
  subset = nFeature_RNA >= MIN_FEATURES &
           nFeature_RNA <= MAX_FEATURES &
           percent.mt   <  MAX_MT_PCT
)
cells_after_qc <- ncol(merged)

message("Cells before QC : ", cells_before_qc)
message("Cells after QC  : ", cells_after_qc,
        " (", round(cells_after_qc / cells_before_qc * 100, 1), "% retained)")

p_post <- VlnPlot(
  merged,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  group.by = "donor", ncol = 3, pt.size = 0
) + plot_annotation(title = "QC after filtering — by donor")
ggsave("../plots/qc_violin_after.png", p_post, width = 16, height = 6, dpi = 150)

# ── 4. Normalise + variable features ─────────────────────────────────────────
message("\n=== Step 4: Normalise + variable features ===\n")
merged <- NormalizeData(merged, normalization.method = "LogNormalize",
                        scale.factor = 10000, verbose = FALSE)
merged <- FindVariableFeatures(merged, selection.method = "vst",
                               nfeatures = N_VARIABLE, verbose = FALSE)

top10 <- head(VariableFeatures(merged), 10)
message("Top 10 variable features: ", paste(top10, collapse = ", "))

p_hvg <- LabelPoints(VariableFeaturePlot(merged), points = top10, repel = TRUE)
ggsave("../plots/variable_features.png", p_hvg, width = 10, height = 6, dpi = 150)

# ── 5. Scale + PCA ───────────────────────────────────────────────────────────
message("\n=== Step 5: Scale + PCA ===\n")
# Scale variable features only (faster; DoHeatmap uses var features anyway)
merged <- ScaleData(merged, features = VariableFeatures(merged), verbose = FALSE)
merged <- RunPCA(merged, npcs = 50, seed.use = SEED, verbose = FALSE)

p_elbow <- ElbowPlot(merged, ndims = 50) +
  geom_vline(xintercept = N_PCS, linetype = "dashed", colour = "red") +
  ggtitle(paste0("Elbow plot — red line = N_PCS used (", N_PCS, ")"))
ggsave("../plots/pca_elbow.png", p_elbow, width = 8, height = 5, dpi = 150)
message("Elbow plot saved. Using N_PCS = ", N_PCS)

p_pca <- DimPlot(merged, reduction = "pca", group.by = "donor") +
  ggtitle("PCA coloured by donor")
ggsave("../plots/pca_donor.png", p_pca, width = 9, height = 7, dpi = 150)

# ── 6. Cluster + UMAP ────────────────────────────────────────────────────────
message("\n=== Step 6: Cluster + UMAP (dims 1:", N_PCS,
        ", res ", RESOLUTION, ") ===\n")
set.seed(SEED)
merged <- FindNeighbors(merged, dims = 1:N_PCS, verbose = FALSE)
merged <- FindClusters(merged, resolution = RESOLUTION,
                       random.seed = SEED, verbose = FALSE)
merged <- RunUMAP(merged, dims = 1:N_PCS, seed.use = SEED, verbose = FALSE)

n_clusters <- length(levels(merged$seurat_clusters))
message("Clusters found: ", n_clusters)

p_clust   <- DimPlot(merged, reduction = "umap", label = TRUE, pt.size = 0.5) +
  ggtitle(paste0(n_clusters, " clusters (res=", RESOLUTION, ")"))
p_disease <- DimPlot(merged, reduction = "umap", group.by = "disease_status",
                     pt.size = 0.4) + ggtitle("Disease status")
p_donor   <- DimPlot(merged, reduction = "umap", group.by = "donor",
                     pt.size = 0.4) + ggtitle("Donor")
p_frac    <- DimPlot(merged, reduction = "umap", group.by = "fraction",
                     pt.size = 0.4) + ggtitle("CD45 fraction")

ggsave("../plots/umap_clusters.png",  p_clust,   width = 9, height = 7, dpi = 150)
ggsave("../plots/umap_disease.png",   p_disease, width = 9, height = 7, dpi = 150)
ggsave("../plots/umap_donor.png",     p_donor,   width = 10, height = 7, dpi = 150)
ggsave("../plots/umap_fraction.png",  p_frac,    width = 9,  height = 7, dpi = 150)

p_panel <- (p_clust | p_disease) / (p_donor | p_frac)
ggsave("../plots/umap_panel.png", p_panel, width = 16, height = 14, dpi = 150)

# ── 7. Save checkpoint BEFORE markers (markers step is slow) ─────────────────
saveRDS(merged, "../results/nafld_seurat_clustered.rds")
message("\nCheckpoint saved: results/nafld_seurat_clustered.rds")
message("Run script 04_find_markers.R next to compute marker genes.\n")

message("========================================")
message("PIPELINE SUMMARY (through clustering)")
message("========================================")
message("Libraries loaded        : ", nrow(meta))
message("Donors                  : 7 (5 healthy | 2 NAFLD cirrhotic)")
message("Cells before QC         : ", cells_before_qc)
message("Cells after QC          : ", cells_after_qc,
        " (", round(cells_after_qc / cells_before_qc * 100, 1), "% retained)")
message("Clusters found          : ", n_clusters)
message("========================================\n")

message("Cells per cluster:")
print(sort(table(merged$seurat_clusters)))
message("\nCells per donor:")
print(sort(table(merged$donor), decreasing = TRUE))

# ── 8. Marker genes (split to separate script for long run) ──────────────────
message("\n=== Step 7: FindAllMarkers ===\n")
all_markers <- FindAllMarkers(
  merged,
  only.pos        = TRUE,
  min.pct         = 0.25,
  logfc.threshold = 0.25,
  test.use        = "wilcox",
  verbose         = FALSE
)

top5 <- all_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5) %>%
  ungroup()

write.csv(all_markers, "../results/all_markers.csv",            row.names = FALSE)
write.csv(top5,        "../results/top5_markers_per_cluster.csv", row.names = FALSE)

message("\nTop 5 markers per cluster (by avg_log2FC):\n")
print(as.data.frame(top5[, c("cluster","gene","avg_log2FC","pct.1","pct.2","p_val_adj")]),
      row.names = FALSE)

top5_genes <- unique(top5$gene)
p_heat <- DoHeatmap(merged, features = top5_genes, size = 3) +
  ggtitle("Top 5 markers per cluster") +
  theme(axis.text.y = element_text(size = 7))
ggsave("../plots/marker_heatmap.png", p_heat,
       width = 16, height = max(10, ceiling(length(top5_genes) * 0.3) + 4),
       dpi = 150)

# ── 8. Save Seurat object ─────────────────────────────────────────────────────
saveRDS(merged, "../results/nafld_seurat_clustered.rds")
message("\nSaved: results/nafld_seurat_clustered.rds")

# ── 9. Summary report ─────────────────────────────────────────────────────────
message("\n========================================")
message("PIPELINE SUMMARY")
message("========================================")
message("Libraries loaded        : ", nrow(meta))
message("Donors                  : 7 (5 healthy | 2 NAFLD cirrhotic)")
message("Cells before QC         : ", cells_before_qc)
message("Cells after QC          : ", cells_after_qc,
        " (", round(cells_after_qc / cells_before_qc * 100, 1), "% retained)")
message("QC thresholds           : nFeature 200–2500, percent.mt < 5%")
message("Variable features       : ", N_VARIABLE)
message("PCs used                : ", N_PCS)
message("Clustering resolution   : ", RESOLUTION)
message("Clusters found          : ", n_clusters)
message("Marker gene rows saved  : ", nrow(all_markers))
message("========================================")
message("\nCells per cluster:")
print(sort(table(merged$seurat_clusters)))
message("\nCells per donor:")
print(sort(table(merged$donor), decreasing = TRUE))
