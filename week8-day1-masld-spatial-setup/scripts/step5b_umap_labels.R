#!/usr/bin/env Rscript
# Step 5b — Add biological labels to Seurat object and regenerate UMAP.
# Runs after step5_biological_annotation.R has written cluster_biological_annotation.csv.
# The only change: use unname() to prevent Seurat v5 named-vector assignment bug.

suppressMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

PROJ_ROOT <- "/Users/senaesen/Desktop/bioinfo-learning/week8-day1-masld-spatial-setup"
RESULTS   <- file.path(PROJ_ROOT, "results")
DOM_DIR   <- file.path(RESULTS, "spatial_domains")

cat("Loading Seurat object...\n")
combined <- readRDS(file.path(RESULTS, "masld_spatial_combined.rds"))
cat(sprintf("  %d spots, %d clusters\n", ncol(combined),
    length(unique(combined$seurat_clusters))))

cat("Loading annotation table...\n")
ann <- read.csv(file.path(DOM_DIR, "cluster_biological_annotation.csv"),
                stringsAsFactors = FALSE)
cat(sprintf("  %d cluster annotations loaded\n", nrow(ann)))

# Build lookup vector, then strip names to avoid Seurat v5's named-vector matching
cl_to_label <- setNames(
  paste0("C", ann$cluster_id, ": ", ann$biological_label),
  as.character(ann$cluster_id)
)
cl_to_short <- setNames(
  ann$biological_label,
  as.character(ann$cluster_id)
)

cat("Adding bio_label metadata...\n")
combined$bio_label       <- unname(cl_to_label[as.character(combined$seurat_clusters)])
combined$bio_label_short <- unname(cl_to_short[as.character(combined$seurat_clusters)])
cat(sprintf("  Unique labels: %d\n", length(unique(combined$bio_label_short))))

cat("Generating UMAP panels...\n")

p_clusters <- DimPlot(combined,
  group.by   = "seurat_clusters",
  label      = TRUE, label.size = 3.5,
  raster     = FALSE, repel      = TRUE) +
  ggtitle("24 Spatial Clusters (numbered)") +
  theme(legend.position = "none")

p_bio <- DimPlot(combined,
  group.by   = "bio_label_short",
  label      = TRUE, label.size = 2.8,
  raster     = FALSE, repel      = TRUE) +
  ggtitle("Biological Annotation") +
  theme(
    legend.position  = "right",
    legend.text      = element_text(size = 6.5),
    legend.key.size  = unit(0.28, "cm"),
    legend.title     = element_blank()
  )

suppressWarnings(ggsave(
  file.path(DOM_DIR, "umap_biological_labels.png"),
  p_clusters | p_bio, width = 22, height = 7, dpi = 150
))
cat("  Saved: umap_biological_labels.png\n")

suppressWarnings(ggsave(
  file.path(DOM_DIR, "umap_bio_labels_only.png"),
  p_bio, width = 14, height = 7, dpi = 150
))
cat("  Saved: umap_bio_labels_only.png\n")

saveRDS(combined, file.path(RESULTS, "masld_spatial_combined.rds"))
cat("  Updated Seurat object saved.\n")
cat("Done.\n")
