#!/usr/bin/env Rscript
# =============================================================================
# Vu et al. 2025 MASLD Spatial Transcriptomics — Full Analysis Pipeline
# Week 8 Day 1 | Steps: Spot QC → SCTransform → PCA/UMAP/Clustering → CARD
# =============================================================================
# Data:      10x Visium CytAssist (FFPE), 8 arrays (from paper's GitHub code)
# Reference: GSE136103 human liver scRNA-seq (Ramachandran 2019, 35,050 cells)
# =============================================================================

suppressMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(hdf5r)
  library(patchwork)
  library(Matrix)
})

# ── Paths ─────────────────────────────────────────────────────────────────────
PROJ_ROOT <- "/Users/senaesen/Desktop/bioinfo-learning/week8-day1-masld-spatial-setup"
DATA_DIR  <- file.path(PROJ_ROOT, "data/raw")
REF_RDS   <- "/Users/senaesen/Desktop/bioinfo-learning/week6-day1-nafld-scrna-gse136103/results/nafld_seurat_annotated.rds"
RESULTS   <- file.path(PROJ_ROOT, "results")

dir.create(file.path(RESULTS, "qc_plots"),       recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(RESULTS, "spatial_domains"), recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(RESULTS, "deconvolution"),   recursive=TRUE, showWarnings=FALSE)

# ── Arrays (8 used in Vu et al. published analysis — VLP115_D and VLP116_A excluded) ──
ARRAYS <- c("VLP115_A", "VLP116_D", "VLP119_A", "VLP119_D",
            "VLP120_A", "VLP120_D", "VLP121_A", "VLP121_D")

# QC thresholds (matching paper's GitHub QC script logic)
MIN_GENES   <- 200   # minimum genes detected per spot
MAX_MITO_PC <- 10    # maximum mitochondrial gene percentage per spot
MIN_CELLS   <- 3     # gene must be detected in ≥3 spots to be retained

cat("==========================================================================\n")
cat("MASLD Spatial Transcriptomics Pipeline | Seurat", as.character(packageVersion("Seurat")), "\n")
cat("==========================================================================\n\n")

# =============================================================================
# STEP 1: Spot QC — load, filter, visualize per array
# =============================================================================
cat("=== STEP 1: Spot QC ===\n")

seurat_list <- list()
qc_summary  <- data.frame()

for (arr in ARRAYS) {
  cat("  Loading", arr, "...\n")

  arr_path <- file.path(DATA_DIR, arr)
  obj <- Load10X_Spatial(
    data     = arr_path,
    filename = "filtered_feature_bc_matrix.h5",
    slice    = arr
  )
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj$array_id <- arr

  pre_spots     <- ncol(obj)
  pre_med_genes <- median(obj$nFeature_Spatial)
  pre_med_umi   <- median(obj$nCount_Spatial)

  # Filter spots: minimum genes and maximum mitochondrial %
  obj <- subset(obj,
    subset = nFeature_Spatial >= MIN_GENES & percent.mt < MAX_MITO_PC
  )

  # Filter genes: keep genes detected in ≥3 spots
  raw_counts  <- GetAssayData(obj, assay="Spatial", layer="counts")
  genes_keep  <- rowSums(raw_counts > 0) >= MIN_CELLS
  obj         <- obj[genes_keep, ]

  post_spots <- ncol(obj)
  cat(sprintf("    %s: %d → %d spots (%.1f%% kept); %d genes retained\n",
      arr, pre_spots, post_spots, 100*post_spots/pre_spots, nrow(obj)))

  qc_summary <- rbind(qc_summary, data.frame(
    array          = arr,
    pre_spots      = pre_spots,
    post_spots     = post_spots,
    pct_kept       = round(100 * post_spots / pre_spots, 1),
    pre_med_genes  = pre_med_genes,
    pre_med_umi    = pre_med_umi,
    post_med_genes = median(obj$nFeature_Spatial),
    post_med_umi   = median(obj$nCount_Spatial),
    post_med_mito  = round(median(obj$percent.mt), 2)
  ))

  # Per-array QC violin plot
  p_qc <- VlnPlot(obj,
    features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"),
    ncol = 3, pt.size = 0
  ) & theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  p_qc <- p_qc + plot_annotation(title = paste(arr, "— Spot QC Metrics (post-filter)"))
  suppressWarnings(ggsave(
    file.path(RESULTS, "qc_plots", paste0(arr, "_qc_violin.png")),
    p_qc, width=10, height=4, dpi=150
  ))

  seurat_list[[arr]] <- obj
}

write.csv(qc_summary, file.path(RESULTS, "qc_plots", "qc_summary.csv"), row.names=FALSE)
cat("\nQC summary:\n")
print(qc_summary[, c("array","pre_spots","post_spots","pct_kept","post_med_genes","post_med_umi")])
cat("\n")

total_post_spots <- sum(qc_summary$post_spots)
cat(sprintf("Total spots after QC across all 8 arrays: %d\n\n", total_post_spots))

# =============================================================================
# STEP 2: SCTransform normalization per array, then merge
# =============================================================================
cat("=== STEP 2: SCTransform normalization ===\n")

seurat_norm <- lapply(names(seurat_list), function(arr) {
  cat("  SCTransform:", arr, "...\n")
  SCTransform(seurat_list[[arr]], assay="Spatial", verbose=FALSE,
              return.only.var.genes=FALSE)
})
names(seurat_norm) <- names(seurat_list)

cat("  Selecting shared variable features across arrays...\n")
var_feats <- SelectIntegrationFeatures(seurat_norm, nfeatures=3000)
cat(sprintf("  %d shared variable features selected\n", length(var_feats)))

cat("  Merging all 8 arrays...\n")
combined <- merge(
  seurat_norm[[1]],
  y            = seurat_norm[-1],
  add.cell.ids = names(seurat_norm),
  merge.dr     = FALSE
)
VariableFeatures(combined) <- var_feats

cat(sprintf("  Combined object: %d spots × %d genes\n\n", ncol(combined), nrow(combined)))

# =============================================================================
# STEP 3: Spatial Domain Discovery — PCA, UMAP, Clustering
# =============================================================================
cat("=== STEP 3: Spatial Domains — PCA, UMAP, Clustering ===\n")

set.seed(42)
combined <- RunPCA(combined,  assay="SCT", npcs=30, verbose=FALSE)
combined <- RunUMAP(combined, dims=1:20,  verbose=FALSE)
combined <- FindNeighbors(combined, dims=1:20, verbose=FALSE)
combined <- FindClusters(combined,  resolution=0.5, verbose=FALSE)

n_clusters <- length(unique(combined$seurat_clusters))
cat(sprintf("  Found %d spatial domain clusters (resolution=0.5)\n", n_clusters))

# UMAP colored by cluster and by array (batch check)
p_cluster <- DimPlot(combined, group.by="seurat_clusters", label=TRUE, label.size=4,
                     raster=FALSE, repel=TRUE) +
  ggtitle(sprintf("UMAP — %d Spatial Domains (res=0.5)", n_clusters)) +
  theme(legend.position="none")

p_array <- DimPlot(combined, group.by="array_id", raster=FALSE) +
  ggtitle("UMAP — By Visium Array") +
  theme(legend.text=element_text(size=8))

suppressWarnings(ggsave(
  file.path(RESULTS, "spatial_domains", "umap_clusters_and_arrays.png"),
  p_cluster | p_array, width=16, height=6, dpi=150
))

# Elbow plot for PCA dimension choice
p_elbow <- ElbowPlot(combined, ndims=30) + ggtitle("PCA Elbow Plot")
ggsave(file.path(RESULTS, "spatial_domains", "elbow_plot.png"),
       p_elbow, width=6, height=4, dpi=150)

# Cluster × array spot count table
cluster_tab <- as.data.frame(table(
  Cluster = combined$seurat_clusters,
  Array   = combined$array_id
))
cluster_tab$pct <- round(100 * cluster_tab$Freq / sum(cluster_tab$Freq), 2)
write.csv(cluster_tab, file.path(RESULTS, "spatial_domains", "cluster_by_array.csv"), row.names=FALSE)

# Cluster size summary
cluster_sizes <- as.data.frame(table(Cluster=combined$seurat_clusters))
cluster_sizes$pct <- round(100 * cluster_sizes$Freq / sum(cluster_sizes$Freq), 1)
cat("  Cluster sizes:\n")
print(cluster_sizes)

write.csv(cluster_sizes, file.path(RESULTS, "spatial_domains", "cluster_sizes.csv"), row.names=FALSE)

# Save combined object
saveRDS(combined, file.path(RESULTS, "masld_spatial_combined.rds"))
cat("\n  Saved: results/masld_spatial_combined.rds\n\n")

# =============================================================================
# STEP 4: Deconvolution
# =============================================================================
cat("=== STEP 4: Cell-type Deconvolution ===\n")

# Load scRNA-seq reference (GSE136103 — Ramachandran 2019, human liver)
cat("  Loading GSE136103 scRNA-seq reference...\n")
ref_obj <- readRDS(REF_RDS)
cat(sprintf("  Reference: %d genes × %d cells, %d cell types\n",
    nrow(ref_obj), ncol(ref_obj), length(unique(ref_obj$cell_type))))

if (requireNamespace("CARD", quietly=TRUE)) {
  # ── CARD deconvolution ──────────────────────────────────────────────────────
  cat("  Using CARD deconvolution (Ma & Zhou 2022, Nat Biotechnol)\n\n")
  library(CARD)

  sc_count <- GetAssayData(ref_obj, assay="RNA", layer="counts")
  sc_meta  <- data.frame(
    cellID     = colnames(ref_obj),
    cellType   = ref_obj$cell_type,
    sampleInfo = ref_obj$orig.ident,
    row.names  = colnames(ref_obj)
  )

  all_props <- list()

  for (arr in ARRAYS) {
    cat("  Processing", arr, "...\n")
    arr_obj <- seurat_list[[arr]]

    spatial_count <- GetAssayData(arr_obj, assay="Spatial", layer="counts")

    # Shared genes: reference ∩ spatial
    shared_genes <- intersect(rownames(sc_count), rownames(spatial_count))
    cat(sprintf("    Shared genes: %d\n", length(shared_genes)))

    # Spatial coordinates (Seurat v5: GetTissueCoordinates returns x/y or imagerow/imagecol)
    coord_df <- tryCatch(
      GetTissueCoordinates(arr_obj),
      error = function(e) {
        # Fallback: extract from Spatial assay metadata
        as.data.frame(arr_obj@images[[arr]]@coordinates[, c("imagerow","imagecol")])
      }
    )
    spatial_loc <- data.frame(
      x = coord_df[, 1],
      y = coord_df[, 2],
      row.names = rownames(coord_df)
    )

    # Align barcodes between count matrix and coordinate table
    common_bc     <- intersect(colnames(spatial_count), rownames(spatial_loc))
    spatial_count_sub <- spatial_count[shared_genes, common_bc]
    spatial_loc_sub   <- spatial_loc[common_bc, ]

    tryCatch({
      CARD_obj <- createCARDObject(
        sc_count         = sc_count[shared_genes, ],
        sc_meta          = sc_meta,
        spatial_count    = spatial_count_sub,
        spatial_location = spatial_loc_sub,
        ct.varname       = "cellType",
        ct.select        = NULL,
        sample.varname   = "sampleInfo",
        minCountGene     = 100,
        minCountSpot     = 5
      )
      CARD_obj   <- CARD_deconvolution(CARD_obj)
      props      <- as.data.frame(CARD_obj@Proportion_CARD)
      props$barcode  <- rownames(props)
      props$array_id <- arr
      write.csv(props,
        file.path(RESULTS, "deconvolution", paste0(arr, "_CARD_proportions.csv")),
        row.names=FALSE
      )
      all_props[[arr]] <- props
      cat(sprintf("    Saved CARD proportions for %s (%d spots)\n", arr, nrow(props)))
    }, error = function(e) {
      cat(sprintf("    CARD error for %s: %s\n", arr, conditionMessage(e)))
    })
  }

  if (length(all_props) > 0) {
    combined_props <- do.call(rbind, all_props)
    write.csv(combined_props,
      file.path(RESULTS, "deconvolution", "all_arrays_CARD_proportions.csv"),
      row.names=FALSE
    )

    # Summary: mean cell-type proportion per array
    ct_cols <- setdiff(colnames(combined_props), c("barcode","array_id"))
    mean_props <- combined_props %>%
      group_by(array_id) %>%
      summarise(across(all_of(ct_cols), mean), .groups="drop")
    write.csv(mean_props,
      file.path(RESULTS, "deconvolution", "mean_celltype_props_per_array.csv"),
      row.names=FALSE
    )

    # Heatmap of mean proportions by array
    mat <- as.matrix(mean_props[, ct_cols])
    rownames(mat) <- mean_props$array_id
    # Reshape for ggplot
    df_long <- tidyr::pivot_longer(
      as.data.frame(mat) %>% mutate(Array=rownames(mat)),
      cols = -Array, names_to="CellType", values_to="Proportion"
    )
    p_heat <- ggplot(df_long, aes(CellType, Array, fill=Proportion)) +
      geom_tile() +
      scale_fill_gradientn(colors=c("white","#2196F3","#1565C0")) +
      theme_bw(base_size=10) +
      theme(axis.text.x = element_text(angle=45, hjust=1, size=7)) +
      labs(title="Mean CARD Cell-Type Proportions per Array",
           x=NULL, y=NULL)
    suppressWarnings(ggsave(
      file.path(RESULTS, "deconvolution", "card_proportions_heatmap.png"),
      p_heat, width=12, height=5, dpi=150
    ))
    cat("\n  CARD deconvolution complete. Results saved.\n")
  }

} else {
  # ── Fallback: Seurat AddModuleScore cell-type scoring ──────────────────────
  cat("  CARD not available — using Seurat AddModuleScore for cell-type scoring\n\n")

  # Build marker gene lists per cell type from reference
  Seurat::DefaultAssay(ref_obj) <- "RNA"
  ref_obj <- NormalizeData(ref_obj, verbose=FALSE)

  # Use top 50 genes by average expression per cell type as marker sets
  ref_counts_norm <- GetAssayData(ref_obj, assay="RNA", layer="data")

  marker_genes <- lapply(unique(ref_obj$cell_type), function(ct) {
    ct_cells <- colnames(ref_obj)[ref_obj$cell_type == ct]
    if (length(ct_cells) < 5) return(NULL)
    avg_exp <- rowMeans(ref_counts_norm[, ct_cells, drop=FALSE])
    names(sort(avg_exp, decreasing=TRUE))[1:min(50, length(avg_exp))]
  })
  names(marker_genes) <- unique(ref_obj$cell_type)
  marker_genes <- Filter(Negate(is.null), marker_genes)

  # Score each cell type in combined spatial object
  # Use SCT data for scoring
  DefaultAssay(combined) <- "SCT"

  for (ct_name in names(marker_genes)) {
    module_name <- gsub("[^A-Za-z0-9]", "_", ct_name)
    markers_avail <- intersect(marker_genes[[ct_name]], rownames(combined))
    if (length(markers_avail) < 5) next
    combined <- AddModuleScore(combined,
      features = list(markers_avail),
      name     = paste0("CT_", module_name),
      verbose  = FALSE
    )
  }

  # Save spot-level scores
  score_cols <- grep("^CT_", colnames(combined@meta.data), value=TRUE)
  ct_scores  <- combined@meta.data[, c("array_id", "seurat_clusters", score_cols)]
  ct_scores$barcode <- rownames(ct_scores)
  write.csv(ct_scores,
    file.path(RESULTS, "deconvolution", "celltype_module_scores.csv"),
    row.names=FALSE
  )

  # Mean scores per cluster
  mean_scores <- ct_scores %>%
    group_by(seurat_clusters) %>%
    summarise(across(all_of(score_cols), mean), .groups="drop")
  write.csv(mean_scores,
    file.path(RESULTS, "deconvolution", "mean_celltype_scores_per_cluster.csv"),
    row.names=FALSE
  )

  # Heatmap of cell-type scores per cluster
  mat <- as.matrix(mean_scores[, score_cols])
  rownames(mat) <- paste0("C", mean_scores$seurat_clusters)
  colnames(mat) <- gsub("^CT_", "", colnames(mat))
  colnames(mat) <- gsub("_1$", "", colnames(mat))

  df_long <- as.data.frame(mat) %>%
    mutate(Cluster=rownames(mat)) %>%
    tidyr::pivot_longer(-Cluster, names_to="CellType", values_to="Score")

  p_heat <- ggplot(df_long, aes(CellType, Cluster, fill=Score)) +
    geom_tile() +
    scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=0) +
    theme_bw(base_size=9) +
    theme(axis.text.x=element_text(angle=45, hjust=1, size=7)) +
    labs(title="Mean Cell-Type Module Scores per Spatial Cluster",
         x=NULL, y="Cluster")
  suppressWarnings(ggsave(
    file.path(RESULTS, "deconvolution", "celltype_scores_heatmap.png"),
    p_heat, width=14, height=6, dpi=150
  ))

  cat("  Cell-type module scores saved.\n")
  cat("  Note: AddModuleScore gives relative cell-type activity scores per spot,\n")
  cat("        not compositional proportions. Install CARD for true deconvolution.\n")
}

# =============================================================================
# Save final combined object (with deconvolution scores if AddModuleScore was used)
# =============================================================================
saveRDS(combined, file.path(RESULTS, "masld_spatial_combined.rds"))
cat("\n  Final Seurat object saved: results/masld_spatial_combined.rds\n")

# =============================================================================
# Session info
# =============================================================================
cat("\n=== Pipeline complete ===\n")
cat(sprintf("Results directory: %s\n", RESULTS))
cat("\nKey outputs:\n")
cat("  results/qc_plots/qc_summary.csv           — per-array spot QC stats\n")
cat("  results/spatial_domains/umap_clusters*.png — UMAP + cluster plots\n")
cat("  results/spatial_domains/cluster_sizes.csv  — cluster size table\n")
cat("  results/deconvolution/                      — cell-type composition\n")
cat("  results/masld_spatial_combined.rds          — final Seurat object\n")
cat("\nSession info:\n")
sessionInfo()
