## ============================================================
## Week 6 Day 1 — GSE136103: Doublet detection (scDblFinder)
## Loads the QC-filtered, annotated Seurat object from script 04.
## Runs scDblFinder per GSM library, removes doublets, then
## re-runs the full downstream pipeline (normalise → PCA →
## cluster → UMAP → markers → annotation) on the clean data.
## Key check: do TREM2, SPP1, COL1A1 assignments change?
## ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(scDblFinder)
  library(SingleCellExperiment)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

# Parameters — identical to script 03
N_VARIABLE <- 2000
N_PCS      <- 20
RESOLUTION <- 0.5
SEED       <- 42
set.seed(SEED)

if (interactive()) {
  tryCatch({
    setwd(dirname(dirname(rstudioapi::getSourceEditorContext()$path)))
  }, error = function(e) invisible(NULL))
}
cat("Working directory:", getwd(), "\n")
dir.create("results", showWarnings = FALSE)
dir.create("plots",   showWarnings = FALSE)

# ============================================================
# 1. Load existing annotated object
# ============================================================
cat("\n=== 1. Loading annotated Seurat object ===\n")
obj <- readRDS("results/nafld_seurat_annotated.rds")
cells_pre <- ncol(obj)
cat(sprintf("Loaded: %d cells across %d genes\n", cells_pre, nrow(obj)))
cat("Cell types present:\n")
print(sort(table(obj$cell_type), decreasing = TRUE))

# ============================================================
# 2. Run scDblFinder per GSM library (per-sample detection)
# ============================================================
cat("\n=== 2. Running scDblFinder per GSM library ===\n")
cat("Converting to SingleCellExperiment...\n")
sce <- as.SingleCellExperiment(obj)

# Per-library detection prevents cross-sample correlations from inflating
# doublet scores — each library's doublets are estimated independently.
set.seed(SEED)
sce <- scDblFinder(sce, samples = "gsm_id", verbose = TRUE)

# Transfer predictions back to Seurat metadata
obj$dbl_class <- sce$scDblFinder.class
obj$dbl_score <- sce$scDblFinder.score
rm(sce); gc()

# Overall doublet counts
n_dbl      <- sum(obj$dbl_class == "doublet")
n_singlet  <- sum(obj$dbl_class == "singlet")
pct_dbl    <- round(100 * n_dbl / cells_pre, 2)

cat(sprintf("\nDoublets:  %d (%.2f%%)\n", n_dbl, pct_dbl))
cat(sprintf("Singlets:  %d (%.2f%%)\n", n_singlet, 100 - pct_dbl))

# Per-sample breakdown
dbl_per_sample <- obj@meta.data %>%
  group_by(gsm_id, donor, group, fraction) %>%
  summarise(
    total   = n(),
    dbl_n   = sum(dbl_class == "doublet"),
    dbl_pct = round(100 * dbl_n / total, 2),
    .groups = "drop"
  ) %>%
  arrange(desc(dbl_pct))

cat("\nDoublet rate per GSM library:\n")
print(as.data.frame(dbl_per_sample))
write.csv(dbl_per_sample, "results/doublet_summary_per_sample.csv", row.names = FALSE)

# Doublet rate per existing cell type
dbl_per_ct <- obj@meta.data %>%
  group_by(cell_type) %>%
  summarise(
    total   = n(),
    dbl_n   = sum(dbl_class == "doublet"),
    dbl_pct = round(100 * dbl_n / total, 2),
    .groups = "drop"
  ) %>%
  arrange(desc(dbl_pct))

cat("\nDoublet rate per cell type (pre-removal annotation):\n")
print(as.data.frame(dbl_per_ct))
write.csv(dbl_per_ct, "results/doublet_rate_per_celltype.csv", row.names = FALSE)

# UMAP: doublets overlaid on original embedding
p_dbl_umap <- DimPlot(
  obj, group.by = "dbl_class", pt.size = 0.3, order = "doublet",
  cols = c(singlet = "grey85", doublet = "#E65100")
) +
  ggtitle(sprintf("scDblFinder predictions — %d doublets (%.2f%%) on %d cells",
                  n_dbl, pct_dbl, cells_pre)) +
  theme(legend.position = "right")

ggsave("plots/umap_doublets.png", p_dbl_umap, width = 9, height = 7, dpi = 150)
cat("Saved: plots/umap_doublets.png\n")

# ============================================================
# 3. Remove doublets
# ============================================================
cat("\n=== 3. Removing doublets ===\n")
obj_clean <- subset(obj, subset = dbl_class == "singlet")
cells_post <- ncol(obj_clean)
cat(sprintf("Cells before: %d  |  After doublet removal: %d  |  Removed: %d\n",
            cells_pre, cells_post, cells_pre - cells_post))

# ============================================================
# 4. Re-run downstream pipeline on clean data
# ============================================================
cat("\n=== 4. Re-running pipeline on doublet-free data ===\n")

# Seurat v5: join layers if needed before any operations
obj_clean <- JoinLayers(obj_clean)

cat("Normalising...\n")
obj_clean <- NormalizeData(obj_clean, normalization.method = "LogNormalize",
                           scale.factor = 10000, verbose = FALSE)
cat("Variable features...\n")
obj_clean <- FindVariableFeatures(obj_clean, selection.method = "vst",
                                  nfeatures = N_VARIABLE, verbose = FALSE)
cat("Scaling...\n")
obj_clean <- ScaleData(obj_clean, features = VariableFeatures(obj_clean), verbose = FALSE)
cat("PCA...\n")
obj_clean <- RunPCA(obj_clean, npcs = 50, seed.use = SEED, verbose = FALSE)
cat("Clustering...\n")
obj_clean <- FindNeighbors(obj_clean, dims = 1:N_PCS, verbose = FALSE)
obj_clean <- FindClusters(obj_clean, resolution = RESOLUTION,
                          random.seed = SEED, verbose = FALSE)
cat("UMAP...\n")
obj_clean <- RunUMAP(obj_clean, dims = 1:N_PCS, seed.use = SEED, verbose = FALSE)

n_clusters_post <- length(levels(obj_clean$seurat_clusters))
cat(sprintf("Clusters (original): 21  |  Clusters (post-doublet): %d\n",
            n_clusters_post))

p_umap_post <- DimPlot(obj_clean, reduction = "umap", label = TRUE,
                        pt.size = 0.4, label.size = 3) +
  ggtitle(sprintf("%d clusters after doublet removal (res=%.1f, n=%d cells)",
                  n_clusters_post, RESOLUTION, cells_post))
ggsave("plots/umap_clusters_postdbl.png", p_umap_post, width = 9, height = 7, dpi = 150)
cat("Saved: plots/umap_clusters_postdbl.png\n")

# ============================================================
# 5. Find marker genes on clean data
# ============================================================
cat("\n=== 5. FindAllMarkers on doublet-free data ===\n")
markers_post <- FindAllMarkers(
  obj_clean,
  only.pos        = TRUE,
  min.pct         = 0.25,
  logfc.threshold = 0.25,
  test.use        = "wilcox",
  verbose         = FALSE
)

top5_post <- markers_post %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5) %>%
  ungroup()

write.csv(markers_post, "results/all_markers_postdbl.csv",       row.names = FALSE)
write.csv(top5_post,    "results/top5_markers_postdbl.csv",      row.names = FALSE)
cat(sprintf("Marker rows (original): 12,313  |  Post-doublet: %d\n", nrow(markers_post)))

cat("\nTop 3 markers per cluster (post-doublet):\n")
print(as.data.frame(
  top5_post %>% group_by(cluster) %>% slice_head(n = 3) %>%
    dplyr::select(cluster, gene, avg_log2FC)
))

# ============================================================
# 6. Re-annotate cell types using same marker logic
# ============================================================
cat("\n=== 6. Re-annotating cell types ===\n")

# Cell-type marker sets (same logic as script 04)
ct_marker_sets <- list(
  "CD4+ T cells"              = c("CD3D","CD3G","CD4","IL7R","LTB","CD40LG"),
  "CD8+ T cells"              = c("CD8A","CD8B","CD3G"),
  "CD8+ T cells (exhausted)"  = c("CD8A","LAG3","CRTAM","TIGIT","HAVCR2"),
  "Naive T cells"             = c("MAL","LEF1","CCR7","TCF7","SELL"),
  "NK/NKT cells"              = c("XCL1","TOX2","KLRC1","KLRC2"),
  "NK cells (cytotoxic)"      = c("GNLY","GZMB","FGFBP2","NKG7","PRF1"),
  "NK cells (liver-resident)" = c("IL2RB","CD160","NCR1","CXCR6"),
  "gdT/NK cells"              = c("GZMH","TRGC2","CX3CR1","FCGR3A"),
  "Monocytes"                 = c("S100A8","S100A9","S100A12","FCN1","LYZ"),
  "Dendritic cells"           = c("CD1C","FCER1A","CLEC10A","HLA-DQA1"),
  "Kupffer cells"             = c("C1QB","C1QC","CD5L","CD163","GPNMB","TIMD4"),
  "B cells"                   = c("CD79A","CD19","MS4A1","BANK1"),
  "Plasma cells"              = c("IGHGP","IGLL5","IGHA2","MZB1","JCHAIN"),
  "Endothelial cells"         = c("GPIHBP1","PODXL","AQP7","PECAM1"),
  "LSEC"                      = c("CLEC4G","FCN2","FCN3","OIT3","STAB2"),
  "Hepatic stellate cells"    = c("DCN","TCF21","RGS5","COL1A1","ACTA2"),
  "Hepatocytes"               = c("UGT2B15","UGT2A3","CYP3A4","ALB","APOB"),
  "Proliferating cells"       = c("CENPA","RRM2","TYMS","MKI67","TOP2A"),
  "Cholangiocytes"            = c("SCT","PTCRA","LRRC26","KRT7","KRT19"),
  "Mast cells"                = c("TPSAB1","TPSB2","CPA3","GATA2")
)

# Score each cluster against each cell-type marker set
clusters_post <- levels(obj_clean$seurat_clusters)
ct_annot_post <- setNames(character(length(clusters_post)), clusters_post)

for (cl in clusters_post) {
  cl_genes <- markers_post$gene[markers_post$cluster == cl &
                                markers_post$avg_log2FC >= 0.4]
  scores <- sapply(ct_marker_sets, function(m) sum(m %in% cl_genes))
  if (max(scores) == 0) {
    ct_annot_post[cl] <- paste0("Unknown_", cl)
  } else {
    ct_annot_post[cl] <- names(which.max(scores))
  }
}

cat("\nCluster → cell type (post-doublet):\n")
for (cl in names(ct_annot_post)) {
  cat(sprintf("  Cluster %2s → %s\n", cl, ct_annot_post[cl]))
}

obj_clean$cell_type_postdbl <- unname(ct_annot_post[as.character(obj_clean$seurat_clusters)])

Idents(obj_clean) <- "cell_type_postdbl"
p_annot_post <- DimPlot(obj_clean, reduction = "umap", label = TRUE,
                         label.size = 3, repel = TRUE, pt.size = 0.3) +
  ggtitle("UMAP — cell types after doublet removal") +
  theme(legend.position = "right")
ggsave("plots/umap_annotated_postdbl.png", p_annot_post,
       width = 14, height = 8, dpi = 150)
cat("Saved: plots/umap_annotated_postdbl.png\n")

# ============================================================
# 7. Key gene check — TREM2, SPP1, COL1A1
# ============================================================
cat("\n=== 7. Key gene assignments: TREM2 / SPP1 / COL1A1 ===\n")

check_genes <- c("TREM2", "SPP1", "COL1A1")

gene_results <- lapply(check_genes, function(gname) {
  hits <- markers_post[markers_post$gene == gname, ]
  if (nrow(hits) == 0) {
    cat(sprintf("  %-10s → not in FindAllMarkers output (below min.pct=0.25)\n", gname))
    return(data.frame(gene = gname, top_cluster = NA,
                      top_celltype = "below min.pct threshold",
                      avg_log2FC = NA, stringsAsFactors = FALSE))
  }
  top <- hits[which.max(hits$avg_log2FC), ]
  top_cl <- as.character(top$cluster)
  top_ct <- ct_annot_post[top_cl]
  cat(sprintf("  %-10s → cluster %s (%s)  log2FC=%.2f  pct.1=%.2f\n",
              gname, top_cl, top_ct, top$avg_log2FC, top$pct.1))
  data.frame(gene = gname, top_cluster = top_cl, top_celltype = top_ct,
             avg_log2FC = round(top$avg_log2FC, 3), stringsAsFactors = FALSE)
})
gene_results_df <- do.call(rbind, gene_results)

# Add original assignments for comparison
gene_results_df$original_celltype <- c(
  "Monocytes + Kupffer cells (within-celltype DE; not a FindAllMarkers hit — below min.pct=0.25)",
  "Hepatocytes (cluster 13 marker log2FC 6.24)",
  "Hepatic stellate cells (cluster 12 marker log2FC 7.89)"
)

write.csv(gene_results_df, "results/key_gene_comparison.csv", row.names = FALSE)
cat("\nKey gene comparison table:\n")
print(gene_results_df[, c("gene", "original_celltype", "top_celltype", "avg_log2FC")])

# Feature plots for the three key genes
Idents(obj_clean) <- "cell_type_postdbl"
p_feat <- FeaturePlot(obj_clean,
  features = c("TREM2", "SPP1", "COL1A1"),
  ncol = 3, pt.size = 0.2, order = TRUE
)
ggsave("plots/featureplot_key_genes_postdbl.png", p_feat,
       width = 18, height = 6, dpi = 150)
cat("Saved: plots/featureplot_key_genes_postdbl.png\n")

# ============================================================
# 8. Save clean object + summary
# ============================================================
saveRDS(obj_clean, "results/nafld_seurat_postdbl.rds")
cat("Saved: results/nafld_seurat_postdbl.rds\n")

cat("\n========================================================\n")
cat("DOUBLET DETECTION SUMMARY\n")
cat("========================================================\n")
cat(sprintf("Cells entering (post-QC):       %d\n", cells_pre))
cat(sprintf("Doublets detected (scDblFinder): %d  (%.2f%%)\n", n_dbl, pct_dbl))
cat(sprintf("Cells after doublet removal:     %d\n", cells_post))
cat(sprintf("Clusters — original:             21\n"))
cat(sprintf("Clusters — post-doublet:         %d\n", n_clusters_post))
cat("========================================================\n")
