## ============================================================
## Week 4 Day 3 — scRNA-seq Human Liver: Healthy vs Cirrhosis
## Dataset   : GSE136103 (Ramachandran et al. 2019, Nature)
## Platform  : GPL20301 (Human, HiSeq 4000) — ignoring GPL21103 (Mouse)
## Fractions : CD45+ (immune) and CD45- (parenchymal) liver fractions
##             Blood samples (blood1–4) excluded; liver only
## Cells     : Up to 20,000 (stratified random sample by condition)
## Pipeline  : QC → Normalise → PCA → UMAP → Clustering → Annotation
## Goal      : Compare macrophage, hepatocyte, stellate, endothelial, T cell
##             proportions between Healthy and Cirrhosis livers
## Run from  : week4-day3-nafld-scrna-seq/
##   Rscript scripts/nafld_scrna_liver.R
## ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(GEOquery)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(RColorBrewer)
  library(cowplot)
  library(scales)
  library(tidyr)
  library(Matrix)
})

set.seed(42)
options(timeout = 7200)

if (!interactive()) pdf(NULL)

dir.create("plots",   showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)
dir.create("data",    showWarnings = FALSE)

COND_COLS <- c(Healthy = "#2196F3", Cirrhosis = "#F44336")

CELLTYPE_COLS <- c(
  Macrophage    = "#E65100",
  Hepatocyte    = "#1B5E20",
  Stellate      = "#4A148C",
  Endothelial   = "#006064",
  "T Cell"      = "#B71C1C",
  Cholangiocyte = "#F57F17",
  "B Cell"      = "#880E4F",
  "NK Cell"     = "#01579B",
  Other         = "#9E9E9E"
)

cat("=== Week 4 Day 3: Human Liver scRNA-seq (GSE136103) ===\n\n")

## ============================================================
## STEP 1 — Download GEO metadata and identify human liver samples
## ============================================================
cat("STEP 1: Downloading GSE136103 metadata...\n")

gse <- getGEO("GSE136103", GSEMatrix = TRUE, destdir = "data/")

platforms <- sapply(gse, annotation)
cat("Platforms:", paste(platforms, collapse = ", "), "\n")

human_idx <- which(platforms == "GPL20301")
eset <- gse[[human_idx[1]]]
meta <- pData(eset)
cat("Total GPL20301 samples:", nrow(meta), "\n")

# Parse condition from characteristics columns
char_cols <- grep("characteristics", colnames(meta), value = TRUE, ignore.case = TRUE)

assign_condition <- function(row) {
  combined <- paste(unlist(row[char_cols]), collapse = " ")
  if (grepl("healthy|normal|unaffected|donor", combined, ignore.case = TRUE))   return("Healthy")
  if (grepl("cirrhosis|cirrhotic|fibrosis", combined, ignore.case = TRUE))       return("Cirrhosis")
  # Fall back to title
  title <- as.character(row[["title"]])
  if (grepl("healthy|normal|donor", title, ignore.case = TRUE))   return("Healthy")
  if (grepl("cirrhosis|cirrhotic", title, ignore.case = TRUE))    return("Cirrhosis")
  return("Unknown")
}

meta$condition <- apply(meta, 1, assign_condition)
meta$title_lc  <- tolower(meta$title)

# Exclude blood samples — we want liver only
meta$tissue <- ifelse(grepl("blood", meta$title_lc), "Blood", "Liver")
liver_meta  <- meta[meta$tissue == "Liver", ]

cat("Liver samples:", nrow(liver_meta), "\n")
cat("\nCondition distribution (liver only):\n")
print(table(liver_meta$condition))

write.csv(liver_meta[, c("title","geo_accession","condition","tissue")],
          "results/sample_metadata.csv", row.names = TRUE)

## ============================================================
## STEP 2 — Files already downloaded; just list them
## ============================================================
cat("\nSTEP 2: Checking downloaded supplementary files...\n")

gsm_ids <- rownames(liver_meta)

# Verify all liver GSM directories exist
for (gsm in gsm_ids) {
  n <- length(list.files(file.path("data", gsm)))
  cat("  ", gsm, "(", liver_meta[gsm, "condition"], "):", n, "files\n")
}

## ============================================================
## STEP 3 — Load count matrices with ReadMtx()
## ============================================================
cat("\nSTEP 3: Loading count matrices with ReadMtx()...\n")

# GSE136103 stores files as {GSM}_{samplename}_{fraction}_{matrix|barcodes|genes}.{tsv|mtx}
# Read10X() expects canonical filenames — use ReadMtx() instead.

load_prefixed_mtx <- function(gsm_dir, gsm_id, condition) {
  mat_files <- list.files(gsm_dir, pattern = "_matrix\\.mtx$", full.names = TRUE)
  if (length(mat_files) == 0) {
    # Try .mtx.gz if .mtx not yet unzipped
    mat_files <- list.files(gsm_dir, pattern = "_matrix\\.mtx\\.gz$", full.names = TRUE)
    if (length(mat_files) == 0) return(NULL)
  }
  mat_file <- mat_files[1]
  prefix   <- sub("_matrix\\.mtx(\\.gz)?$", "", mat_file)

  bar_file  <- paste0(prefix, "_barcodes.tsv")
  gene_file <- paste0(prefix, "_genes.tsv")
  if (!file.exists(bar_file))  bar_file  <- paste0(bar_file, ".gz")
  if (!file.exists(gene_file)) gene_file <- paste0(gene_file, ".gz")

  if (!file.exists(bar_file) || !file.exists(gene_file)) return(NULL)

  counts <- ReadMtx(
    mtx            = mat_file,
    cells          = bar_file,
    features       = gene_file,
    feature.column = 2   # col 1 = Ensembl ID, col 2 = gene symbol
  )
  counts
}

seurat_list <- list()

for (gsm in gsm_ids) {
  cond    <- liver_meta[gsm, "condition"]
  title   <- liver_meta[gsm, "title"]
  gsm_dir <- file.path("data", gsm)

  counts <- load_prefixed_mtx(gsm_dir, gsm, cond)
  if (is.null(counts)) {
    cat("  ", gsm, ": no matrix found — skipping\n")
    next
  }

  so <- CreateSeuratObject(
    counts       = counts,
    project      = gsm,
    min.cells    = 3,
    min.features = 200
  )
  so$sample     <- gsm
  so$condition  <- cond
  so$title      <- title
  # Donor ID: extract healthy1, healthy2, cirrhotic1 etc. from title
  so$donor      <- sub(".*(healthy\\d+|cirrhotic\\d+).*", "\\1", title,
                        ignore.case = TRUE, perl = TRUE)
  so$fraction   <- ifelse(grepl("cd45\\+", title, ignore.case = TRUE), "CD45+",
                   ifelse(grepl("cd45-",  title, ignore.case = TRUE), "CD45-", "Unknown"))

  seurat_list[[gsm]] <- so
  cat("  ", gsm, "(", cond, "—", so$fraction[1], "):",
      ncol(so), "cells,", nrow(so), "genes\n")
}

cat("\nLoaded", length(seurat_list), "samples successfully.\n")

if (length(seurat_list) == 0) stop("No samples loaded.")

## ============================================================
## STEP 4 — Merge all liver samples into one object
## ============================================================
cat("\nSTEP 4: Merging Seurat objects...\n")

liver <- merge(
  x          = seurat_list[[1]],
  y          = seurat_list[-1],
  add.cell.ids = names(seurat_list),
  project    = "GSE136103_HumanLiver"
)
liver <- JoinLayers(liver)

cat("Merged:", ncol(liver), "cells,", nrow(liver), "genes\n")
cat("Condition breakdown:\n")
print(table(liver$condition))

## ============================================================
## STEP 5 — QC filtering
## ============================================================
cat("\nSTEP 5: QC filtering...\n")

liver[["percent.mt"]] <- PercentageFeatureSet(liver, pattern = "^MT-")
liver[["percent.rb"]] <- PercentageFeatureSet(liver, pattern = "^RP[SL]")

p_qc_before <- VlnPlot(
  liver,
  features  = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  group.by  = "condition",
  pt.size   = 0,
  cols      = COND_COLS,
  ncol      = 3
) & theme(legend.position = "none", axis.title.x = element_blank())

ggsave("plots/01_qc_before_filtering.png", p_qc_before,
       width = 12, height = 5, dpi = 150)
cat("  Saved plots/01_qc_before_filtering.png\n")

cat("  Pre-filter:", ncol(liver), "cells |",
    "median nFeature:", round(median(liver$nFeature_RNA)),
    "| median %MT:", round(median(liver$percent.mt), 1), "%\n")

liver <- subset(
  liver,
  subset = nFeature_RNA > 200  &
           nFeature_RNA < 5000 &
           nCount_RNA   > 500  &
           nCount_RNA   < 25000 &
           percent.mt   < 20
)
cat("  Post-filter:", ncol(liver), "cells\n")

## ============================================================
## STEP 6 — Stratified downsampling to ≤ 20,000 cells
## ============================================================
cat("\nSTEP 6: Stratified downsampling to max 20,000 cells...\n")

MAX_CELLS  <- 20000
total_pre  <- ncol(liver)
cond_table <- table(liver$condition)

cat("  Pre-downsample condition counts:\n")
print(cond_table)

if (total_pre > MAX_CELLS) {
  # Proportional stratified sampling per condition
  # Pre-compute per-condition target counts then sample
  cell_df <- data.frame(
    cell_id   = colnames(liver),
    condition = liver$condition,
    stringsAsFactors = FALSE
  )
  cond_sizes <- cell_df %>%
    count(condition) %>%
    mutate(n_keep = round(MAX_CELLS * n / total_pre))

  sampled <- cell_df %>%
    left_join(cond_sizes[, c("condition","n_keep")], by = "condition") %>%
    group_by(condition) %>%
    group_split() %>%
    lapply(function(df) slice_sample(df, n = unique(df$n_keep))) %>%
    bind_rows() %>%
    pull(cell_id)

  liver <- subset(liver, cells = sampled)
  cat("  Downsampled to", ncol(liver), "cells\n")
} else {
  cat("  Total", total_pre, "cells ≤ 20,000 — no downsampling needed\n")
}

cat("  Post-downsample condition counts:\n")
print(table(liver$condition))

ds_df <- liver@meta.data %>%
  count(condition, sample, name = "n_cells")
write.csv(ds_df, "results/downsampling_summary.csv", row.names = FALSE)

## ============================================================
## STEP 7 — Normalisation and HVG selection
## ============================================================
cat("\nSTEP 7: Normalisation + HVG selection...\n")

liver <- NormalizeData(liver, normalization.method = "LogNormalize",
                       scale.factor = 10000, verbose = FALSE)
liver <- FindVariableFeatures(liver, selection.method = "vst",
                              nfeatures = 2000, verbose = FALSE)

top10 <- head(VariableFeatures(liver), 10)
cat("  Top 10 HVGs:", paste(top10, collapse = ", "), "\n")

p_hvg <- VariableFeaturePlot(liver) +
  LabelPoints(plot = VariableFeaturePlot(liver), points = top10,
              repel = TRUE, xnudge = 0, ynudge = 0) +
  theme_classic() +
  ggtitle("Highly Variable Genes (top 2000)")

ggsave("plots/02_variable_features.png", p_hvg, width = 8, height = 5, dpi = 150)
cat("  Saved plots/02_variable_features.png\n")

## ============================================================
## STEP 8 — Scaling and PCA
## ============================================================
cat("\nSTEP 8: Scaling and PCA...\n")

liver <- ScaleData(liver, features = rownames(liver), verbose = FALSE)
liver <- RunPCA(liver, npcs = 50, verbose = FALSE)

p_elbow <- ElbowPlot(liver, ndims = 50) +
  geom_vline(xintercept = 30, linetype = "dashed", colour = "red") +
  ggtitle("PCA Elbow Plot (red dashed = 30 PCs used)")

ggsave("plots/03_pca_elbow.png", p_elbow, width = 7, height = 4, dpi = 150)
cat("  Saved plots/03_pca_elbow.png\n")

N_PCS <- 30

## ============================================================
## STEP 9 — UMAP and clustering
## ============================================================
cat("\nSTEP 9: UMAP and clustering (", N_PCS, "PCs)...\n")

liver <- FindNeighbors(liver, dims = 1:N_PCS, verbose = FALSE)
liver <- FindClusters(liver,  resolution = 0.5, verbose = FALSE)
liver <- RunUMAP(liver,       dims = 1:N_PCS, verbose = FALSE)

n_clust <- length(levels(liver$seurat_clusters))
cat("  Clusters found:", n_clust, "\n")

p_clust <- DimPlot(liver, reduction = "umap", group.by = "seurat_clusters",
                   label = TRUE, label.size = 4, repel = TRUE) +
  theme_classic() + ggtitle("UMAP — Clusters")

p_cond  <- DimPlot(liver, reduction = "umap", group.by = "condition",
                   cols = COND_COLS, pt.size = 0.3) +
  theme_classic() + ggtitle("UMAP — Condition")

p_frac  <- DimPlot(liver, reduction = "umap", group.by = "fraction",
                   pt.size = 0.3) +
  theme_classic() + ggtitle("UMAP — CD45 Fraction")

ggsave("plots/04_umap_clusters_condition.png",
       p_clust | p_cond | p_frac, width = 18, height = 6, dpi = 150)
cat("  Saved plots/04_umap_clusters_condition.png\n")

## ============================================================
## STEP 10 — Cell type annotation
## ============================================================
cat("\nSTEP 10: Cell type annotation...\n")

# Canonical liver markers
# Hepatocyte:    ALB, APOE, CYP3A4, TF, AFP
# Macrophage:    TREM2, CD68, MARCO, VSIG4, CLEC4F
# Stellate (HSC): ACTA2, COL1A1, COL1A2, PDGFRB, LUM
# Endothelial:   PECAM1, VWF, CLEC4M, LYVE1
# T Cell:        CD3D, CD3E, CD8A, CD4, TRAC
# B Cell:        CD79A, MS4A1, CD19
# NK Cell:       GNLY, NKG7, KLRD1
# Cholangiocyte: EPCAM, KRT19, KRT7

marker_sets <- list(
  Hepatocyte    = c("ALB","APOE","CYP3A4","TF"),
  Macrophage    = c("TREM2","CD68","MARCO","VSIG4"),
  Stellate      = c("ACTA2","COL1A1","COL1A2","PDGFRB"),
  Endothelial   = c("PECAM1","VWF","CLEC4M","LYVE1"),
  "T Cell"      = c("CD3D","CD3E","CD8A","CD4"),
  "B Cell"      = c("CD79A","MS4A1","CD19"),
  "NK Cell"     = c("GNLY","NKG7","KLRD1"),
  Cholangiocyte = c("EPCAM","KRT19","KRT7")
)

markers_plot <- unique(unlist(marker_sets))
markers_plot <- markers_plot[markers_plot %in% rownames(liver)]
cat("  Markers present in data:", length(markers_plot), "/", length(unique(unlist(marker_sets))), "\n")

# Dot plot
p_dot <- DotPlot(liver, features = markers_plot, group.by = "seurat_clusters") +
  RotatedAxis() +
  scale_colour_gradient2(low = "steelblue", mid = "white", high = "red3", midpoint = 0) +
  theme_classic() +
  ggtitle("Marker Expression per Cluster") +
  xlab("") + ylab("Cluster")

ggsave("plots/05_dotplot_markers.png", p_dot,
       width = max(10, length(markers_plot) * 0.65), height = 7, dpi = 150)
cat("  Saved plots/05_dotplot_markers.png\n")

# Feature plots for key markers
feat_markers <- c("ALB","TREM2","CD68","ACTA2","PECAM1","CD3D")
feat_markers <- feat_markers[feat_markers %in% rownames(liver)]

p_feat <- FeaturePlot(liver, features = feat_markers, ncol = 3,
                      order = TRUE, cols = c("lightgrey", "red3")) &
  theme_classic() &
  theme(plot.title = element_text(face = "bold"))

ggsave("plots/06_feature_plots_markers.png", p_feat,
       width = 15, height = ceiling(length(feat_markers) / 3) * 5, dpi = 150)
cat("  Saved plots/06_feature_plots_markers.png\n")

# ── Score each cluster against each cell type marker set ─────
avg_expr <- AverageExpression(
  liver, features = markers_plot,
  group.by = "seurat_clusters",
  assays = "RNA", layer = "data"
)$RNA

score_cluster <- function(cluster_col) {
  sapply(marker_sets, function(genes) {
    g <- intersect(genes, rownames(avg_expr))
    if (length(g) == 0) return(0)
    mean(cluster_col[g])
  })
}

scores_mat <- apply(avg_expr, 2, score_cluster)  # rows = cell types, cols = clusters
celltype_map <- apply(scores_mat, 2, function(col) {
  winner <- names(which.max(col))
  if (max(col) < 0.05) "Other" else winner
})

# AverageExpression prepends "g" to numeric cluster names — strip it so names
# match the actual seurat_clusters levels (e.g. "g0" → "0")
names(celltype_map) <- sub("^g(\\d+)$", "\\1", names(celltype_map))

cat("\n  Cluster → Cell Type assignment:\n")
for (cl in names(celltype_map)) cat("    Cluster", cl, "→", celltype_map[[cl]], "\n")

liver$cell_type <- unname(celltype_map[as.character(liver$seurat_clusters)])

# Colour handling — ensure all present cell types have a colour
present_ct <- unique(liver$cell_type)
ct_cols_used <- CELLTYPE_COLS[names(CELLTYPE_COLS) %in% present_ct]
missing_ct <- setdiff(present_ct, names(ct_cols_used))
if (length(missing_ct) > 0) {
  extra <- setNames(
    colorRampPalette(brewer.pal(8, "Dark2"))(length(missing_ct)),
    missing_ct
  )
  ct_cols_used <- c(ct_cols_used, extra)
}

# Save mapping
write.csv(
  data.frame(cluster = names(celltype_map), cell_type = unname(celltype_map)),
  "results/cluster_celltype_map.csv", row.names = FALSE
)

# UMAP coloured by cell type
p_ct <- DimPlot(liver, reduction = "umap", group.by = "cell_type",
                cols = ct_cols_used, label = TRUE, label.size = 3,
                repel = TRUE, pt.size = 0.3) +
  theme_classic() + ggtitle("UMAP — Cell Type Annotation")

ggsave("plots/07_umap_cell_types.png", p_ct, width = 11, height = 7, dpi = 150)
cat("  Saved plots/07_umap_cell_types.png\n")

# Split UMAP by condition
p_split <- DimPlot(liver, reduction = "umap", group.by = "cell_type",
                   split.by = "condition", cols = ct_cols_used,
                   pt.size = 0.3, label = FALSE) +
  theme_classic() + ggtitle("Cell Types by Condition") +
  theme(legend.position = "bottom")

ggsave("plots/08_umap_split_condition.png", p_split, width = 14, height = 6, dpi = 150)
cat("  Saved plots/08_umap_split_condition.png\n")

## ============================================================
## STEP 11 — Healthy vs Cirrhosis cell type proportion analysis
## ============================================================
cat("\nSTEP 11: Cell type proportion comparison...\n")

prop_df <- liver@meta.data %>%
  group_by(condition, cell_type) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(condition) %>%
  mutate(proportion = n / sum(n) * 100) %>%
  ungroup()

prop_sample_df <- liver@meta.data %>%
  group_by(condition, donor, cell_type) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(condition, donor) %>%
  mutate(proportion = n / sum(n) * 100) %>%
  ungroup()

write.csv(prop_df,        "results/celltype_proportions.csv",        row.names = FALSE)
write.csv(prop_sample_df, "results/celltype_proportions_per_donor.csv", row.names = FALSE)

# Stacked bar
ct_order <- prop_df %>%
  filter(condition == "Healthy") %>%
  arrange(desc(proportion)) %>%
  pull(cell_type)
ct_order <- c(ct_order, setdiff(unique(prop_df$cell_type), ct_order))
prop_df$cell_type <- factor(prop_df$cell_type, levels = ct_order)

p_bar <- ggplot(prop_df, aes(x = condition, y = proportion, fill = cell_type)) +
  geom_bar(stat = "identity", width = 0.6, colour = "white", linewidth = 0.3) +
  scale_fill_manual(values = ct_cols_used, name = "Cell Type") +
  scale_y_continuous(expand = c(0, 0), labels = function(x) paste0(x, "%")) +
  labs(title = "Cell Type Composition: Healthy vs Cirrhosis",
       x = NULL, y = "Proportion of Cells") +
  theme_classic(base_size = 13) +
  theme(axis.text.x = element_text(size = 12, face = "bold"),
        legend.position = "right")

ggsave("plots/09_celltype_proportion_stacked.png", p_bar, width = 8, height = 6, dpi = 150)
cat("  Saved plots/09_celltype_proportion_stacked.png\n")

# Grouped bar
p_grp <- ggplot(prop_df, aes(x = cell_type, y = proportion, fill = condition)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  scale_fill_manual(values = COND_COLS) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  labs(title = "Cell Type Proportions: Healthy vs Cirrhosis",
       x = NULL, y = "Proportion (%)") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1),
        legend.position = "top")

ggsave("plots/10_celltype_proportion_grouped.png", p_grp, width = 10, height = 6, dpi = 150)
cat("  Saved plots/10_celltype_proportion_grouped.png\n")

# TREM2 / CD68 split by condition — scar-associated macrophages in cirrhosis
if (all(c("TREM2","CD68") %in% rownames(liver))) {
  p_mac <- FeaturePlot(liver, features = c("TREM2","CD68"),
                       split.by = "condition", ncol = 2,
                       cols = c("lightgrey","darkred"), order = TRUE) &
    theme_classic() & theme(plot.title = element_text(face = "bold"))
  ggsave("plots/11_trem2_cd68_by_condition.png", p_mac,
         width = 12, height = 6, dpi = 150)
  cat("  Saved plots/11_trem2_cd68_by_condition.png\n")
}

## ============================================================
## STEP 12 — Cluster marker genes
## ============================================================
cat("\nSTEP 12: Finding cluster marker genes...\n")

Idents(liver) <- "cell_type"
markers_all <- FindAllMarkers(
  liver,
  only.pos        = TRUE,
  min.pct         = 0.25,
  logfc.threshold = 0.5,
  test.use        = "wilcox",
  verbose         = FALSE
)

top5 <- markers_all %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 5)

write.csv(markers_all, "results/celltype_marker_genes.csv",      row.names = FALSE)
write.csv(top5,        "results/top5_markers_per_celltype.csv",  row.names = FALSE)
cat("  Marker gene tables saved\n")

# Heatmap: top 3 markers per cell type
top3_genes <- markers_all %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 3) %>%
  pull(gene) %>%
  unique() %>%
  intersect(rownames(liver))

if (length(top3_genes) > 0) {
  p_heat <- DoHeatmap(
    liver, features = top3_genes, group.by = "cell_type",
    group.colors = ct_cols_used, label = TRUE, size = 3
  ) +
    scale_fill_gradient2(low = "steelblue", mid = "white", high = "red3",
                         name = "Scaled\nExpression") +
    ggtitle("Top 3 Marker Genes per Cell Type")

  ggsave("plots/12_heatmap_top_markers.png", p_heat,
         width = max(14, length(top3_genes) * 0.4 + 4), height = 9, dpi = 150)
  cat("  Saved plots/12_heatmap_top_markers.png\n")
}

## ============================================================
## STEP 13 — Save Seurat object and summary
## ============================================================
cat("\nSTEP 13: Saving results...\n")

summary_tbl <- liver@meta.data %>%
  group_by(condition, cell_type) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  pivot_wider(names_from = condition, values_from = n_cells, values_fill = 0)
write.csv(summary_tbl, "results/cell_counts_summary.csv", row.names = FALSE)

saveRDS(liver, "data/liver_seurat_annotated.rds")
cat("  Saved data/liver_seurat_annotated.rds\n")

## ============================================================
## FINAL SUMMARY
## ============================================================
cat("\n========================================================\n")
cat("  ANALYSIS COMPLETE — GSE136103 Human Liver scRNA-seq\n")
cat("========================================================\n")
cat("  Cells:      ", ncol(liver), "\n")
cat("  Conditions: ", paste(unique(liver$condition), collapse = " vs "), "\n")
cat("  Samples:    ", length(unique(liver$sample)), "(liver fractions)\n")
cat("  Clusters:   ", length(unique(liver$seurat_clusters)), "\n")
cat("  Cell types: ", length(unique(liver$cell_type)), "\n")
cat("\n  Cell type proportions:\n")
prop_wide <- prop_df %>%
  select(cell_type, condition, proportion) %>%
  pivot_wider(names_from = condition, values_from = proportion, values_fill = 0) %>%
  mutate(across(where(is.numeric), ~round(.x, 1))) %>%
  arrange(desc(Healthy))
print(as.data.frame(prop_wide), row.names = FALSE)
cat("========================================================\n")
cat("  Plots  → plots/\n")
cat("  Tables → results/\n")
cat("  Object → data/liver_seurat_annotated.rds\n")
cat("========================================================\n")
