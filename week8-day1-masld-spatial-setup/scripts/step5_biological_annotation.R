#!/usr/bin/env Rscript
# =============================================================================
# Week 8 Day 1 — Step 5: Biological Annotation of 24 Spatial Clusters
# FindAllMarkers → biological label assignment → UMAP with labels
# =============================================================================

suppressMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

PROJ_ROOT <- "/Users/senaesen/Desktop/bioinfo-learning/week8-day1-masld-spatial-setup"
RESULTS   <- file.path(PROJ_ROOT, "results")
DOM_DIR   <- file.path(RESULTS, "spatial_domains")

cat("Loading combined Seurat object...\n")
combined <- readRDS(file.path(RESULTS, "masld_spatial_combined.rds"))
cat(sprintf("  Object: %d spots × %d genes, %d clusters\n",
    ncol(combined), nrow(combined),
    length(unique(combined$seurat_clusters))))

# =============================================================================
# STEP 5.1 — FindAllMarkers
# =============================================================================
cat("\n=== Step 5.1: FindAllMarkers ===\n")

DefaultAssay(combined) <- "SCT"

# Seurat v5: JoinLayers on the Spatial assay (Assay5) before PrepSCTFindMarkers.
# SCTAssay does not support JoinLayers — only the raw counts assay needs joining.
cat("  JoinLayers on Spatial assay (Seurat v5 requirement after merge)...\n")
combined <- JoinLayers(combined, assay = "Spatial")

# PrepSCTFindMarkers: recalculates SCT residuals for cross-sample DE testing
cat("  PrepSCTFindMarkers...\n")
combined <- PrepSCTFindMarkers(combined, verbose = FALSE)

cat("  FindAllMarkers (Wilcoxon, only.pos=TRUE, min.pct=0.1, logfc.threshold=0.25)...\n")
markers_all <- FindAllMarkers(
  combined,
  assay           = "SCT",
  only.pos        = TRUE,
  min.pct         = 0.1,
  logfc.threshold = 0.25,
  test.use        = "wilcox",
  verbose         = FALSE
)

write.csv(markers_all,
  file.path(DOM_DIR, "cluster_markers_all.csv"),
  row.names = FALSE)
cat(sprintf("  Saved: %d marker genes total\n", nrow(markers_all)))

# Top 10 significant markers per cluster (by avg_log2FC)
top10 <- markers_all %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 10) %>%
  arrange(cluster, desc(avg_log2FC))

write.csv(top10,
  file.path(DOM_DIR, "cluster_markers_top10.csv"),
  row.names = FALSE)
cat(sprintf("  Top 10 per cluster saved (%d rows)\n", nrow(top10)))

# =============================================================================
# STEP 5.2 — Biological annotation
# =============================================================================
cat("\n=== Step 5.2: Biological annotation ===\n")

# Liver-relevant marker gene sets — chosen to distinguish spatial hepatocyte
# zones and non-parenchymal cell niches in MASLD tissue
MARKERS <- list(
  periportal    = c("ALB","ASS1","CPS1","PCK1","GLS2","HAL","APOB","APOC3",
                    "FABP1","ARG1","OTC","TAT","HPD","ACOX1","SLC1A2","UGP2",
                    "PYGB","PGM1","GBE1","AGL","PFKL"),
  pericentral   = c("CYP2E1","GLUL","CYP1A2","CYP3A4","CYP7A1","ADH4",
                    "ALDH3A2","UGT1A4","ACSS2","VNN1","CYP2C8","CYP2C9",
                    "CYP1A1","ALDH1L1","SULT1A1"),
  midzonal      = c("APOE","PPARA","HMGCR","LDLR","INSIG1","FDPS","SQLE",
                    "MVK","DHCR7","SC5D","EBP","TM7SF2","NSDHL","HSD17B7",
                    "MSMO1","LSS","IDI1"),
  stellate      = c("COL1A1","COL3A1","COL1A2","ACTA2","TAGLN","VIM","THY1",
                    "PDGFRB","FN1","POSTN","LOXL2","MMP2","MMP14","TIMP1",
                    "TIMP2","PDGFRA","DCN","LUM","SPARC","MFAP4"),
  kupffer_macro = c("TREM2","SPP1","CD68","MARCO","C1QB","C1QA","AIF1","VSIG4",
                    "CD163","MNDA","MRC1","FCGR3A","CCL2","HMOX1","APOC1",
                    "LILRB5","GPNMB","FOLR2","CLEC7A","MS4A7"),
  lsec          = c("LYVE1","STAB2","RAMP3","FCN3","FCGR2B","STAB1","CRHBP",
                    "AQP1","PTPRB","PECAM1","ENG","THBD","CLEC4M","OIT3",
                    "CALCRL","EMCN","F8","VWF","KDR","ESAM"),
  biliary       = c("KRT7","KRT19","SOX9","EPCAM","CFTR","ANXA4","TFF1",
                    "TFF2","MUC1","CLDN4","CLDN7","CLDN3","FXYD2","DEFB1"),
  lipogenic     = c("FASN","ACACA","ACLY","SCD","GPAM","SREBF1","MLXIPL",
                    "ACSL5","ELOVL6","PLIN2","CIDEC","DGAT2","MOGAT1",
                    "G6PC","GCK","PKLR"),
  inflammatory  = c("IL6","TNF","CXCL10","CCL20","CXCL1","ICAM1","NFKB1",
                    "IL1B","LCN2","SAA1","SAA2","CRP","SERPINA3","ORM1",
                    "ORM2","HP","APCS"),
  senescent     = c("CDKN2A","CDKN1A","TP53","SERPINE1","IGFBP3","IGFBP7",
                    "THBS1","CCN1","GDF15","CDKN2B","CCND1","GLB1","LMNB1",
                    "MKI67","CDK4")
)

# Human-readable label for each category
LABEL_MAP <- c(
  periportal    = "Periportal hepatocyte zone",
  pericentral   = "Pericentral hepatocyte zone",
  midzonal      = "Midzonal hepatocyte zone",
  stellate      = "Fibrotic/stellate-activated region",
  kupffer_macro = "Inflammatory/Kupffer-macrophage region",
  lsec          = "Sinusoidal endothelial (LSEC) region",
  biliary       = "Biliary/portal tract region",
  lipogenic     = "Lipogenic hepatocytes",
  inflammatory  = "Pro-inflammatory/acute-phase region",
  senescent     = "Senescence-associated region"
)

# Genes previously validated in this project's bulk RNA-seq and scRNA-seq
KEY_GENES <- c("TREM2", "SPP1", "COL1A1", "FASN")

# Score cluster top markers against each biological category
score_markers <- function(gene_vec) {
  g_up <- toupper(gene_vec)
  sapply(MARKERS, function(gset) sum(g_up %in% toupper(gset)))
}

all_clusters <- sort(as.integer(as.character(unique(combined$seurat_clusters))))

annotation_rows <- lapply(all_clusters, function(cl) {
  cl_str <- as.character(cl)

  sig_markers <- markers_all %>%
    filter(cluster == cl_str, p_val_adj < 0.05) %>%
    arrange(desc(avg_log2FC))

  cl_genes_20 <- head(sig_markers$gene, 20)
  cl_genes_5  <- head(sig_markers$gene, 5)
  n_sig       <- nrow(sig_markers)

  n_spots  <- sum(combined$seurat_clusters == cl_str)
  pct_spots <- round(100 * n_spots / ncol(combined), 1)
  top5_str <- paste(cl_genes_5, collapse = ", ")

  # Key gene hits in ALL significant markers (not just top 20)
  key_hits <- intersect(toupper(sig_markers$gene), toupper(KEY_GENES))
  key_str  <- if (length(key_hits) > 0) paste(key_hits, collapse = ", ") else ""

  if (n_sig < 3) {
    label <- "Indeterminate (insufficient significant markers)"
  } else {
    scores    <- score_markers(cl_genes_20)
    top_cat   <- names(which.max(scores))
    top_score <- max(scores)

    if (top_score < 2) {
      label <- "Hepatocyte parenchyma (zone indeterminate)"
    } else {
      label <- LABEL_MAP[top_cat]

      # Check for secondary enrichment (within 1 of top score) and append
      scores_sorted <- sort(scores, decreasing = TRUE)
      second_cat   <- names(scores_sorted)[2]
      second_score <- scores_sorted[2]
      if (second_score >= 2 && second_score >= top_score - 1 && second_cat != top_cat) {
        second_label <- LABEL_MAP[second_cat]
        label <- paste0(label, " / ", second_label)
      }
    }
  }

  data.frame(
    cluster_id        = cl,
    biological_label  = label,
    top5_markers      = top5_str,
    key_gene_hits     = key_str,
    n_sig_markers     = n_sig,
    n_spots           = n_spots,
    pct_spots         = pct_spots,
    stringsAsFactors  = FALSE
  )
})

annotation_df <- do.call(rbind, annotation_rows)

write.csv(annotation_df,
  file.path(DOM_DIR, "cluster_biological_annotation.csv"),
  row.names = FALSE)

cat("\nBiological annotation summary:\n")
cat(sprintf("%-5s %-52s %-35s %s\n", "Cl", "Label", "Top 5 markers", "Key genes"))
cat(strrep("-", 110), "\n")
for (i in seq_len(nrow(annotation_df))) {
  r <- annotation_df[i, ]
  label_str <- substr(r$biological_label, 1, 50)
  markers_str <- substr(r$top5_markers, 1, 33)
  cat(sprintf("C%-4d %-52s %-35s %s\n",
      r$cluster_id, label_str, markers_str,
      if (r$key_gene_hits == "") "(none)" else r$key_gene_hits))
}

# Report cross-modality confirmation
cat("\n--- Cross-modality confirmation (TREM2/SPP1/COL1A1/FASN) ---\n")
hits_df <- annotation_df %>% filter(key_gene_hits != "")
if (nrow(hits_df) > 0) {
  for (i in seq_len(nrow(hits_df))) {
    cat(sprintf("  C%d (%s): %s\n",
        hits_df$cluster_id[i],
        hits_df$biological_label[i],
        hits_df$key_gene_hits[i]))
  }
} else {
  cat("  None of TREM2/SPP1/COL1A1/FASN appeared as cluster markers (may reflect probe panel gaps)\n")
}

# =============================================================================
# STEP 5.3 — Regenerate UMAP with biological labels
# =============================================================================
cat("\n=== Step 5.3: UMAP with biological labels ===\n")

# Add labels to Seurat metadata
cl_to_label <- setNames(
  paste0("C", annotation_df$cluster_id, ": ", annotation_df$biological_label),
  as.character(annotation_df$cluster_id)
)
cl_to_short <- setNames(
  annotation_df$biological_label,
  as.character(annotation_df$cluster_id)
)

combined$bio_label      <- cl_to_label[as.character(combined$seurat_clusters)]
combined$bio_label_short <- cl_to_short[as.character(combined$seurat_clusters)]

# Panel 1: numbered clusters
p_clusters <- DimPlot(combined,
  group.by   = "seurat_clusters",
  label      = TRUE, label.size = 3.5,
  raster     = FALSE, repel      = TRUE) +
  ggtitle("24 Spatial Clusters (numbered)") +
  theme(legend.position = "none")

# Panel 2: biological labels
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
  p_clusters | p_bio,
  width = 22, height = 7, dpi = 150
))
cat("  Saved: umap_biological_labels.png\n")

suppressWarnings(ggsave(
  file.path(DOM_DIR, "umap_bio_labels_only.png"),
  p_bio,
  width = 14, height = 7, dpi = 150
))
cat("  Saved: umap_bio_labels_only.png\n")

# =============================================================================
# STEP 5.4 — Summary table (cluster_annotation_summary.csv)
# =============================================================================
summary_df <- annotation_df %>%
  select(cluster_id, biological_label, top5_markers, key_gene_hits, n_spots, pct_spots) %>%
  rename(
    `Cluster ID`       = cluster_id,
    `Biological label` = biological_label,
    `Top 5 markers`    = top5_markers,
    `Key gene hits`    = key_gene_hits,
    `Spot count`       = n_spots,
    `% of total`       = pct_spots
  )

write.csv(summary_df,
  file.path(DOM_DIR, "cluster_annotation_summary.csv"),
  row.names = FALSE)
cat("  Saved: cluster_annotation_summary.csv\n")

# =============================================================================
# STEP 5.5 — Save updated object
# =============================================================================
saveRDS(combined, file.path(RESULTS, "masld_spatial_combined.rds"))
cat("  Updated Seurat object (bio_label added to metadata) saved.\n")

cat("\n=== Step 5 complete ===\n")
cat("Outputs in results/spatial_domains/:\n")
cat("  cluster_markers_all.csv           — all significant markers per cluster\n")
cat("  cluster_markers_top10.csv         — top 10 markers per cluster\n")
cat("  cluster_biological_annotation.csv — label + top5 + key gene hits\n")
cat("  cluster_annotation_summary.csv    — formatted summary table\n")
cat("  umap_biological_labels.png        — dual UMAP (numbered | annotated)\n")
cat("  umap_bio_labels_only.png          — biological labels only\n")
cat("\n")
sessionInfo()
