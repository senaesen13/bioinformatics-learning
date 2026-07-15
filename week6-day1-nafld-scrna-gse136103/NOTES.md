# Week 6 Day 1 — NAFLD scRNA-seq: GSE136103 Setup and Seurat Pipeline

Human liver single-cell RNA-seq data from Ramachandran et al. 2019 (Nature).
This folder covers dataset download, NAFLD-specific subsetting, and the full
Seurat pipeline through clustering and marker gene identification. Cell-type
annotation is in a subsequent session.

---

## Dataset Overview

**GEO Accession:** GSE136103  
**Paper:** Ramachandran et al. 2019, *Nature* — "Resolving the fibrotic niche of
human liver cirrhosis at single-cell level"  
**Organism:** Homo sapiens only (all 24 GEO samples; no mouse data in this accession)

### CD45 sorting design

Cells were sorted before sequencing into two fractions:
- **CD45+** (Leucocytes): immune / inflammatory cells
- **CD45−** (Other NPC): hepatic stellate cells, endothelial cells, cholangiocytes

This means each biological donor contributes 1–3 sequencing libraries (one per sort
fraction), so library count > donor count.

---

## NAFLD-Specific Subset (scripts 02 + 03)

From 24 human libraries (14 donors), the analysis retains only NAFLD-relevant liver
tissue:

| Group | Donors | Libraries | Etiology |
|---|---|---|---|
| Healthy liver | 5 (Healthy1–5) | 11 | None |
| NAFLD cirrhosis | 2 (Cirrhotic1, Cirrhotic4) | 4 | NAFLD |
| **Total** | **7** | **15** | |

**Excluded (9 libraries):**
- Cirrhotic2 (GSM4041164–165) — Alcohol cirrhosis
- Cirrhotic3 (GSM4041166–167) — Alcohol cirrhosis
- Cirrhotic5 (GSM4041169) — PBC
- Blood1–4 (GSM4041170–173) — PBMC, not liver tissue

### Cirrhotic4 CD45− gap

Cirrhotic4 (GSM4041168) has only one library (CD45+ fraction). No CD45− library was
submitted to GEO for this donor. This means Cirrhotic4 contributes immune/leucocyte
cells to the merged object but **no hepatic stellate cells, endothelial cells, or
cholangiocytes**. Cell-type composition will therefore be asymmetric between Cirrhotic1
(3 libraries: CD45+, CD45−A, CD45−B) and Cirrhotic4 (1 library: CD45+ only). This
caveat applies to any cluster-level healthy vs. cirrhosis comparison that involves
non-immune liver cell types.

---

## Seurat Pipeline (script 03)

### Step 1 — Load libraries

`ReadMtx()` used instead of `Read10X()` because GEO filenames are non-standard
(e.g. `GSM4041150_healthy1_cd45+_matrix.mtx.gz`). `feature.column = 2` selects
gene symbol over Ensembl ID. Feature names with underscores are auto-converted to
dashes by Seurat (cosmetic warning only).

**Seurat v5 note:** After `merge()`, count data is stored in per-sample layers.
`JoinLayers()` is called immediately after merge so that `FindAllMarkers` and other
DE functions can see a single unified count matrix.

### Step 2 — Merge

| | |
|---|---|
| Libraries merged | 15 |
| Total cells (raw) | 46,331 |
| Total features | 22,517 |

`JoinLayers()` called immediately after merge (Seurat v5 requirement for DE testing).

Per-library cell counts:

| GSM | Donor | Fraction | Cells |
|---|---|---|---|
| GSM4041150 | Healthy1 | CD45pos | 1,860 |
| GSM4041151 | Healthy1 | CD45neg | 969 |
| GSM4041152 | Healthy1 | CD45neg | 436 |
| GSM4041153 | Healthy2 | CD45pos | 6,571 |
| GSM4041154 | Healthy2 | CD45neg | 4,684 |
| GSM4041155 | Healthy3 | CD45pos | 3,295 |
| GSM4041156 | Healthy3 | CD45neg | 3,207 |
| GSM4041157 | Healthy3 | CD45neg | 1,131 |
| GSM4041158 | Healthy4 | CD45pos | 4,932 |
| GSM4041159 | Healthy4 | CD45neg | 3,306 |
| GSM4041160 | Healthy5 | CD45pos | 5,121 |
| GSM4041161 | Cirrhotic1 | CD45pos | 2,344 |
| GSM4041162 | Cirrhotic1 | CD45neg | 1,652 |
| GSM4041163 | Cirrhotic1 | CD45neg | 2,004 |
| GSM4041168 | Cirrhotic4 | CD45pos | 4,819 |

### Step 3 — QC

Thresholds (matching Week 3 PBMC 3k analysis):
- `nFeature_RNA`: 200 – 2,500
- `percent.mt`: < 5 %

| | |
|---|---|
| Cells before QC | 46,331 |
| Cells after QC | 35,050 (75.7% retained) |
| Cells removed | 11,281 (24.3%) |

Pre-filter median percent.mt was 3.2%; max was 95.5%, indicating a small fraction of
dying/damaged cells that were caught by the 5% threshold.

**Threshold note:** The 2,500 nFeature cap is conservative for liver NPC cells (hepatic
stellate cells and endothelial cells can have high feature counts). If cluster
inspection suggests doublet-like clusters at the top of the feature distribution,
the cap can be raised to 4,000 in a subsequent run.

### Step 4 — Normalisation + variable features

- Method: LogNormalize, scale factor 10,000
- Variable features: top 2,000 by VST

### Step 5 — Scale + PCA

- Scaled features: variable features only (2,000) for speed
- PCA: 50 PCs computed; N_PCS = 20 used for clustering (see elbow plot)

### Step 6 — Clustering + UMAP

| Parameter | Value |
|---|---|
| Dims | 1:20 |
| Resolution | 0.5 |
| Clusters found | **21** |
| UMAP method | UWOT (R-native cosine metric; Seurat v5 default) |

Cells per cluster (sorted ascending):

| Cluster | Cells | Cluster | Cells |
|---|---|---|---|
| 20 | 77 | 5 | 2,658 |
| 19 | 189 | 4 | 3,155 |
| 18 | 238 | 3 | 3,185 |
| 17 | 345 | 2 | 3,706 |
| 16 | 362 | 1 | 5,166 |
| 15 | 373 | 0 | 5,561 |
| 14 | 463 | | |
| 13 | 736 | **Total** | **35,050** |
| 12 | 932 | | |
| 11 | 1,027 | | |
| 10 | 1,052 | | |
| 9 | 1,061 | | |
| 8 | 1,132 | | |
| 7 | 1,162 | | |
| 6 | 2,470 | | |

### Step 7 — Marker genes

`FindAllMarkers` with Wilcoxon rank-sum test:
- `only.pos = TRUE`
- `min.pct = 0.25`
- `logfc.threshold = 0.25`
- Total significant marker rows: **12,313**

Full marker table: `results/all_markers.csv`  
Top 5 per cluster: `results/top5_markers_per_cluster.csv`

Top 5 markers per cluster (top gene by avg_log2FC):

| Cluster | Top markers | Likely identity (not yet annotated) |
|---|---|---|
| 0 | CD40LG, IL7R, LTB, AC092580.4, GPR171 | CD4+ T cells |
| 1 | TOX2, XCL1, XCL2, KLRC1, IL2RB | NK / NKT cells |
| 2 | CD8B, CD3G, GZMH, CD8A, CLEC2D | CD8+ T cells |
| 3 | CD40LG, IL7R, AQP3, TRAC, CXCR4 | CD4+ T cells (memory) |
| 4 | CD8B, CD8A, CRTAM, LAG3, TRGC2 | CD8+ T cells (cytotoxic/exhausted) |
| 5 | FGFBP2, GNLY, GZMB, CX3CR1, SPON2 | NK cells (cytotoxic) |
| 6 | S100A12, S100A9, S100A8, FCN1, VCAN | Inflammatory monocytes |
| 7 | AQP7, AIF1L, GPIHBP1, PODXL, HSPA12B | Endothelial cells |
| 8 | VPREB3, CD19, LINC00926, FCRLA, CD79A | B cells |
| 9 | CLEC1B, OIT3, FCN2, CLEC4G, FCN3 | Liver sinusoidal EC (LSEC) |
| 10 | FCER1A, IL1R2, CD1C, CLEC10A, PKIB | cDC2 dendritic cells |
| 11 | CD5L, SDC3, C1QC, CD163, C1QB | Kupffer cells |
| 12 | NTF3, TCF21, RGS5, DCN, FOXS1 | Hepatic stellate cells / pericytes |
| 13 | UGT2B15, UGT2A3, DCDC2, FXYD2 | Hepatocytes |
| 14 | GZMH, FGFBP2, TRGC2, CX3CR1 | NK cells (second subset) |
| 15 | IGHGP, IGLL5, IGHA2, AMPD1, IGHG2 | Plasma cells |
| 16 | MAL, LEF1, CCR7 | Naive / central memory T cells |
| 17 | SPC25, ASPM, CENPA, RRM2, TYMS | Proliferating cells |
| 18 | IL2RB, CD160, NCR1, XCL1 | Liver-resident NK cells |
| 19 | SCT, PTCRA, LRRC26 | Rare epithelial / cholangiocyte? |
| 20 | TPSD1, TPSAB1, TPSB2, CPA3 | Mast cells |

**Note:** "Likely identity" column is preliminary read of top markers only — formal annotation step is separate.

---

## Next Steps

1. Review elbow plot (`plots/pca_elbow.png`) — confirm N_PCS = 20 is appropriate
2. Review UMAP panel (`plots/umap_panel.png`) — check for donor batch effects
3. Review top markers per cluster and annotate cell types (hepatocytes,
   Kupffer cells, HSC, endothelial, NK/T, B cells, plasma cells, etc.)
4. If strong donor batch effects visible: re-run with Harmony integration
5. NAFLD differential analysis: compare healthy vs NAFLD cirrhosis per cell type

---

## Output Files

| File | Description |
|---|---|
| `scripts/01_download_setup.R` | GEO metadata download + human organism filter |
| `scripts/02_nafld_subset.R` | NAFLD-specific sample subsetting |
| `scripts/03_seurat_pipeline.R` | Full Seurat pipeline (QC → clustering → markers) |
| `results/human_sample_metadata.csv` | All 24 human GSM metadata |
| `results/nafld_subset_metadata.csv` | 15-library NAFLD subset metadata |
| `results/nafld_seurat_clustered.rds` | Clustered Seurat object (not in git) |
| `results/all_markers.csv` | All FindAllMarkers output |
| `results/top5_markers_per_cluster.csv` | Top 5 markers per cluster by log2FC |
| `plots/qc_violin_before.png` | QC violin plots before filtering |
| `plots/qc_violin_after.png` | QC violin plots after filtering |
| `plots/variable_features.png` | Top variable features plot |
| `plots/pca_elbow.png` | Elbow plot (PC selection) |
| `plots/pca_donor.png` | PCA coloured by donor |
| `plots/umap_clusters.png` | UMAP: cluster labels |
| `plots/umap_disease.png` | UMAP: healthy vs NAFLD cirrhosis |
| `plots/umap_donor.png` | UMAP: coloured by donor |
| `plots/umap_fraction.png` | UMAP: CD45+ vs CD45− fraction |
| `plots/umap_panel.png` | 2×2 UMAP panel |
| `plots/marker_heatmap.png` | Heatmap of top 5 markers per cluster |
| `data/GSM404XXXX/` | Raw 10x MEX files — not tracked in git (380 MB) |
