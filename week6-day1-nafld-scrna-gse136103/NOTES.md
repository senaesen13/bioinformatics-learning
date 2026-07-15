# Week 6 Day 1 — NAFLD scRNA-seq: GSE136103 Setup, Seurat Pipeline, and Composition/DE Analysis

Human liver single-cell RNA-seq data from Ramachandran et al. 2019 (Nature).
This folder covers dataset download, NAFLD-specific subsetting, the full Seurat
pipeline through clustering, cell-type annotation, cell composition comparison,
and within-cell-type differential expression (healthy vs NAFLD cirrhosis).

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

Formal annotation map used in script 04 (cluster → cell type):

| Cluster | Cell type label | Key markers used |
|---|---|---|
| 0, 3 | CD4+ T cells | CD40LG, IL7R, LTB, TRAC |
| 2 | CD8+ T cells | CD8A, CD8B, CD3G |
| 4 | CD8+ T cells (exhausted) | CD8A, LAG3, CRTAM |
| 16 | Naive T cells | MAL, LEF1, CCR7 |
| 1 | NK/NKT cells | TOX2, XCL1, XCL2, KLRC1 |
| 5 | NK cells (cytotoxic) | GNLY, GZMB, FGFBP2, CX3CR1 |
| 14 | gdT/NK cells | GZMH, TRGC2, CX3CR1 |
| 18 | NK cells (liver-resident) | IL2RB, CD160, NCR1 |
| 6 | Monocytes | S100A8/9/12, FCN1, VCAN |
| 10 | Dendritic cells | CD1C, FCER1A, CLEC10A |
| 11 | Kupffer cells | CD5L, C1QA/B/C, CD163, MARCO |
| 8 | B cells | CD79A, CD19, LINC00926 |
| 15 | Plasma cells | IGHGP, IGLL5, IGHA2 |
| 7 | Endothelial cells | GPIHBP1, PODXL, AQP7 |
| 9 | LSEC | CLEC4G, FCN2/3, OIT3, CLEC1B |
| 12 | Hepatic stellate cells | DCN, TCF21, RGS5, NTF3 |
| 13 | Hepatocytes | UGT2B15, UGT2A3, FXYD2 |
| 17 | Proliferating cells | CENPA, RRM2, TYMS, ASPM |
| 19 | Cholangiocytes | SCT, PTCRA, LRRC26 |
| 20 | Mast cells | TPSAB1, TPSB2, TPSD1, CPA3 |

---

## Cell Type Composition: Healthy vs NAFLD Cirrhosis (script 04)

Script: `scripts/04_composition_and_de.R`  
Results: `results/cell_composition_major.csv`, `plots/composition_healthy_vs_nafld.png`

Mean cell type percentages per donor, averaged within each disease group:

| Major cell type | Healthy (mean %) | NAFLD cirrhosis (mean %) | Difference | Fold |
|---|---|---|---|---|
| T cells | 43.83 | 46.41 | +2.6 | 1.06 |
| NK/NKT cells | **24.36** | **16.35** | −8.0 | 0.67 |
| B/Plasma cells | 3.66 | 7.63 | +4.0 | 2.08 |
| Monocytes | 7.61 | 10.70 | +3.1 | 1.40 |
| Endothelial/LSEC | 8.51 | 8.52 | ~0 | 1.00 |
| Kupffer cells | 3.37 | 2.65 | −0.7 | 0.78 |
| Hepatic stellate cells | 2.74 | 1.46 | −1.3 | 0.53 |
| Dendritic cells | 3.12 | 1.67 | −1.5 | 0.54 |
| Hepatocytes | 1.26 | 2.32 | +1.1 | 1.84 |
| Proliferating cells | 0.80 | 0.82 | ~0 | 1.02 |
| Cholangiocytes | 0.53 | 0.99 | +0.5 | 1.86 |
| Mast cells | 0.19 | 0.47 | +0.3 | 2.51 |

### Composition caveats

**Cirrhotic4 CD45− gap:** Cirrhotic4 (GSM4041168) has only a CD45+ library. This
means its contribution to non-immune cell types (hepatocytes, HSC, LSEC, endothelial
cells, cholangiocytes) is near zero. This artificially lowers the proportions of
stromal/parenchymal cells in the NAFLD group and inflates immune cell percentages
in Cirrhotic4. Treat the hepatocyte, HSC, and endothelial composition numbers as
lower bounds for NAFLD rather than accurate estimates.

**Two NAFLD donors:** Composition means and SDs are computed across only 2 donors
(Cirrhotic1 and Cirrhotic4). The ±SD on the bar plot is not statistically
meaningful at n=2. The direction of differences (especially NK/NKT decrease and
monocyte/B cell increase) is consistent with literature on liver cirrhosis, but
cannot be treated as statistically validated.

---

## Within-Cell-Type DE: Healthy vs NAFLD Cirrhosis (script 04)

Full results: `results/within_celltype_de_sig.csv` (8,913 genes, padj < 0.05)  
Top 10 up/down per type: `results/top10_de_per_celltype.csv`

18 of 20 cell types had sufficient cells for DE (≥20 per group). Liver-resident NK
cells and Naive T cells were too rare in NAFLD to test.

### Headline findings

**TREM2 upregulated in NAFLD Monocytes** (avg log2FC = 3.51, padj = 1.16e-26,
pct.1 = 0.116 in NAFLD vs 0.012 in healthy). TREM2+ monocyte-derived macrophages
are the scar-associated macrophage (SAM) population described in the Ramachandran
2019 paper. Their detection here in the Monocyte cluster (rather than Kupffer cell
cluster) is consistent with recruited monocytes differentiating into SAMs in
fibrosis. TREM2 was not a significant cluster-level marker (too few cells per
cluster exceeded the 25% threshold), but is clearly disease-enriched.

**CXCL13 massively upregulated in CD8+ T cells (exhausted)** (log2FC = 12.5 in
NAFLD, padj = 6.5e-83). CXCL13-expressing exhausted CD8+ T cells are a well-described
feature of advanced liver fibrosis and hepatocellular carcinoma risk. Co-expression
of HAVCR2 (TIM-3, log2FC = 4.0) in the same cluster confirms this is a genuine
T cell exhaustion signature.

**HSC fibrosis markers upregulated in NAFLD:**
- ADAMTS2 (procollagen N-proteinase, ECM remodeling) log2FC = 4.07
- RARRES1 (activated stellate cell marker) log2FC = 5.52
- MFAP2 (microfibril-associated protein, fibrosis ECM) log2FC = 4.61

**Hepatocyte stress in NAFLD:**
- NNMT (nicotinamide N-methyltransferase) log2FC = 3.66 — NNMT is upregulated in
  NASH hepatocytes and promotes lipid accumulation via one-carbon metabolism.
- XIST log2FC = 7.05 in hepatocytes — see sex-composition caveat below.

**Kupffer cell C3 upregulation** (log2FC = 4.01): complement component C3 is
upregulated in NAFLD Kupffer cells, consistent with complement activation in
NAFLD/NASH.

**Endothelial POSTN upregulation** (log2FC = 4.18): Periostin (POSTN) is a
fibrosis-associated extracellular matrix protein upregulated in portal fibroblasts
and activated endothelium.

### XIST caveat

XIST (X-inactive specific transcript, the X-chromosome inactivation gene) appears
upregulated in NAFLD across many cell types (hepatocytes, HSC, endothelial, HSC).
XIST is expressed only in female cells (for X-inactivation). Its apparent upregulation
in NAFLD likely reflects a difference in **sex composition between the donors**: if
the NAFLD donors (Cirrhotic1, Cirrhotic4) contain more female donors than the healthy
group, XIST will appear as a NAFLD-upregulated gene spuriously. The GEO metadata does
not record sex for all samples. Treat any XIST-driven signal as a potential sex-
composition confound, not a genuine NAFLD biology signal.

---

## Next Steps

1. Check donor sex metadata (available in pData characteristics columns) to confirm
   XIST/sex-composition confound
2. Re-run DE excluding sex-linked genes (XIST, RPS4Y1, DDX3Y, EIF1AY) or controlling
   for sex as a covariate
3. Harmony integration to remove donor batch effects before DE
4. Focused analysis: TREM2+ macrophage sub-cluster within Monocytes/Kupffer cells
5. Visualise key genes with FeaturePlot: TREM2, CXCL13, HAVCR2, NNMT, POSTN

---

## Output Files

| File | Description |
|---|---|
| `scripts/01_download_setup.R` | GEO metadata download + human organism filter |
| `scripts/02_nafld_subset.R` | NAFLD-specific sample subsetting |
| `scripts/03_seurat_pipeline.R` | Seurat pipeline: QC → clustering → markers |
| `scripts/04_composition_and_de.R` | Cell-type annotation, composition, within-cell-type DE |
| `results/human_sample_metadata.csv` | All 24 human GSM metadata |
| `results/nafld_subset_metadata.csv` | 15-library NAFLD subset metadata |
| `results/nafld_seurat_clustered.rds` | Clustered Seurat object (not in git) |
| `results/nafld_seurat_annotated.rds` | Annotated Seurat object (not in git) |
| `results/all_markers.csv` | FindAllMarkers output (cluster markers) |
| `results/top5_markers_per_cluster.csv` | Top 5 markers per cluster by log2FC |
| `results/cell_composition_major.csv` | Mean % per major cell type per disease group |
| `results/within_celltype_de_all.csv` | All DE genes tested (healthy vs NAFLD), 18 cell types |
| `results/within_celltype_de_sig.csv` | Significant DE genes only (padj < 0.05) |
| `results/top10_de_per_celltype.csv` | Top 10 up + down per cell type |
| `plots/umap_cell_types.png` | UMAP: annotated cell type labels |
| `plots/composition_stacked_by_donor.png` | Stacked bar: per-donor cell type composition |
| `plots/composition_healthy_vs_nafld.png` | Dodged bar: mean % healthy vs NAFLD + SD |
| `plots/nafld_genes_dotplot.png` | Dot plot: key NAFLD genes across cell types |
| `plots/qc_violin_before.png` | QC violin plots before filtering |
| `plots/qc_violin_after.png` | QC violin plots after filtering |
| `plots/variable_features.png` | Top variable features |
| `plots/pca_elbow.png` | Elbow plot (PC selection) |
| `plots/pca_donor.png` | PCA coloured by donor |
| `plots/umap_clusters.png` | UMAP: cluster numbers |
| `plots/umap_disease.png` | UMAP: healthy vs NAFLD |
| `plots/umap_donor.png` | UMAP: coloured by donor |
| `plots/umap_fraction.png` | UMAP: CD45+/CD45− fraction |
| `plots/umap_panel.png` | 2×2 UMAP panel |
| `plots/marker_heatmap.png` | Heatmap: top 5 markers per cluster |
| `data/GSM404XXXX/` | Raw 10x MEX files — not tracked in git (380 MB) |
