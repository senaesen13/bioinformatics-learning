# Week 3 Day 3 — Spatial Transcriptomics: Human Heart (10x Visium)

> **Script:** `scripts/spatial_heart.R`
> **Dataset:** V1_Human_Heart — 10x Genomics Visium, healthy human heart section
> **Run from:** `week3-day3-spatial-transcriptomics-heart/` as `Rscript scripts/spatial_heart.R`

---

## What is spatial transcriptomics?

Bulk RNA-seq gives one expression value per gene per sample. scRNA-seq gives one profile per cell but loses where in the tissue each cell was. Spatial transcriptomics keeps both: you get gene expression **at known physical coordinates** across a tissue section.

10x Visium works by placing a tissue section on a slide printed with ~5,000 spots. Each spot (~55 µm diameter) captures RNA from the 1–10 cells sitting on top of it. You get a count matrix of genes × spots, plus a photograph of the tissue, plus the (x, y) coordinates of every spot. The result: a gene expression map overlaid on the tissue image.

---

## Dataset

**V1_Human_Heart** — 10x Genomics Visium spatial gene expression

- Tissue: healthy human heart, cryosection
- Spots: 4,247 detected (4,233 after QC)
- Genes: 36,601
- Expression: UMI counts per spot (not per cell)
- Files: `filtered_feature_bc_matrix/` (count matrix) + `spatial/` (image + coordinates)

---

## Pipeline steps

| Step | What it does |
|------|-------------|
| **1. Load** | `Read10X()` reads the count matrix; `Read10X_Image()` loads the tissue photograph and spot coordinates; combined into a Seurat object |
| **2. QC** | The total spot count is recorded before any filtering. Spots are then kept if nCount > 200, nFeature > 100, and percent.mt < 60%. Heart tissue is naturally ~30–60% mitochondrial — a strict MT threshold would remove real cardiac spots. Result: 4,247 → 4,233 spots |
| **3. SCTransform** | Normalises UMI counts for sequencing depth differences between spots using a regularised negative binomial regression. Better than log-normalisation for spatial data |
| **4. PCA + UMAP** | 30 PCs computed on SCT-normalised data; UMAP for 2D visualisation |
| **5. Clustering** | KNN graph (30 PCs) → Louvain at resolution 0.5 → 7 spatial clusters |
| **6. Spatially variable genes** | `FindSpatiallyVariableFeatures()` with Moran's I: finds genes whose expression is spatially autocorrelated (varies across the tissue in a non-random pattern) |
| **7. Marker genes** | `FindAllMarkers()` Wilcoxon test: genes that distinguish each cluster from all others |
| **8. Cardiac markers** | Spatial maps of known cell-type markers: TNNT2, MYH7, MYH6, VIM, DCN, PECAM1, ACTA2, CD68 |

---

## Key results

**7 spatial clusters** across 4,233 QC-passed spots:

| Cluster | Spots | % |
|---------|-------|---|
| 0 | 1,187 | 28% |
| 1 | 1,155 | 27% |
| 2 | 1,124 | 27% |
| 3 | 596 | 14% |
| 4 | 69 | 2% |
| 5 | 65 | 2% |
| 6 | 37 | 1% |

**Top spatially variable genes** (highest Moran's I — most spatially structured expression):

`ACTA1, MYL2, MYL7, MB, MYH7, MT-ND3, MT-CO2, IGLC2, MT-CO3, MT-CO1, FOS, JUNB`

Biologically coherent: sarcomeric proteins (ACTA1, MYL2, MYL7, MYH7), myoglobin (MB), mitochondrial respiratory chain (MT-CO1/2/3, MT-ND3), immune (IGLC2), and immediate-early stress response (FOS, JUNB).

**All 8 canonical cardiac markers found in the dataset:** TNNT2, MYH7, MYH6, VIM, DCN, PECAM1, ACTA2, CD68 — allowing spatial mapping of cardiomyocytes, fibroblasts, endothelium, smooth muscle, and macrophages across the section.

---

## Output files

| File | Contents |
|------|---------|
| `plots/01_qc_spatial.png` | UMI count, gene count, % mitochondrial per spot on tissue |
| `plots/02_qc_violin.png` | Violin plots of QC metrics |
| `plots/03_pca_elbow.png` | Variance explained per PC |
| `plots/04_umap_clusters.png` | UMAP coloured by cluster |
| `plots/05_spatial_clusters.png` | Cluster map on tissue image |
| `plots/06_spatially_variable_genes.png` | Top 6 spatially variable genes on tissue |
| `plots/07_marker_dotplot.png` | Top 2 marker genes per cluster |
| `plots/08_top_marker_spatial.png` | Spatial maps of top cluster markers |
| `plots/09_cardiac_markers_spatial.png` | TNNT2, MYH7, MYH6, VIM, DCN, PECAM1, ACTA2, CD68 on tissue |
| `plots/10_cardiac_markers_violin.png` | Cardiac markers across clusters |
| `results/spatially_variable_genes.csv` | All SVGs ranked by Moran's I |
| `results/cluster_markers_all.csv` | All marker genes with statistics |
| `results/cluster_markers_top5.csv` | Top 5 markers per cluster |
| `results/cluster_sizes.csv` | Spot counts per cluster |

---

## Connection to previous weeks

- **vs scRNA-seq (Week 3 Day 1):** Seurat pipeline is nearly identical — same QC, normalisation, PCA, clustering steps. Key difference: you also plot onto the tissue image, and `FindSpatiallyVariableFeatures` replaces `FindAllMarkers` as the primary discovery tool.
- **vs bulk RNA-seq (Week 2):** Bulk would average all 4,233 spots into one profile. Spatial preserves the geography — you can see that MYH7 (ventricular myosin) is enriched in one region while MYH6 (atrial) is in another.
- **vs drug repositioning (Week 2–3):** Spatially variable genes like ACTA1 and MYL2 could be used as disease-relevant inputs for GSEA/CMap analysis to find drugs targeting specific cardiac regions.

---

## References

1. **Dataset** — 10x Genomics V1_Human_Heart Visium dataset.
   https://www.10xgenomics.com/datasets/human-heart-1-standard-1-1-0

2. **Seurat spatial vignette** — Satija Lab.
   https://satijalab.org/seurat/articles/spatial_vignette

3. **Moran's I** — Moran PAP. (1950). "Notes on Continuous Stochastic Phenomena." *Biometrika* 37(1):17–23.

4. **SCTransform** — Hafemeister C, Satija R. (2019). "Normalization and variance stabilization of single-cell RNA-seq data using regularized negative binomial regression." *Genome Biology* 20:296.
   https://doi.org/10.1186/s13059-019-1874-1

---

*Week 3, Day 3 — Spatial transcriptomics pipeline*
*Follows from: Week 3 Day 1 (Seurat scRNA-seq PBMC) and Week 3 Day 2 (ccRCC drug repositioning)*
