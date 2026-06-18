# Week 3 Day 1 — scRNA-seq: Seurat PBMC 3k

> **Script:** `scripts/seurat_pbmc3k.R`
> **Dataset:** 10x Genomics PBMC 3k — 2,700 peripheral blood mononuclear cells from a healthy donor
> **Run from:** `week3-day1-scrna-seq/` as `Rscript scripts/seurat_pbmc3k.R`

---

## What is this dataset?

The PBMC 3k dataset is the standard learning dataset for scRNA-seq — the equivalent
of the airway dataset in bulk RNA-seq. It was released by 10x Genomics and is used
in virtually every Seurat tutorial and scRNA-seq course.

**PBMC** = Peripheral Blood Mononuclear Cells. These are the immune cells in blood
that have a round nucleus — mainly T cells, B cells, NK cells, and monocytes.
Blood is a good starting point for scRNA-seq because:
- Easy to collect (unlike heart or brain tissue)
- Well-characterised cell types with known markers
- No tissue dissociation stress (cells are already in suspension)

---

## Pipeline steps

### Step 1 — Load data
`Read10X()` reads the three files that make up a 10x count matrix:
- `barcodes.tsv` — one cell barcode per row
- `genes.tsv` — one gene per row
- `matrix.mtx` — the sparse UMI count matrix

`CreateSeuratObject()` wraps these into a Seurat object with initial filters:
`min.cells = 3` (keep genes seen in ≥3 cells) and `min.features = 200` (keep
cells with ≥200 detected genes).

**Result:** 2,700 cells × 13,714 genes

---

### Step 2 — Quality Control (QC)

Three metrics identify bad cells:

| Metric | Too low | Too high |
|--------|---------|----------|
| `nFeature_RNA` (genes per cell) | Empty droplet / dead cell | Doublet (two cells) |
| `nCount_RNA` (UMIs per cell) | Poor capture | Doublet |
| `percent.mt` (% mitochondrial) | — | Dead/dying cell |

**Why mitochondrial %?** When a cell's membrane is damaged, cytoplasmic RNA leaks
out. Mitochondria (membrane-enclosed organelles) stay intact. Dead cells therefore
have disproportionately high mitochondrial RNA relative to total RNA.

Thresholds applied: `nFeature_RNA` 200–2500, `percent.mt` < 5%.

**Result:** 2,700 → **2,638 cells** after QC (62 low-quality cells removed)

---

### Step 3 — Normalisation

Different cells capture different amounts of RNA due to technical noise, not biology.
`NormalizeData()` corrects this:

1. Divide each cell's counts by its total UMIs × 10,000 (= counts per 10k, like CPM)
2. Apply log1p: `log(value + 1)`

After normalisation, cells are on the same scale and the data is less skewed.

---

### Step 4 — Highly Variable Genes (HVGs)

Not all 13,714 genes are informative. `FindVariableFeatures()` selects the
**2,000 genes** that vary most across cells relative to their expression level.
Housekeeping genes (similar in every cell) are excluded.

Only these 2,000 HVGs are used for PCA and clustering.

Top HVGs: `PPBP, LYZ, S100A9, IGLL5, GNLY, FTL, PF4, FTH1, GNG11, S100A8`

---

### Step 5 — Scaling

`ScaleData()` z-scores each gene across cells (mean = 0, variance = 1).
This prevents high-expression genes from dominating PCA simply because of their
magnitude rather than their biological variation.

---

### Step 6 — PCA

`RunPCA()` compresses 2,000 HVGs into 50 principal components. Each PC captures
a dimension of variation:
- PC1 separates monocytes from lymphocytes (top genes: CST3, TYROBP, LST1)
- PC2 separates T cells from B cells
- etc.

The **elbow plot** (`plots/04_elbow_plot.png`) shows variance explained per PC.
The curve flattens around PC10 — meaning the first 10 PCs capture most of the
real signal. PCs beyond ~15 mainly capture noise.

**We use dims = 1:10 for all downstream steps.**

---

### Step 7 — Clustering

Two sub-steps:

1. **`FindNeighbors(dims = 1:10)`** — builds a K-nearest-neighbour graph in PCA space.
   Each cell is connected to its 20 most similar neighbours. Similar cells form
   dense, connected communities in this graph.

2. **`FindClusters(resolution = 0.5)`** — Louvain community detection finds groups
   of cells that are more connected internally than to the rest of the graph.
   Resolution 0.5 gives 9 clusters for this dataset — a good granularity for
   identifying the major PBMC cell types.

**Result:** 9 clusters (sizes: 684, 481, 476, 344, 291, 162, 155, 32, 13 cells)

---

### Step 8 — UMAP

`RunUMAP(dims = 1:10)` projects the 10-dimensional PCA embedding into 2D for
visualisation. Cells that are neighbours in PCA space end up close on the UMAP.

**Important:** UMAP is only for visualisation. Cluster distances on UMAP are not
meaningful for statistics — always use PCA space or normalised counts for any
quantitative analysis.

---

### Step 9 — Marker Genes

`FindAllMarkers()` runs a Wilcoxon rank-sum test comparing each cluster against
all others. Parameters:
- `only.pos = TRUE` — only upregulated markers
- `min.pct = 0.25` — gene expressed in ≥25% of cluster cells
- `logfc.threshold = 0.25` — minimum fold change

The dot plot (`plots/08_marker_dotplot.png`) shows canonical PBMC markers across
clusters and is the key plot for deciding cell type assignments.

---

### Step 10 — Cell Type Annotation

Clusters are labelled using canonical PBMC marker genes:

| Cluster | n cells | Key markers | Cell type |
|---------|---------|-------------|-----------|
| 0 | 684 | IL7R+, CCR7+ | CD4+ Naive T |
| 1 | 481 | CD14+, LYZ+, CST3+ | CD14+ Monocytes |
| 2 | 476 | IL7R+, S100A4+ | CD4+ Memory T |
| 3 | 344 | MS4A1+, CD79A+ | B cells |
| 4 | 291 | CD8A+, NKG7+ | CD8+ T |
| 5 | 162 | FCGR3A+, MS4A7+ | CD16+ Monocytes |
| 6 | 155 | GNLY++, NKG7++ | NK cells |
| 7 | 32 | FCER1A+, FCER1A+, CST3+ | Dendritic cells |
| 8 | 13 | PPBP++ | Platelets |

---

## Final cell type composition

| Cell type | Cells | % of total |
|-----------|-------|-----------|
| CD4+ Naive T | 684 | 25.9% |
| CD14+ Monocytes | 481 | 18.2% |
| CD4+ Memory T | 476 | 18.0% |
| B cells | 344 | 13.0% |
| CD8+ T | 291 | 11.0% |
| CD16+ Monocytes | 162 | 6.1% |
| NK cells | 155 | 5.9% |
| Dendritic cells | 32 | 1.2% |
| Platelets | 13 | 0.5% |
| **Total** | **2,638** | |

---

## Output files

| File | Contents |
|------|---------|
| `plots/01_qc_violin.png` | Violin plots of QC metrics per cell |
| `plots/02_qc_scatter.png` | nCount vs percent.mt / nCount vs nFeature |
| `plots/03_variable_features.png` | HVG selection — top 10 labelled |
| `plots/04_elbow_plot.png` | PCA variance explained, elbow at ~PC10 |
| `plots/05_pca_heatmap.png` | Top genes driving PCs 1–6 |
| `plots/06_umap_clusters.png` | UMAP coloured by cluster number |
| `plots/07_marker_heatmap.png` | Heatmap of top 5 markers per cluster |
| `plots/08_marker_dotplot.png` | Dot plot of canonical PBMC markers |
| `plots/09_umap_annotated.png` | UMAP coloured by cell type |
| `plots/10_umap_comparison.png` | Side-by-side: clusters vs annotation |
| `results/cell_type_counts.csv` | Cells per annotated cell type |
| `results/cluster_markers_all.csv` | All marker genes with statistics |
| `results/cluster_markers_top5.csv` | Top 5 markers per cluster |
| `results/pbmc3k_seurat.rds` | Full Seurat object (reload with `readRDS()`) |

---

## Connection to Week 2 bulk RNA-seq

In Week 2, we analysed bulk RNA-seq from mouse heart tissue (MI vs sham) — one
expression value per gene per sample. Here, we have one expression profile per
cell. The key difference in practice:

- DESeq2 (bulk): compare two conditions, find differentially expressed genes
- Seurat (scRNA-seq): cluster cells, find what makes each cluster unique

The tools overlap more than they seem:
- Both use normalisation before comparison
- Both use log-transformation
- Both produce marker genes (DESeq2 via differential expression; FindAllMarkers via
  Wilcoxon test)
- The GSEA pathway analysis from Week 2 Day 2 can be applied to scRNA-seq marker
  genes to understand what each cell type is biologically doing

---

*Week 3, Day 1 — Seurat PBMC 3k on 10x Genomics data*
*Follows from: Week 2 Day 5 scRNA-seq theory notes*
