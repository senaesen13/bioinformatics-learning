# Coexpression Analysis — 139 Shared NAFLD Genes

This folder contains a coexpression analysis of the 139 genes that are significantly
dysregulated in both Dataset 1 (GSE162694, discovery) and Dataset 2 (GSE135251,
validation). These 139 genes are the primary output of the cross-cohort comparison in
`week4-day4-nafld-cross-cohort-comparison/results/pairwise_comparison.md` (Comparison C).

The goal is to ask: among these 139 genes, which ones move together across samples?
Genes that co-vary across patients likely share regulatory programs or cell-of-origin.

---

## Method

**Expression data:** VST-normalised counts from Dataset 1 (GSE162694, n=143 samples).
Dataset 1 is used because it is the larger discovery cohort (31 normal + 112 NAFLD)
and has a good dynamic range for coexpression. VST (variance-stabilising transformation)
from DESeq2 makes the expression values suitable for correlation analysis.

**Steps:**
1. Compute the 139-gene overlap (gene_symbol intersection of Day 1 and Day 2
   significant_genes.csv; threshold: padj < 0.05 AND |MLE log2FC| > 1).
2. Extract VST-normalised expression for these 139 genes across all 143 samples.
3. Compute a 139×139 Pearson correlation matrix (gene-by-gene, across samples).
4. Build a distance matrix: distance = 1 − Pearson r. This puts highly correlated
   gene pairs close together and anti-correlated gene pairs far apart.
5. Hierarchical clustering (average linkage) on the distance matrix.
6. Cut the dendrogram into k clusters. k = 7 was chosen by silhouette analysis
   (best average silhouette width = 0.311 across k = 2–8).

**Note on k = 7:** The 139 genes are all part of the same broad NAFLD signature,
so they are generally co-upregulated and do not split into many distinct programs.
At k = 7, the two main modules (M1 and M2, together 132 genes) are the meaningful
structure. The remaining 7 genes fall into singleton or small clusters (M3–M7) that
do not cleanly belong to either main group.

---

## Results

### Module Sizes

| Module | Genes | Mean intra-module r | Biological character |
|---|---|---|---|
| **M1** | **62** | **0.618** | ECM / fibrosis / matrix remodelling |
| **M2** | **70** | **0.507** | Macrophage / immune / stress response |
| M3 | 1 | — | Outlier |
| M4 | 3 | 0.367 | Small outlier cluster |
| M5–M7 | 1 each | — | Outliers |

M1 and M2 together contain 132 of the 139 genes. They are the focus of the analysis.

### M1 — ECM / Fibrosis Module (62 genes, mean r = 0.618)

Selected key members: **COL1A1**, COL6A2, EMILIN1, ADAMTSL2, LAMA5, LTBP4, MEGF6,
FASN, CLEC11A, CTSD, CTSG, DERL3, NOTUM, NLRP6, S1PR2, MIF, CCND1.

COL1A1 is the canonical fibrosis marker and collagen-I gene. COL6A2, EMILIN1, LAMA5,
LTBP4 are all extracellular matrix structural or remodelling proteins. The tight
co-clustering of these genes (mean r = 0.618) reflects their shared source — hepatic
stellate cells, which are the primary producers of ECM in NAFLD fibrosis. This module
likely represents the stellate cell / fibrosis transcriptional program that is
consistently activated across NAFLD patients in Dataset 1.

### M2 — Macrophage / Immune / Stress-Response Module (70 genes, mean r = 0.507)

Selected key members: **TREM2**, **SPP1**, CCL2, CCL24, FCN1, C1QC, CXCR4, CHI3L1,
MMP9, THBS1, SERPINE1, FOS, FOSB, JUN, JUNB, ATF3, DUSP1, DUSP2, EGR2, NR4A1,
KLF4, MAFF, SOCS3, CDKN1A, GADD45B, GADD45G, HBEGF, MYC.

TREM2 and SPP1 co-cluster in M2, consistent with their established roles as co-expressed
markers of lipid-associated macrophages (LAMs) in NAFLD. CCL2 and CCL24 are macrophage
recruitment chemokines. FCN1 and C1QC are complement pathway components expressed by
liver-resident Kupffer cells. CHI3L1 is a hepatic inflammation marker.

The stress-response transcription factors (FOS, FOSB, JUN, JUNB, ATF3, DUSP1, DUSP2,
EGR2, NR4A1, KLF4) co-cluster with the macrophage genes. This reflects the known
biology: macrophage activation in NAFLD is driven by danger-signal and lipotoxic stress
signalling, which induces these immediate-early transcription factors. These genes likely
co-vary because they are all responsive to the same inflammatory/lipotoxic signals in
the liver microenvironment.

### Key Gene Summary

| Gene | Module | Interpretation |
|---|---|---|
| TREM2 | **M2** | Macrophage/immune program — co-clusters with SPP1 |
| SPP1 | **M2** | Macrophage/immune program — co-clusters with TREM2 |
| COL1A1 | **M1** | Fibrosis/ECM program — separate from macrophage genes |
| GPNMB | *not in overlap* | Significant in Day 1 but not Day 2 (LFC +1.19, padj 2.3e-09 in Day 1; not significant in Day 2) |

**TREM2 and SPP1 are in the same module.** This is consistent with scRNA-seq data
(see `week6-day1-nafld-scrna-gse136103/`) showing they are co-expressed in the
TREM2+ lipid-associated macrophage subpopulation. The bulk coexpression analysis
independently recovers this cell-type-specific co-regulation signal.

**COL1A1 is in a separate module from TREM2/SPP1.** This reflects their different
cellular origins: TREM2/SPP1 come from macrophages, COL1A1 comes from hepatic stellate
cells. These two programs both contribute to NAFLD pathology but through different
cell types and regulatory mechanisms.

---

## Interpretation

The 139 shared NAFLD genes split into two main functional programs:

1. **Fibrosis/ECM program (M1):** Driven by stellate cell activation and matrix
   production. COL1A1 is the anchor gene.

2. **Macrophage/immune program (M2):** Driven by macrophage infiltration and activation,
   plus lipotoxic stress-response signalling. TREM2 and SPP1 are the anchor genes.

This two-program structure makes biological sense. NAFLD progression involves both
immune infiltration and fibrotic scarring, and these two processes involve different
cell types. The coexpression analysis shows that the 139 reproducible genes capture
both programs, and that they are largely separable at the expression level.

The mean intra-module correlation is moderate (r = 0.62 for M1, r = 0.51 for M2)
rather than high. This is expected for bulk RNA-seq data: the signal from each cell
type is diluted by the many other cell types in the bulk liver biopsy. Single-cell
data would give much tighter co-regulation.

---

## Output Files

| File | Description |
|---|---|
| `scripts/coexpression_analysis.R` | Full analysis pipeline |
| `results/vst_matrix_139genes.csv` | VST expression of 139 genes × 143 samples |
| `results/correlation_matrix.csv` | 139 × 139 Pearson correlation matrix |
| `results/module_assignments.csv` | Gene → module mapping (k = 7) |
| `plots/coexpression_heatmap.png` | Correlation heatmap, genes ordered by module |
