# Week 4 Day 3 — NAFLD scRNA-seq: Human Liver Healthy vs Cirrhosis

## Dataset

**GEO Accession**: GSE136103  
**Publication**: Ramachandran et al. (2019) "Resolving the fibrotic niche of human liver cirrhosis at single-cell level." *Nature*, 575, 512–518.  
**Platform**: GPL20301 — Illumina HiSeq 4000 (Homo sapiens). GPL21103 (mouse) excluded.  
**Technology**: 10x Genomics Chromium single-cell RNA-seq  
**Tissue**: Human liver, sorted into CD45+ (immune-enriched) and CD45− (parenchymal-enriched) fractions  

### Samples used (liver only)

| Condition | Donors | GSM IDs | Fractions |
|-----------|--------|---------|-----------|
| Healthy   | 5      | GSM4041150 – GSM4041160 | CD45+, CD45−A, CD45−B per donor |
| Cirrhosis | 5      | GSM4041161 – GSM4041169 | CD45+, CD45− per patient |

**Excluded**: GSM4041170–GSM4041173 (peripheral blood, 4 samples) — this analysis focuses on liver tissue only.

## Downsampling Strategy

After merging all 20 liver fractions the object contained **61,900 cells** (35,512 Healthy / 26,388 Cirrhosis). After QC filtering (nFeature 200–5000, nCount 500–25000, percent.MT < 20%), **59,845 cells** passed.

**Stratified random downsampling** was applied to cap at **20,000 cells** while preserving condition proportions:

- Healthy: 34,705 → **11,598 cells** (57.99%)  
- Cirrhosis: 25,140 → **8,402 cells** (42.01%)

Sampling was done proportionally using `slice_sample()` per condition group (`set.seed(42)`).

## Pipeline Summary

| Step | Tool | Parameters |
|------|------|------------|
| Load | `ReadMtx()` | feature.column=2 (gene symbol), min.cells=3, min.features=200 |
| Normalise | `LogNormalize` | scale.factor = 10,000 |
| HVG | `vst` | nfeatures = 2,000 |
| PCA | `RunPCA` | npcs = 50 |
| Neighbours | `FindNeighbors` | dims = 1:30 |
| Clustering | `FindClusters` | resolution = 0.5 → 24 clusters |
| UMAP | `RunUMAP` | dims = 1:30, UWOT cosine |

## Cell Type Annotation

Clusters were annotated by scoring mean normalised expression of canonical liver markers across clusters:

| Cell Type | Markers Used |
|-----------|--------------|
| Hepatocyte | ALB, APOE, CYP3A4, TF |
| Macrophage | TREM2, CD68, MARCO, VSIG4 |
| Hepatic Stellate | ACTA2, COL1A1, COL1A2, PDGFRB |
| Endothelial | PECAM1, VWF, CLEC4M, LYVE1 |
| T Cell | CD3D, CD3E, CD8A, CD4 |
| B Cell | CD79A, MS4A1, CD19 |
| NK Cell | GNLY, NKG7, KLRD1 |
| Cholangiocyte | EPCAM, KRT19, KRT7 |

All 29 queried markers were found in the dataset. The 24 Seurat clusters resolved into 8 cell types (some types assigned to multiple clusters).

## Key Findings: Healthy vs Cirrhosis

### Cell type proportions

| Cell Type | Healthy (%) | Cirrhosis (%) | Direction |
|-----------|-------------|---------------|-----------|
| NK Cell | **42.9** | 21.9 | ↓↓ depleted in cirrhosis |
| T Cell | 18.6 | **22.7** | ↑ expanded in cirrhosis |
| Hepatocyte | 13.4 | **19.4** | ↑ (enrichment in parenchymal fraction) |
| Endothelial | 9.0 | **19.6** | ↑↑ expanded — sinusoidal capillarisation |
| Macrophage | 8.6 | 7.8 | ~ similar |
| B Cell | 2.9 | **5.6** | ↑ biliary/portal infiltration |
| Stellate | **4.7** | 2.7 | ↓ (but ACTA2+ activated subset present) |
| Cholangiocyte | 0.0 | **0.3** | only detected in cirrhosis |

### Biological interpretation

1. **NK cell depletion in cirrhosis** (~43% → ~22%): NK cells are the dominant immune cell in healthy liver. Their depletion in cirrhosis is well established and linked to hepatic immune tolerance and progressive fibrosis (Bogdanos et al. 2012; Ramachandran 2019).

2. **Endothelial expansion** (~9% → ~20%): Suggests sinusoidal capillarisation and neoangiogenesis, hallmarks of liver fibrosis where liver sinusoidal endothelial cells (LSECs) lose fenestrae and become capillarised.

3. **T cell increase**: Expansion of intrahepatic T cells in cirrhosis consistent with chronic inflammatory infiltrate. Contains both cytotoxic (CD8+) and helper (CD4+) subsets.

4. **TREM2+ macrophages**: Although overall macrophage proportion is similar, TREM2 expression (a key marker of scar-associated macrophages, SAMs) is enriched in cirrhosis — visible in the split UMAP and feature plots. SAMs are a key fibrogenic cell population discovered in Ramachandran 2019.

5. **Stellate cell activation**: ACTA2+ hepatic stellate cells (HSCs) are the primary fibrogenic cell type. Although proportions appear lower in cirrhosis (sampling artifact from CD45 enrichment), the ACTA2/COL1A1 expressing stellate clusters are predominantly in cirrhotic samples.

6. **Cholangiocyte appearance in cirrhosis**: Ductular reaction (biliary epithelial expansion) is a hallmark of advanced liver disease — EPCAM+ cholangiocytes appear only in cirrhotic samples.

## Files Generated

### Scripts
- `scripts/nafld_scrna_liver.R` — complete Seurat pipeline

### Plots
| File | Description |
|------|-------------|
| `01_qc_before_filtering.png` | Violin plots: nFeature, nCount, %MT by condition |
| `02_variable_features.png` | HVG plot with top 10 labelled |
| `03_pca_elbow.png` | PCA variance explained (30 PCs used) |
| `04_umap_clusters_condition.png` | UMAP by cluster / condition / CD45 fraction |
| `05_dotplot_markers.png` | DotPlot: 29 markers × 24 clusters |
| `06_feature_plots_markers.png` | FeaturePlots: ALB, TREM2, CD68, ACTA2, PECAM1, CD3D |
| `07_umap_cell_types.png` | UMAP coloured by annotated cell type |
| `08_umap_split_condition.png` | UMAP split by Healthy vs Cirrhosis |
| `09_celltype_proportion_stacked.png` | Stacked bar: cell type composition by condition |
| `10_celltype_proportion_grouped.png` | Grouped bar: direct proportion comparison |
| `11_trem2_cd68_by_condition.png` | TREM2 / CD68 split by condition |
| `12_heatmap_top_markers.png` | Heatmap: top 3 markers per cell type |

### Results Tables
| File | Description |
|------|-------------|
| `sample_metadata.csv` | GEO sample metadata with condition labels |
| `downsampling_summary.csv` | Cell counts per condition/sample after downsampling |
| `cluster_celltype_map.csv` | Cluster → cell type annotation |
| `celltype_proportions.csv` | Cell type % per condition |
| `celltype_proportions_per_donor.csv` | Cell type % per individual donor |
| `celltype_marker_genes.csv` | All marker genes per cell type (Wilcoxon) |
| `top5_markers_per_celltype.csv` | Top 5 markers per cell type by log2FC |
| `cell_counts_summary.csv` | Cell counts per cell type per condition |

### Data
- `data/liver_seurat_annotated.rds` — final Seurat object with UMAP, clusters, and cell type labels

## Notes on Data Structure

- GEO submission uses **non-standard filename prefixes** (e.g. `GSM4041150_healthy1_cd45+_matrix.mtx`), requiring `ReadMtx()` instead of `Read10X()`.
- Gene names contain underscores which Seurat converts to dashes automatically.
- Each donor has multiple fractions (CD45+/CD45−) which are merged before analysis.
- Blood samples (GSM4041170–4041173, `blood1`–`blood4`) were excluded from all analyses.

## Reference

Ramachandran, P. et al. (2019). Resolving the fibrotic niche of human liver cirrhosis at single-cell level. *Nature*, 575(7783), 512–518. https://doi.org/10.1038/s41586-019-1631-3
