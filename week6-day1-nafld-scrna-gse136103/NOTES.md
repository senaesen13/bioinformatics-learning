# Week 6 Day 1 — NAFLD scRNA-seq: GSE136103

Single-cell RNA-seq of human liver from Ramachandran et al. 2019 (*Nature*, GSE136103), analysed here for 5 healthy donors and 2 NAFLD-cirrhosis donors (Cirrhotic1 + Cirrhotic4) across 35,050 cells and 21 clusters. The headline finding: TREM2 and SPP1/GPNMB — genes elevated in our bulk NAFLD cohorts — are specifically expressed in Kupffer cells and monocytes at single-cell resolution, and COL1A1 localises to hepatic stellate cells, directly confirming and extending the bulk RNA-seq signal.

---

## Dataset

**GEO:** GSE136103 | **Paper:** Ramachandran et al. 2019, *Nature* | **Species:** Homo sapiens only

| Group | Donors | Libraries | GSM range |
|---|---|---|---|
| Healthy liver | Healthy1–5 | 11 (CD45+/CD45−) | GSM4041150–160 |
| NAFLD cirrhosis | Cirrhotic1, Cirrhotic4 | 4 | GSM4041161–163, GSM4041168 |
| **Total** | **7** | **15** | |

Excluded: Cirrhotic2/3 (alcohol), Cirrhotic5 (PBC), Blood1–4 (PBMC) — 9 libraries dropped before analysis.

---

## Pipeline

| Step | What was done | Key numbers |
|---|---|---|
| Load | `ReadMtx()` × 15 libraries, `JoinLayers()` after merge | 46,331 cells, 22,517 features |
| QC | nFeature 200–2,500, percent.mt < 5% | 35,050 retained (75.7%) |
| Normalise | LogNormalize, scale factor 10,000 | — |
| Variable features | VST, top 2,000 | — |
| PCA | 50 PCs computed, 20 used | Elbow plot: `plots/pca_elbow.png` |
| Cluster + UMAP | FindNeighbors/Clusters res 0.5, UWOT | **21 clusters** |
| Markers | FindAllMarkers, Wilcoxon, min.pct 0.25 | 12,313 significant rows |

Scripts: `scripts/03_seurat_pipeline.R` (clustering) · `scripts/04_composition_and_de.R` (annotation, composition, DE)

---

## Cell Types Found (21 clusters → 20 labels)

Key plots: `plots/umap_annotated.png` · `plots/umap_healthy_vs_nafld.png`

| Cluster(s) | Cell type | Key markers |
|---|---|---|
| 0, 3 | CD4+ T cells | CD40LG, IL7R, LTB |
| 2 | CD8+ T cells | CD8A, CD8B, CD3G |
| 4 | CD8+ T cells (exhausted) | CD8A, LAG3, CRTAM |
| 16 | Naive T cells | MAL, LEF1, CCR7 |
| 1 | NK/NKT cells | TOX2, XCL1, KLRC1 |
| 5 | NK cells (cytotoxic) | GNLY, GZMB, FGFBP2 |
| 14 | gdT/NK cells | GZMH, TRGC2, CX3CR1 |
| 18 | NK cells (liver-resident) | IL2RB, CD160, NCR1 |
| 6 | Monocytes | S100A8/9/12, FCN1 |
| 10 | Dendritic cells | CD1C, FCER1A, CLEC10A |
| 11 | Kupffer cells | CD5L, C1QB/C, CD163, GPNMB |
| 8 | B cells | CD79A, CD19 |
| 15 | Plasma cells | IGHGP, IGLL5, IGHA2 |
| 7 | Endothelial cells | GPIHBP1, PODXL, AQP7 |
| 9 | LSEC | CLEC4G, FCN2/3, OIT3 |
| 12 | Hepatic stellate cells | DCN, TCF21, RGS5, COL1A1 |
| 13 | Hepatocytes | UGT2B15, UGT2A3, SPP1 |
| 17 | Proliferating cells | CENPA, RRM2, TYMS |
| 19 | Cholangiocytes | SCT, PTCRA, LRRC26 |
| 20 | Mast cells | TPSAB1, TPSB2, CPA3 |

---

## Bulk RNA-seq Genes Resolved to Cell Types

Key plots: `plots/dotplot_key_genes.png` · `plots/featureplot_trem2_cd9.png` · `plots/featureplot_col1a1_gpnmb.png`

Genes elevated in NAFLD across all three bulk cohorts (Week 4) are now traceable to specific cell types:

| Gene | Bulk finding | scRNA cluster marker | NAFLD vs healthy DE (within cell type) |
|---|---|---|---|
| TREM2 | Upregulated NAFLD | Not a cluster marker (expressed in <25% of any cluster) | Monocytes log2FC 3.51 padj 1.2e-26; Kupffer cells log2FC 2.56 padj 8.0e-24 |
| SPP1 | Upregulated NAFLD | Cluster 13 (Hepatocytes) log2FC 6.24 | Kupffer cells log2FC 2.62 padj 1.9e-10; Hepatocytes log2FC 2.22 padj 5.8e-30 |
| GPNMB | Upregulated NAFLD | Cluster 11 (Kupffer cells) log2FC 5.79 | HSC log2FC 2.44 padj 2.5e-02 |
| COL1A1 | Upregulated NAFLD | Cluster 12 (HSC) log2FC 7.89 padj 0 | HSC log2FC 3.55 padj 6.1e-16 |

**Reading the table:** TREM2 and GPNMB are macrophage-lineage signals (Kupffer cells + recruited monocytes). SPP1 marks hepatocytes at baseline but is additionally disease-upregulated in Kupffer cells. COL1A1 is a clean HSC fibrosis signal at both levels.

Additional NAFLD DE findings: CXCL13 log2FC 12.5 in exhausted CD8+ T cells (exhaustion signature with HAVCR2/TIM-3 log2FC 4.0); NNMT log2FC 3.66 in hepatocytes (metabolic stress); POSTN log2FC 4.18 in endothelial cells (fibrosis ECM).

Full DE results: `results/within_celltype_de_sig.csv` (8,913 genes, padj < 0.05 across 18 cell types)

---

## Cell Type Composition: Healthy vs NAFLD Cirrhosis

Key plot: `plots/composition_barplot.png` · data: `results/cell_composition_major.csv`

Mean % of total cells per donor group:

| Major cell type | Healthy (n=5) | NAFLD (n=2) | Fold |
|---|---|---|---|
| T cells | 43.8% | 46.4% | 1.06 |
| NK/NKT cells | 24.4% | 16.4% | 0.67 |
| B/Plasma cells | 3.7% | 7.6% | 2.08 |
| Monocytes | 7.6% | 10.7% | 1.40 |
| Endothelial/LSEC | 8.5% | 8.5% | 1.00 |
| Kupffer cells | 3.4% | 2.7% | 0.78 |
| Hepatic stellate cells | 2.7% | 1.5% | 0.53† |
| Dendritic cells | 3.1% | 1.7% | 0.54 |
| Hepatocytes | 1.3% | 2.3% | 1.84† |
| Proliferating cells | 0.8% | 0.8% | 1.02 |
| Cholangiocytes | 0.5% | 1.0% | 1.86† |
| Mast cells | 0.2% | 0.5% | 2.51 |

†Rows marked with † are distorted by the Cirrhotic4 CD45− gap (see Caveats).

---

## Caveats

**Cirrhotic4 CD45− gap.** GSM4041168 (Cirrhotic4) has only a CD45+ library — no CD45− fraction was deposited. Cirrhotic4 therefore contributes no hepatic stellate cells, hepatocytes, LSEC, or cholangiocytes to the merged object. Proportions of these cell types in the NAFLD group are underestimates.

**TREM2 cluster-marker threshold.** TREM2 does not appear in `all_markers.csv` because fewer than 25% of cells in any single cluster express it — the `min.pct = 0.25` cutoff in FindAllMarkers excludes it. It is confirmed disease-enriched via within-cell-type DE in both Monocytes and Kupffer cells (see table above), and should be visualised with `FeaturePlot` for spatial confirmation.

**XIST artefact.** XIST appears upregulated in NAFLD across multiple cell types, likely because the NAFLD donors include more female samples than the healthy group. Exclude sex-linked genes (XIST, RPS4Y1, DDX3Y, EIF1AY) before any downstream pathway analysis.

---

## Conclusion

scRNA-seq on this 7-donor NAFLD subset resolves the bulk signal: TREM2 and GPNMB mark macrophages (Kupffer cells + recruited monocytes), SPP1 marks both hepatocytes and NAFLD-activated macrophages, and COL1A1 cleanly localises to hepatic stellate cells — consistent with the fibrotic niche described in Ramachandran 2019.

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
| `plots/umap_annotated.png` | UMAP: annotated cell types, colour-coded palette (script 05) |
| `plots/umap_healthy_vs_nafld.png` | UMAP: split panels — healthy vs NAFLD cirrhosis (script 05) |
| `plots/featureplot_trem2_cd9.png` | FeaturePlot: TREM2 and CD9 overlaid on UMAP (script 05) |
| `plots/featureplot_col1a1_gpnmb.png` | FeaturePlot: COL1A1 and GPNMB overlaid on UMAP (script 05) |
| `plots/composition_barplot.png` | Bar plot: cell type proportions, healthy vs NAFLD (script 05) |
| `plots/dotplot_key_genes.png` | Dot plot: TREM2/SPP1/GPNMB/COL1A1 across all cell types (script 05) |
| `plots/umap_cell_types.png` | UMAP: annotated cell type labels (script 04, earlier version) |
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
