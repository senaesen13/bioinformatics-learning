# Week 8 Day 1: MASLD Spatial Transcriptomics — Full Analysis Pipeline

## Paper

**Title:** Progressive fibrosis in human MASLD is associated with spatially linked transcriptomic signatures of metabolic reprogramming and senescence  
**Authors:** Hani Vu, Yuliangzi Sun, Zherui Xiong, Xiao Tan, et al. (Quan H Nguyen, Elizabeth E Powell group, University of Queensland)  
**Journal:** JHEP Reports, 2025  
**PMID:** 41541503 | **DOI:** 10.1016/j.jhepr.2025.101657  
**Dataset DOI:** 10.48610/e95155f (UQeSpace) | **Code:** github.com/BiomedicalMachineLearning/Liver

---

## Dataset Overview

**Platform:** 10x Genomics Visium CytAssist (FFPE-compatible, probe-based)  
**Target genes:** 18,085 (Visium Human Transcriptome Probe Set v2.0 — not full transcriptome)  
**Design:** 4–5 MASLD liver biopsies co-spotted per Visium capture area (multiplexed)  
**Arrays used:** 8 (matching Vu et al. published analysis — excluding VLP115_D and VLP116_A)  
**Samples:** 33 liver biopsies from F0–F4 fibrosis stages; 3-group scheme: early (F0–F1, n=14), intermediate (F2–F3a, n=5), late (F3b–F4, n=14)

---

## Step 1 — Spot QC

**Approach:** Filtered per array using: ≥200 genes/spot, <10% mitochondrial reads, gene in ≥3 spots.

**Thresholds rationale:**
- 200 genes/spot: standard lower bound for spatial Visium data; CytAssist FFPE has lower UMI depth than fresh-frozen so a less strict cutoff than single-cell (200 vs 500) is appropriate
- 10% mito: standard liver tissue threshold; higher mito% indicates compromised/dying cells or FFPE extraction artifact
- 3 cells/gene: removes probe artefacts and very lowly detected transcripts

**QC results:**

| Array    | Before QC | After QC | % kept | Med genes/spot | Med UMIs/spot |
|----------|-----------|----------|--------|----------------|---------------|
| VLP115_A | 1,744     | 1,415    | 81.1%  | 1,379          | 2,487         |
| VLP116_D | 2,173     | 2,163    | 99.5%  | 2,600          | 5,562         |
| VLP119_A | 1,691     | 1,689    | 99.9%  | 2,493          | 5,727         |
| VLP119_D | 1,277     | 1,277    | 100%   | 2,531          | 5,519         |
| VLP120_A | 1,884     | 1,866    | 99.0%  | 4,024          | 11,377        |
| VLP120_D | 1,657     | 1,656    | 99.9%  | 2,721          | 6,233         |
| VLP121_A | 1,530     | 1,492    | 97.5%  | 3,970          | 11,500        |
| VLP121_D | 1,683     | 1,681    | 99.9%  | 2,563          | 5,751         |
| **Total**| **13,659**| **13,239**| **96.9%**| —           | —             |

**Key observation:** VLP115_A retained only 81.1% of spots (329 removed) and has the lowest median UMI depth (2,487/spot) across the 8 active arrays. This is the "borderline quality" array noted in the dataset documentation — it was included in the paper's GitHub code but may have contributed noisier data. All other arrays lost <3% of spots, confirming high overall quality.

The two excluded arrays (VLP115_D: 952 genes/spot; VLP116_A: not in published code) were correctly omitted.

---

## Step 2 — Normalization: SCTransform + Merge

**Approach:** Applied SCTransform (Hafemeister & Satija 2019) independently to each of the 8 arrays, then merged into a single combined Seurat object.

**Why SCTransform per array, not after merging?**  
SCTransform fits a regularized negative binomial regression per gene to remove sequencing depth effects. Fitting it per array avoids cross-array batch bias in model estimation. After per-array normalization, a shared set of 3,000 variable features was selected using `SelectIntegrationFeatures` (choosing genes that are highly variable across samples), and all 8 arrays were merged.

**Why not Harmony/CCA integration?**  
The paper's own workflow used SCTransform without batch correction beyond per-sample normalization. We follow the same approach. The multiplexed slide design (mixed fibrosis stages per array) was deliberately chosen to minimize batch effects from technical sources — biological variation is expected to dominate.

**Result:** Combined object = **13,239 spots × 17,910 genes**, 3,000 variable features selected.

Note: glmGamPoi (faster SCTransform backend) was not installed during the initial run, causing the "slower native implementation" fallback. glmGamPoi has since been installed for future use.

---

## Step 3 — Spatial Domain Discovery

**Approach:** Standard Seurat dimensionality reduction + unsupervised clustering on the merged SCT data.
1. **PCA:** 30 principal components computed on 3,000 variable features
2. **UMAP:** 20 PCs used (elbow plot showed eigenvalue leveling ~PC15–18)  
3. **Leiden/Louvain clustering:** `FindNeighbors` (dims=1:20) + `FindClusters` (resolution=0.5)

**Result: 24 spatial domain clusters** across 13,239 spots

| Cluster | Spots | % of total |
|---------|-------|------------|
| C0      | 2,521 | 19.0%      |
| C1      | 1,752 | 13.2%      |
| C2      | 1,061 | 8.0%       |
| C3      | 1,006 | 7.6%       |
| C4      | 833   | 6.3%       |
| C5–C10  | 2,857 | 21.6%      |
| C11–C20 | 2,560 | 19.3%      |
| C21–C23 | 260   | 2.0%       |

**What these clusters likely represent:**  
In a liver biopsy at the resolution of Visium (~55 µm spot diameter, ~3–10 cells per spot), the 24 clusters do NOT represent distinct cell types (liver is ~70–80% hepatocytes by cell number and ~90%+ by volume). Instead, they represent **different transcriptional states of the hepatocyte parenchyma**, reflecting:
- Metabolic zonation (periportal vs centrilobular hepatocytes)
- Fibrosis-associated activation programs
- Inflammatory niche microenvironments
- Biliary/portal vs parenchymal spatial zones

This matches what Vu et al. found — their spatial domains corresponded to distinct metabolic and senescence signatures rather than cell-type boundaries.

**Batch assessment:** The UMAP by array (saved in `results/spatial_domains/umap_clusters_and_arrays.png`) shows how spots from different arrays distribute across the 24 clusters — this indicates whether any cluster is dominated by a single array (technical batch effect) vs. mixed across arrays (biological signal).

---

## Step 4 — Cell-type Deconvolution

**Reference:** GSE136103 scRNA-seq (Ramachandran et al. 2019, human cirrhotic liver) — 22,517 genes × 35,050 cells, annotated into 20 cell types.

**Target method:** CARD (Ma & Zhou 2022, Nat Biotechnol) — a Bayesian deconvolution framework that models spatial autocorrelation to infer cell-type composition per spot.

**What happened:** CARD could not be compiled on this macOS system (C++ linker failure with clang). The pipeline fell back to **Seurat AddModuleScore**, which scores each spot for cell-type activity based on the average expression of top 50 marker genes per cell type from the reference.

**Marker gene overlap (reference → Visium probe panel):**

| Cell Type              | Markers found | Qualified |
|------------------------|---------------|-----------|
| Hepatocytes            | 11/50         | ✓         |
| Hepatic stellate cells | 15/50         | ✓         |
| Kupffer cells          | 15/50         | ✓         |
| LSEC                   | 13/50         | ✓         |
| Plasma cells           | 15/50         | ✓         |
| Monocytes              | 12/50         | ✓         |
| NK/NKT cells           | 11/50         | ✓         |
| Endothelial cells      | 11/50         | ✓         |
| B cells                | <5/50         | ✗ excluded |
| CD4+ T cells           | <5/50         | ✗ excluded |
| Naive T cells          | <5/50         | ✗ excluded |

**Why low overlap?** The Visium CytAssist probe panel targets 18,085 selected genes. Many immune cell surface markers (CD4, CD19, TCR genes) are not probe-targeted or are not expressed at detectable levels in bulk liver tissue sections. This is a known limitation of probe-based spatial transcriptomics for immune cell detection.

**Deconvolution finding:** All 24 spatial clusters showed "Hepatocytes" as the dominant cell type by AddModuleScore. This is biologically expected: in liver parenchyma, hepatocytes dominate both by number and transcriptional output, so the top-expressed reference genes are all hepatocytic. This does NOT mean the 24 clusters are identical — it means the variation between them reflects different **hepatocyte programs**, not cell-type composition changes.

**Limitation of AddModuleScore vs CARD:**  
AddModuleScore uses average expression (not cell-type-specific differential expression) to build marker sets. A proper deconvolution method like CARD, RCTD, or SPOTlight would:
1. Estimate fractional composition per spot (summing to 1.0 across cell types)
2. Account for spatial autocorrelation
3. Provide statistically-grounded proportions rather than relative activity scores

**Next step if CARD is needed:** The CARD compilation failure is a macOS clang linker issue, likely fixable by using a gcc-compiled R or ensuring the MCMCpack/armadillo libraries link correctly. Alternatively, RCTD (from the robust2seDE package) can be installed from Bioconductor and provides similar functionality.

---

## Key Outputs

```
results/
├── qc_plots/
│   ├── qc_summary.csv                       — per-array spot counts and QC metrics
│   └── VLP*_qc_violin.png                   — 8 per-array QC violin plots
├── spatial_domains/
│   ├── umap_clusters_and_arrays.png         — UMAP colored by cluster + by array
│   ├── elbow_plot.png                       — PCA variance explained
│   ├── cluster_sizes.csv                    — spot count per cluster
│   └── cluster_by_array.csv                 — cluster × array cross-table
├── deconvolution/
│   ├── spot_celltype_scores.csv             — 17 cell-type scores per spot (13,239 rows)
│   ├── mean_celltype_scores_per_cluster.csv — mean score per spatial cluster
│   ├── mean_celltype_scores_per_array.csv   — mean score per Visium array
│   ├── cluster_dominant_celltype.csv        — dominant cell type per cluster
│   ├── celltype_scores_heatmap_by_cluster.png — Z-scored heatmap (clusters × cell types)
│   └── celltype_scores_heatmap_by_array.png   — Z-scored heatmap (arrays × cell types)
└── masld_spatial_combined.rds              — final Seurat object (13,239 spots, 24 clusters, + scores)
```

---

## Interpretation

### What we found
1. **13,239 high-quality spots** across 8 arrays; VLP115_A is lower-quality but retained
2. **24 transcriptionally distinct spatial domains** in MASLD liver — these represent different **hepatocyte states** (metabolic zones, fibrosis programs) rather than cell-type boundaries, consistent with Vu et al.'s own findings
3. **Hepatocyte-dominated signal** throughout all clusters, as expected for liver parenchyma; immune and stromal cell contributions are present but secondary
4. **Low probe panel coverage for immune markers** limits our ability to detect T-cell, B-cell, and NK-cell spatial enrichment — this is a platform-level constraint of Visium CytAssist targeting ~18k genes

### What's missing for full biological interpretation
- **Barcode-to-patient mapping:** We don't have the `VLP*_samples.csv` files that map each spot barcode to a specific patient/biopsy. Without this, we cannot overlay fibrosis stage (F0–F4) onto individual spots. This mapping file was not included in the public UQeSpace release and is on the HPC at UQ.
- **Pathologist annotations:** Per-spot tissue compartment annotations (portal zone, lobular, fibrotic septa) are similarly unavailable from the public release.
- **True deconvolution proportions:** CARD/RCTD-based compositional estimates per spot would be more biologically informative than AddModuleScore activity scores.

### Why this still matters
The 24 clusters and their gene expression patterns can be biologically annotated using:
- Marker gene lists (FindMarkers per cluster)
- Pathway enrichment of cluster-defining genes
- Spatial visualization of known MASLD marker genes (TREM2, SPP1, ACTA2, COL1A1, PCK1, HMGCR) across the tissue sections

These analyses are what Vu et al. built their paper's conclusions on — the spatial co-localization of metabolic reprogramming and cellular senescence signatures in fibrotic tissue.

---

## Software and Methods

**R packages used:**
- Seurat 5.5.1 (spatial data loading, SCTransform, clustering)
- sctransform 0.4.3 (normalization model)
- hdf5r 1.3.12 (reading .h5 expression matrices)
- ggplot2 4.0.3, patchwork 1.3.2 (plotting)

**Deconvolution reference:**
- GSE136103 (Ramachandran et al. 2019, Nat Med) — human cirrhotic liver scRNA-seq
- 20 annotated cell types (cell_type column in nafld_seurat_annotated.rds)

**Script:** `scripts/spatial_masld_pipeline.R` (Steps 1–4 integrated) and `scripts/step4_deconvolution.R` (enhanced Step 4)

---

## Step 5 — Biological Annotation of Spatial Clusters

**Scripts:** `scripts/step5_biological_annotation.R`, `scripts/step5b_umap_labels.R`

**Approach:** `FindAllMarkers` (Seurat Wilcoxon, only.pos=TRUE, min.pct=0.1, logFC≥0.25) on the SCT assay after `PrepSCTFindMarkers`. Top marker genes per cluster were scored against 10 biological gene sets (periportal, pericentral, midzonal hepatocytes; stellate/fibrotic; Kupffer/macrophage; LSEC; biliary; lipogenic; inflammatory; senescent). The 4 genes previously validated in this project's bulk RNA-seq and scRNA-seq analyses (TREM2, SPP1, COL1A1, FASN) were explicitly tracked across all cluster markers.

**Result:** 27,416 significant marker genes (adj. p < 0.05) across 24 clusters.

### Cluster Biological Annotation

| Cluster | Biological label | Top 5 markers | Key gene hits | Spots | % |
|---------|-----------------|---------------|---------------|-------|---|
| C0  | Hepatocyte parenchyma (zone indeterminate) | MT-ND1, ZNF549, MT-ND5, ZNF550, GSTA2 | — | 2,521 | 19.0% |
| C1  | Hepatocyte parenchyma — portal/stromal features | MYH11, CCL19, LMOD1, CFTR, NOTCH3 | SPP1, COL1A1 | 1,752 | 13.2% |
| C2  | Hepatocyte parenchyma (zone indeterminate) | ASCL1, TFF3, UPP2, RET, PLCH2 | FASN | 1,061 | 8.0% |
| C3  | Inflammatory/Kupffer-macrophage region | EEF1A2, SPATA21, CXCL10, UBD, FABP4 | COL1A1, SPP1 | 1,006 | 7.6% |
| C4  | Hepatocyte parenchyma (zone indeterminate) | HYDIN, ZNF385D, DHRS2, SLC45A2, ACSL4 | — | 833 | 6.3% |
| C5  | Hepatocyte parenchyma (zone indeterminate) | APOA4, ACSL4, LOXL4, AKR1B10, HKDC1 | COL1A1, SPP1, FASN | 632 | 4.8% |
| C6  | Hepatocyte parenchyma (zone indeterminate) | SLCO1A2, PTH2R, CYP2B6, SULT1C2, CYP4F3 | FASN | 566 | 4.3% |
| C7  | Hepatocyte parenchyma (zone indeterminate) | SYT7, IFI27, HSPA5, MANF, OLFM2 | FASN | 452 | 3.4% |
| C8  | Hepatocyte parenchyma (zone indeterminate) | NUDT8, IGFBP1, GALK1, MT-ND4L, MAP4 | FASN | 451 | 3.4% |
| C9  | Hepatocyte parenchyma (zone indeterminate) | COL7A1, ESPL1, NPIPB15, CYP4X1, MROH7 | — | 407 | 3.1% |
| C10 | Pericentral hepatocyte zone | CYP3A4, IL1RAP, SLC16A1, CFHR4, EEF1B2 | — | 399 | 3.0% |
| C11 | Hepatocyte parenchyma (zone indeterminate) | PRSS51, LINGO4, NECAB2, HAMP, EYA4 | FASN | 396 | 3.0% |
| C12 | Hepatocyte parenchyma (zone indeterminate) | IGFBP1, MX1, HMCN2, MROH7, CMPK2 | — | 382 | 2.9% |
| C13 | Hepatocyte parenchyma (zone indeterminate) | SPINK1, PZP, SMIM24, HAL, BCO2 | — | 361 | 2.7% |
| C14 | Hepatocyte parenchyma (zone indeterminate) | EPB41L1, GABRB3, SLC5A12, A2M, CPN2 | — | 314 | 2.4% |
| C15 | Hepatocyte parenchyma (zone indeterminate) | TAT, MDN1, ZBTB16, GNLY, A2M | — | 293 | 2.2% |
| C16 | Hepatocyte parenchyma — interferon-stimulated | IFI6, AKR1C2, MX1, IFI44L, OAS1 | — | 286 | 2.2% |
| C17 | Pericentral / Midzonal hepatocyte zone | NUP155, CYP51A1, PLCH2, MSMO1, ADFP | FASN | 275 | 2.1% |
| C18 | Hepatocyte parenchyma (zone indeterminate) | UPK3B, FOXN4, AOC1, MEP1B, NPIPB15 | FASN | 212 | 1.6% |
| C19 | Hepatocyte parenchyma (zone indeterminate) | SLC44A5, XPNPEP2, IGFBP1, MYOM1, GPC6 | — | 191 | 1.4% |
| C20 | Hepatocyte parenchyma (zone indeterminate) | SLC29A4, DNM1, MT1H, SCAMP5, MT1G | — | 189 | 1.4% |
| C21 | Hepatocyte parenchyma (zone indeterminate) | UNC93A, ENPP3, SDS, ATAD3C, MAB21L1 | — | 117 | 0.9% |
| C22 | Pericentral / Midzonal hepatocyte zone | SPINK1, PLA2G2A, HMGCS1, LGALS4, PIK3C2G | FASN | 78 | 0.6% |
| C23 | Midzonal hepatocyte zone | GSTM1, FAM151A, ACAT2, CYP2B6, GSTM3 | FASN | 65 | 0.5% |

### Annotation rationale for key clusters

**C3 — Inflammatory/Kupffer-macrophage region:** Strongest non-hepatocyte signal in the dataset. Top markers include CXCL10 (interferon-γ-induced chemokine), UBD (ubiquitin D, interferon-stimulated), FABP4 (fatty acid binding protein 4, canonical macrophage marker), CAPG (macrophage-capping protein), GPNMB (glycoprotein NMB, enriched in activated Kupffer cells in MASLD), and EMILIN2 (extracellular matrix). COL1A1 and SPP1 are also significant markers, consistent with an inflammatory fibrotic niche where activated macrophages and myofibroblasts co-localise.

**C1 — Portal/stromal features:** Scored as "indeterminate" by the automated annotation (top markers don't match canonical hepatocyte sets) but MYH11 (smooth muscle myosin heavy chain), LMOD1 (leiomodin, smooth muscle), CCL19, CCL21 (lymphoid chemokines expressed in portal fibroblasts), NOTCH3 (portal vascular smooth muscle), and FBLN1 (fibulin, ECM) together suggest a **portal tract / smooth muscle / stromal zone**. The presence of SPP1 and COL1A1 among its significant markers confirms fibrosis-associated activity. C1 is the second-largest cluster (1,752 spots, 13.2%) and almost certainly includes portal tract regions across all biopsies.

**C10 — Pericentral hepatocyte zone:** CYP3A4 is the definitive pericentral (zone 3) marker — it is tightly restricted to centrilobular hepatocytes by the Wnt gradient. SLC16A1 (monocarboxylate transporter, zone 3 enriched) further supports this. C10 is the clearest zonation cluster in the dataset.

**C17, C22 — Pericentral/midzonal lipid synthesis:** Both clusters are enriched for cholesterol biosynthesis pathway genes: CYP51A1 (sterol demethylase), MSMO1 (methylsterol monooxygenase), HMGCS1 (HMG-CoA synthase), ACAT2 (acetyl-CoA acetyltransferase). FASN is a significant marker in both. These clusters represent **metabolically active hepatocytes engaged in de novo lipogenesis and sterol synthesis** — directly relevant to MASLD pathophysiology.

**C16 — Interferon-stimulated hepatocytes:** IFI6, MX1, IFI44L, and OAS1 are all classic interferon-stimulated genes (ISGs). This cluster (286 spots, 2.2%) likely represents hepatocytes responding to an active innate immune signal — possibly viral or cytokine-driven. This is a biologically distinct state not captured by the deconvolution but clearly identifiable through marker genes.

**C0 — High mitochondrial gene expression:** The top markers are mitochondrial-encoded genes (MT-ND1, MT-ND5, MT-ND2, MT-ATP6). This could reflect hepatocytes with genuinely elevated mitochondrial transcriptional activity, or it may indicate spots with residual red blood cell contamination (HBA2 is also in the top markers). This cluster should be interpreted cautiously.

**C20 — Metal stress response:** Top markers include MT1H and MT1G (metallothionein 1H and 1G), which are zinc/copper stress-response proteins induced in hepatocytes under oxidative stress or heavy metal exposure. May reflect a subpopulation of hepatocytes under proteotoxic stress.

### Cross-modality confirmation (TREM2 / SPP1 / COL1A1 / FASN)

| Gene | Found in spatial? | Clusters | Interpretation |
|------|-------------------|----------|----------------|
| **TREM2** | No | — | TREM2 was not identified as a significant marker in any of the 24 clusters. This is consistent with the probe panel limitation noted in Step 4: TREM2 is primarily expressed in activated macrophages/Kupffer cells, which represent a small fraction of cells in each Visium spot. It may be below the min.pct threshold even in the macrophage-enriched C3. |
| **SPP1** | Yes | C1, C3, C5 | SPP1 (osteopontin) is a significant marker in 3 clusters. In week 6 scRNA-seq analysis (GSE136103), SPP1 was the top marker for the TREM2+ scar-associated macrophage population. Here, its spatial enrichment in C1 (portal/stromal), C3 (Kupffer/inflammatory), and C5 (hepatocyte parenchyma with stromal features) confirms its role marking fibrotic niches — a direct cross-modality confirmation. |
| **COL1A1** | Yes | C1, C3, C5 | COL1A1 co-localises exactly with SPP1 across the same 3 clusters. In week 4 bulk RNA-seq, COL1A1 was significantly upregulated in fibrotic NAFLD. The spatial data confirms that COL1A1-expressing regions cluster with inflammatory macrophage and portal/stromal transcriptional states — not diffusely throughout hepatocyte parenchyma. |
| **FASN** | Yes | C2, C5, C6, C7, C8, C11, C17, C18, C22, C23 (10/24 clusters) | FASN (fatty acid synthase) is a marker in 10 of 24 clusters, spanning 2,239 spots (16.9% of tissue). This widespread spatial distribution of FASN expression confirms the de novo lipogenesis signature seen in bulk RNA-seq (week 4) as a genuinely pervasive hepatocyte-level phenomenon rather than an artifact of bulk averaging. FASN is particularly enriched in C17 and C22 (cholesterol synthesis clusters) alongside CYP51A1 and HMGCS1, pointing to a lipogenic-plus-cholesterogenic hepatocyte subpopulation. |

**Summary:** SPP1 and COL1A1 show concordant spatial co-localisation in the same 3 fibrosis/inflammation-enriched clusters — a genuine cross-modality confirmation of their role as fibrosis markers. FASN's spatial spread across 10 clusters underscores that de novo lipogenesis is a near-universal hepatocyte state in MASLD tissue, not confined to a specialist lipogenic zone. TREM2's absence as a spatial marker reflects probe panel sensitivity limits, not its biological importance in MASLD.

### Why so many "indeterminate" clusters?

16 of 24 clusters are labelled "Hepatocyte parenchyma (zone indeterminate)." This is honest, not a failure. Visium spots (~55 µm diameter, ~3–10 cells) in hepatocyte-dense liver tissue will mostly produce hepatocyte-dominated transcriptomes. Within this hepatocyte signal, many clusters differ by:
- Different expression of specific metabolic enzymes (CYP isoforms, transporters)
- Stress responses (ER stress in C7: HSPA5/MANF; interferon in C16)
- Proliferative state
- Subtle differences in lobular position not captured by our current marker sets

Without the barcode-to-patient mapping (which would allow overlaying fibrosis stage F0–F4 per spot), we cannot determine whether these indeterminate clusters preferentially associate with early vs late fibrosis stages — which is the key question Vu et al. address in their paper.

---

## Updated Key Outputs

```
results/
├── spatial_domains/
│   ├── cluster_markers_all.csv          — 27,416 marker genes (all clusters, adj.p<0.05)
│   ├── cluster_markers_top10.csv        — top 10 markers per cluster by avg_log2FC
│   ├── cluster_biological_annotation.csv — biological label + top5 + key gene hits
│   ├── cluster_annotation_summary.csv   — formatted summary table
│   ├── umap_biological_labels.png       — dual UMAP: numbered clusters | biological labels
│   └── umap_bio_labels_only.png         — biological labels only
│   (existing: cluster_sizes.csv, elbow_plot.png, umap_clusters_and_arrays.png)
└── masld_spatial_combined.rds           — Seurat object with bio_label metadata added
```

---

## Software and Methods (updated)

**R packages used:**
- Seurat 5.5.1 (spatial data loading, SCTransform, clustering, FindAllMarkers)
- sctransform 0.4.3 (normalization model)
- hdf5r 1.3.12 (reading .h5 expression matrices)
- ggplot2 4.0.3, patchwork 1.3.2 (plotting)

**Deconvolution reference:**
- GSE136103 (Ramachandran et al. 2019, Nat Med) — human cirrhotic liver scRNA-seq
- 20 annotated cell types (cell_type column in nafld_seurat_annotated.rds)

**Scripts:**
- `scripts/spatial_masld_pipeline.R` — Steps 1–4 (QC, normalization, clustering, deconvolution)
- `scripts/step4_deconvolution.R` — enhanced Step 4 (AddModuleScore cell-type scoring)
- `scripts/step5_biological_annotation.R` — Step 5 (FindAllMarkers + biological annotation)
- `scripts/step5b_umap_labels.R` — Step 5b (UMAP regeneration with biological labels)

---

## Next Steps (Week 8 Day 2+)

1. **Known MASLD gene visualization** — SpatialFeaturePlot for TREM2, SPP1, ACTA2, COL1A1, CYP2E1 to visualise spatial distribution directly on tissue
2. **COMPASS/METAFlux** — metabolic flux modeling per spatial cluster (separate task as instructed)
3. **RCTD deconvolution** (if CARD remains uninstallable) — robust cell-type decomposition from Bioconductor
4. **Contact UQ authors** for barcode-to-patient mapping files to enable per-spot fibrosis staging
