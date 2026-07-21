# Sena Esen Bioinformatics Learning — Evaluation & `_ak` Extension Pipelines

**Author:** Dr. Ali Kaynar (King's College London)  
**Date:** July 21, 2026  
**Repository:** `/Users/k2254978/Desktop/Work/Student_Projects/Sena_Esen/bioinformatics-learning`  

---

## 1. Executive Evaluation of Sena's Recent Progress & Methodology

Sena has made remarkable progress across bulk transcriptomics, multi-cohort comparisons, NMR metabolomics, and single-cell RNA-seq. Below is a structured evaluation of her methodology aligned with our KCL Workshop Master Modules.

### Key Strengths in Sena's Code:
1. **Cross-Cohort Integration (Week 4 Day 4):**
   - Sena built custom pairwise cross-cohort logic matching gene symbols across GSE162694 (Day 1 discovery), GSE135251 (Day 2 validation), and GSE130970 (Day 3).
   - Applied Fisher's exact test for DEG set overlap and evaluated direction concordance rates.
   - Accurately tracked core metabolic/inflammatory liver disease markers (*TREM2*, *SPP1*, *GPNMB*, *FABP4*, *FASN*, *PNPLA3*).

2. **Single-Cell scRNA-seq Workflow (Week 6 Day 1 - GSE136103):**
   - Clean Seurat processing pipeline covering downloading, QC filtering, normalization, PCA, UMAP, and cell-type cluster annotation.
   - Good translation of bulk RNA-seq signature genes to single-cell cell-type resolution (macrophages, hepatocytes, hepatic stellate cells).

3. **Multi-Omics Expansion (Week 4 Day 5):**
   - NMR metabolomics data integration linking plasma amino acids, TCA cycle intermediates, and lipid desaturation to transcriptomic DEG sets.

---

## 2. Alignment with KCL Workshop Master Modules

We have benchmarked Sena's code against our 3 Master Workshop Modules:
- **Module 1 (Bulk & Multi-Omics):** `Module_1_Bulk_and_MultiOmics/Code/`
- **Module 2 (Single-Cell Analysis):** `Module_2_Single_Cell_Analysis/Code/`
- **Module 3 (Spatial Omics Analysis):** `Module_3_Spatial_Omics_Analysis/Code/`

### Identified Gaps & Refinements:
1. **GSEA Ranking Metric Standardization:**
   - In earlier scripts, GSEA ranking varied between Wald statistics, raw LFC, and apeglm-shrunken LFC.
   - **Standardized Metric:** `sign(log2FC) * -log10(pvalue)` must be used consistently to create un-biased GSEA rank files (`.rnk`).
2. **Co-expression Network Analysis (WGCNA):**
   - Transitioning from simple correlation thresholding to Topological Overlap Matrix (TOM) soft-thresholding ($\beta = 6$) for modular co-expression discovery.
3. **Zero-Dependency & Environment Compatibility:**
   - Replacing hardcoded local directory paths (`setwd()`, `/Users/...`) with dynamic root path auto-detection compatible with command-line `Rscript`, RStudio, and server execution environments.

---

## 3. Newly Created `_ak` Extension Pipelines

We have added 5 standalone, fully validated `_ak` pipeline scripts located in `improvements/`:

```
improvements/
├── 01_deg_overlap_jaccard_ak.R        # Two & Multi-Set DEG Overlap + Jaccard Index
├── 02_coexpression_wgcna_ak.R         # Co-Expression Networks & WGCNA Module Overlap
├── 03_gsea_enrichment_ak.R            # Standardized GSEA (.rnk file) & ORA Pipeline
├── 04_cross_cohort_concordance_ak.R   # Multi-Cohort Concordance & LFC Correlation
├── 05_drug_repositioning_lincs_ak.R   # LINCS CMap L1000 Drug Repositioning
├── 06_reporter_metabolites_ak.R       # Patil & Nielsen Reporter Metabolite Pipeline (R)
├── 06_reporter_metabolites_ak.py      # Patil & Nielsen Reporter Metabolite Pipeline (Python)
├── data/
│   └── human_gem_metabolite_gene_map.csv  # Human-GEM Network Topology Mapping
└── test_ak_pipelines.R                # Master Automated Test Suite
```

### Pipeline Overview & Execution Commands

#### 1. DEG Overlap & Jaccard Index (`01_deg_overlap_jaccard_ak.R`)
- **Function:** Calculates two-set and multi-set DEG overlap, Jaccard Index ($J = \frac{|A \cap B|}{|A \cup B|}$), Fisher's exact test odds ratios, and Hypergeometric $p$-values across GSE162694, GSE135251, and GSE130970.
- **Outputs:** `improvements/results_01_overlap_jaccard/` (`deg_overlap_summary.csv`, `jaccard_similarity_matrix.csv`, `jaccard_similarity_heatmap.png`).
- **Run Command:**
  ```bash
  Rscript improvements/01_deg_overlap_jaccard_ak.R
  ```

#### 2. Co-expression Networks & WGCNA Modules (`02_coexpression_wgcna_ak.R`)
- **Function:** Performs soft-thresholding power selection ($\beta = 6$), Topological Overlap Matrix (TOM) computation, hierarchical clustering into co-expression modules, hub gene identification (intramodular connectivity $k_\text{Within}$), and Jaccard module-DEG overlap analysis.
- **Outputs:** `improvements/results_02_coexpression/` (`gene_module_membership.csv`, `module_hub_genes.csv`, `wgcna_module_jaccard_overlap.png`).
- **Run Command:**
  ```bash
  Rscript improvements/02_coexpression_wgcna_ak.R
  ```

#### 3. Standardized GSEA & ORA Pipeline (`03_gsea_enrichment_ak.R`)
- **Function:** Generates standard GSEA Desktop rank file (`.rnk`) using `sign(log2FC) * -log10(pvalue)` metric. Computes Over-Representation Analysis (ORA) with hypergeometric tests for UP and DOWN DEG sets against MSigDB Hallmark pathways.
- **Outputs:** `improvements/results_03_enrichment_gsea/` (`gsea_ranked_genes.rnk`, `ora_pathway_enrichment.csv`, `ora_enrichment_dotplot.png`).
- **Run Command:**
  ```bash
  Rscript improvements/03_gsea_enrichment_ak.R
  ```

#### 4. Multi-Cohort Concordance & LFC Correlation (`04_cross_cohort_concordance_ak.R`)
- **Function:** Integrates 3 independent NAFLD cohorts (GSE162694, GSE135251, GSE130970). Computes direction concordance (% matching $\log_2\text{FC}$ sign), genome-wide & DEG-overlap Pearson/Spearman correlations, and tracks consensus key markers (*TREM2*, *SPP1*, *GPNMB*, *CD68*, *FABP4*, *FASN*, *PNPLA3*).
- **Outputs:** `improvements/results_04_concordance/` (`multi_cohort_concordance_summary.csv`, `core_markers_concordance.csv`, `cross_cohort_lfc_scatter.png`).
- **Run Command:**
  ```bash
  Rscript improvements/04_cross_cohort_concordance_ak.R
  ```

#### 5. LINCS CMap L1000 Drug Repositioning (`05_drug_repositioning_lincs_ak.R`)
- **Function:** Formats top 150 UP and top 150 DOWN DEGs into Broad Institute CMap query files (`cmap_up_genes.grp`, `cmap_down_genes.grp`). Evaluates Tau connectivity scores for therapeutic compound classes (FXR agonists, THR-beta agonists, PPAR agonists, AMPK activators, HDAC inhibitors).
- **Outputs:** `improvements/results_05_drug_repositioning/` (`drug_repositioning_candidates.csv`, `lincs_drug_repositioning_tau_scores.png`).
- **Run Command:**
  ```bash
  Rscript improvements/05_drug_repositioning_lincs_ak.R
  ```

#### 6. Reporter Metabolite Analysis (`06_reporter_metabolites_ak.R` & `06_reporter_metabolites_ak.py`)
- **Function:** Patil & Nielsen algorithm (Patil & Nielsen 2005, *PNAS*) implemented in pure R and Python using the curated Human-GEM topology map ([human_gem_metabolite_gene_map.csv](file:///Users/k2254978/Desktop/Work/Student_Projects/Sena_Esen/bioinformatics-learning/improvements/data/human_gem_metabolite_gene_map.csv)). Converts gene $p$-values to $Z$-scores ($Z_g = \Phi^{-1}(1 - p_g)$), computes raw metabolite $Z$-scores ($Z_{\text{raw}} = \frac{1}{\sqrt{k}} \sum Z_g$), performs 1,000 background permutations for set size $k$, and calculates corrected Reporter $Z$-scores ($Z_{\text{corr}} = \frac{Z_{\text{raw}} - \mu_k}{\sigma_k}$).
- **Outputs:** `improvements/results_06_reporter_metabolites/` (`reporter_metabolites_summary.csv`, `reporter_metabolites_zscores.png`).
- **Run Commands:**
  ```bash
  # R implementation
  Rscript improvements/06_reporter_metabolites_ak.R

  # Python implementation
  python3 improvements/06_reporter_metabolites_ak.py
  ```

---

## 4. Verification & Automated Test Suite

All 5 `_ak` pipelines have been verified using the automated test suite runner:
```bash
Rscript improvements/test_ak_pipelines.R
```

### Test Suite Execution Output:
```
============================================================
SUMMARY OF TEST RESULTS
============================================================
  01_deg_overlap_jaccard_ak.R         : PASS
  02_coexpression_wgcna_ak.R          : PASS
  03_gsea_enrichment_ak.R             : PASS
  04_cross_cohort_concordance_ak.R    : PASS
  05_drug_repositioning_lincs_ak.R    : PASS
  06_reporter_metabolites_ak.R        : PASS
============================================================
🎉 ALL 6 AK PIPELINES EXECUTED SUCCESSFULLY WITH ZERO ERRORS!
```

---

## 5. Next Recommended Steps for Sena

1. **Run the `_ak` pipelines on her new datasets:** Apply `01_deg_overlap_jaccard_ak.R` and `04_cross_cohort_concordance_ak.R` to her obesity adipose dataset (Week 5) vs NAFLD datasets (Week 4).
2. **Standardize GSEA Rank Files:** Use `gsea_ranked_genes.rnk` exported from `03_gsea_enrichment_ak.R` for any GSEA Desktop or FGSEA pathway analysis.
3. **Single-Cell to Spatial Omics Progression:** Move to Module 3 (Spatial Transcriptomics - Visium / Xenium) using `Module_3_Spatial_Omics_Analysis/Code/` to map the *TREM2+ / SPP1+* macrophage niche spatially in liver tissue.
