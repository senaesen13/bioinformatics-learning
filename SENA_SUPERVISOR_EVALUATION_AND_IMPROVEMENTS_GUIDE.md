# Detailed Evaluation, Supervisor Improvements, and `_ak` Code Guide for Sena Esen

**Student:** Sena Esen  
**Supervisor:** Dr. Ali Kaynar (Postdoctoral Research Associate, Systems Medicine / Microbiome / Metabolic Modelling, King's College London)  
**Date:** July 21, 2026  
**Repository:** [github.com/senaesen13/bioinformatics-learning](https://github.com/senaesen13/bioinformatics-learning)  

---

## Executive Summary & Supervisor Note

Dear Sena,

You have made outstanding progress in building a multi-omics bioinformatics portfolio. Across bulk RNA-seq, multi-cohort validation, NMR metabolomics, and single-cell RNA-seq, your code demonstrates a strong understanding of computational biology and translational disease modeling.

To help elevate your methodology to publication standards and ensure seamless compatibility across research environments, I have conducted a thorough evaluation of your work, added methodological corrections, and created **6 standardized supervisor pipelines (`_ak` versions)** in your repository under `improvements/`.

All 6 pipelines have been tested with an automated test suite (`improvements/test_ak_pipelines.R`) and pass with 100% clean, zero-error execution on your real datasets.

---

## 1. Detailed Week-by-Week Evaluation & Methodological Enhancements

### Week 1 & Week 2: R Basics, Bulk RNA-seq, DESeq2, and GSEA
- **Sena's Work:** Strong implementation of basic data structures, DESeq2 differential expression, and GSEA pathway analysis.
- **Supervisor Enhancement 1 — GSEA Ranking Metric Standardization:**
  - *Problem:* In earlier scripts, GSEA ranking statistics varied between Wald statistics, raw $\log_2\text{FC}$, and apeglm-shrunken $\log_2\text{FC}$. Shrunken LFC compresses low-count genes toward zero, introducing rank bias.
  - *Fix:* Implemented a standardized ranking metric everywhere:
    $$\text{Ranking Metric} = \text{sign}(\log_2\text{FC}) \times -\log_{10}(p\text{-value})$$
  - *Implementation:* Created `improvements/03_gsea_enrichment_ak.R` which exports standard GSEA Desktop compatible `.rnk` files (`gsea_ranked_genes.rnk`).

### Week 3: Single-Cell RNA-seq, Spatial Transcriptomics & Drug Repositioning (ccRCC)
- **Sena's Work:** Good Seurat workflow on PBMC3k, 10x Visium heart spatial transcriptomics, and LINCS CMap drug repositioning.
- **Supervisor Enhancement 2 — Kaplan–Meier Optimal-Cutpoint False-Positive Correction (`km_cutpoint_corrected.R`):**
  - *Problem:* Scanning cutpoints across the 20th–80th percentiles and reporting the minimum $p$-value uncorrected inflates false-positive rates (~5.5× inflation under null hypothesis).
  - *Fix:* Created `improvements/km_cutpoint_corrected.R` applying **Lausen & Schumacher (1992) maxstat correction** alongside pre-specified median splits.
- **Supervisor Enhancement 3 — Single-Cell Doublet Removal:**
  - *Fix:* Added `scDblFinder` (scverse recommended) to `week3-day1-scrna-seq/scripts/seurat_pbmc3k.R` to eliminate heterotypic cell doublets.
- **Supervisor Enhancement 4 — Spatial QC "Before" Count Bug:**
  - *Fix:* Corrected the spot count tracking in `week3-day3-spatial-transcriptomics-heart/scripts/spatial_heart.R` to capture pre-QC spot counts accurately.

### Week 4: NAFLD Bulk RNA-seq, 3-State DESeq2, Cross-Cohort Comparison & Metabolomics
- **Sena's Work:** Excellent multi-cohort pairwise comparisons matching gene symbols across **GSE162694** (Day 1 discovery), **GSE135251** (Day 2 validation), and **GSE130970** (Day 3).
- **Supervisor Enhancement 5 — Multi-Cohort Cross-Cohort Concordance:**
  - *Implementation:* Built `improvements/04_cross_cohort_concordance_ak.R` calculating Fisher's exact test, direction concordance rate (% matching $\log_2\text{FC}$ sign), and tracking core consensus markers (*TREM2*, *SPP1*, *GPNMB*, *CD68*, *FABP4*, *FASN*, *PNPLA3*).
- **Supervisor Enhancement 6 — NAFLD NMR Metabolomics Data Integration:**
  - Extracted supplementary data from the MDPI NASH metabolomics paper (DOI: 10.3390/biomedicines10071669) and created `week4-day5-nafld-metabolomics/Table_S1_Metabolite_Levels.csv` linking plasma amino acids, TCA cycle intermediates, and lipid desaturation with transcriptomic DEG sets.

### Week 5: Obesity Adipose Tissue RNA-seq
- **Sena's Work:** Extended transcriptomic workflow to human visceral/subcutaneous adipose tissue (GSE152991).
- **Supervisor Enhancement 7 — Location-Independent Dynamic Path Detection:**
  - *Fix:* Replaced hardcoded local paths (`/Users/senaesen/...`) with `get_root_dir()` dynamic path resolution compatible across macOS, Linux, and Windows environments.

### Week 6: scRNA-seq GSE136103 Human Liver
- **Sena's Work:** Cell-type clustering (macrophages, hepatocytes, hepatic stellate cells) and translation of bulk NAFLD DEG signatures to cell-type specific resolution.

---

## 2. Alignment with KCL Master Workshop Modules

Your work aligns directly with our 3 Master Workshop Modules located at `/Users/k2254978/Desktop/Work/Workshops/`:

```
Workshops/
├── Module_1_Bulk_and_MultiOmics/
│   └── Code/ (01a–01b DESeq2, 02 Variant Calling, 03 Epigenomics, 04 Proteomics,
│             05 Metabolomics, 06 Metallomics, 07 DEG Overlap, 08 WGCNA,
│             09 GSEA/clusterProfiler, 10 Reporter Metabolites, 11 LINCS CMap,
│             12 Survival, 13 MOFA2, 14 Cohort Concordance, 15 Disease Staging)
├── Module_2_Single_Cell_Analysis/
│   └── Code/ (01 QC/Doublets, 02 Normalization/Integration, 03 Clustering, 04 CellChat)
└── Module_3_Spatial_Omics_Analysis/
    └── Code/ (01 Spot QC, 02 Deconvolution, 03 Spatial Domains, 04 COMPASS/METAFlux)
```

---

## 3. Comprehensive Guide to the 6 Standardized `_ak` Pipelines

All supervisor pipelines are saved in your repository under `improvements/`:

```
improvements/
├── 01_deg_overlap_jaccard_ak.R        # Two & Multi-Set DEG Overlap + Jaccard Index
├── 02_coexpression_wgcna_ak.R         # Co-Expression Networks & WGCNA Module Overlap
├── 03_gsea_enrichment_ak.R            # Standardized GSEA (.rnk file) & ORA Pipeline
├── 04_cross_cohort_concordance_ak.R   # Multi-Cohort Concordance & LFC Correlation
├── 05_drug_repositioning_lincs_ak.R   # LINCS CMap L1000 Drug Repositioning
├── 06_reporter_metabolites_ak.R       # Patil & Nielsen GEM Reporter Metabolite (R)
├── 06_reporter_metabolites_ak.py      # Patil & Nielsen GEM Reporter Metabolite (Python)
├── parse_gem.py                       # GEM SBML XML Parser
├── data/
│   ├── Human-GEM.xml                  # Full Human-GEM 2.0 SBML Model
│   ├── human_gem_topology_network.csv # Parsed 38,527 Edge GEM Network Mapping
│   └── human_gem_metabolite_gene_map.csv # Curated Intermediate Mapping
└── test_ak_pipelines.R                # Master Automated Test Suite
```

### Overview of Each `_ak` Pipeline

#### Pipeline 01: DEG Overlap & Jaccard Index (`01_deg_overlap_jaccard_ak.R`)
- **What it does:** Calculates pairwise and 3-way DEG set overlap, Jaccard Index ($J = \frac{|A \cap B|}{|A \cup B|}$), Fisher's exact test odds ratios, and Hypergeometric $p$-values across your NAFLD cohorts.
- **Run Command:**
  ```bash
  Rscript improvements/01_deg_overlap_jaccard_ak.R
  ```
- **Outputs:** `improvements/results_01_overlap_jaccard/` (`deg_overlap_summary.csv`, `jaccard_similarity_heatmap.png`).

#### Pipeline 02: Co-expression Networks & WGCNA Modules (`02_coexpression_wgcna_ak.R`)
- **What it does:** Calculates Spearman correlation, soft-thresholding power selection ($\beta = 6$), Topological Overlap Matrix (TOM), hierarchical module clustering, intramodular hub gene connectivity ($k_\text{Within}$), and module-DEG overlap analysis.
- **Run Command:**
  ```bash
  Rscript improvements/02_coexpression_wgcna_ak.R
  ```
- **Outputs:** `improvements/results_02_coexpression/` (`gene_module_membership.csv`, `module_hub_genes.csv`, `wgcna_module_jaccard_overlap.png`).

#### Pipeline 03: Standardized GSEA & ORA Pipeline (`03_gsea_enrichment_ak.R`)
- **What it does:** Uses $\text{sign}(\log_2\text{FC}) \times -\log_{10}(p\text{-value})$ to generate standard GSEA Desktop rank files (`gsea_ranked_genes.rnk`). Calculates Over-Representation Analysis (ORA) with hypergeometric tests for UP and DOWN DEGs.
- **Run Command:**
  ```bash
  Rscript improvements/03_gsea_enrichment_ak.R
  ```
- **Outputs:** `improvements/results_03_enrichment_gsea/` (`gsea_ranked_genes.rnk`, `ora_pathway_enrichment.csv`, `ora_enrichment_dotplot.png`).

#### Pipeline 04: Multi-Cohort Concordance & LFC Correlation (`04_cross_cohort_concordance_ak.R`)
- **What it does:** Compares DESeq2 results across GSE162694, GSE135251, and GSE130970. Calculates direction concordance (% matching $\log_2\text{FC}$ sign), genome-wide Pearson/Spearman correlations, and tracks consensus key markers (*TREM2*, *SPP1*, *GPNMB*, *CD68*, *FABP4*, *FASN*, *PNPLA3*).
- **Run Command:**
  ```bash
  Rscript improvements/04_cross_cohort_concordance_ak.R
  ```
- **Outputs:** `improvements/results_04_concordance/` (`multi_cohort_concordance_summary.csv`, `core_markers_concordance.csv`, `cross_cohort_lfc_scatter.png`).

#### Pipeline 05: LINCS CMap L1000 Drug Repositioning (`05_drug_repositioning_lincs_ak.R`)
- **What it does:** Formats top 150 UP and top 150 DOWN DEGs into Broad Institute CMap query files (`cmap_up_genes.grp`, `cmap_down_genes.grp`). Evaluates Tau connectivity scores for candidate drug classes (FXR agonists, THR-beta agonists, PPAR agonists, AMPK activators, HDAC inhibitors).
- **Run Command:**
  ```bash
  Rscript improvements/05_drug_repositioning_lincs_ak.R
  ```
- **Outputs:** `improvements/results_05_drug_repositioning/` (`drug_repositioning_candidates.csv`, `lincs_drug_repositioning_tau_scores.png`).

#### Pipeline 06: Genome-Scale Metabolic Model Reporter Metabolite Analysis (`06_reporter_metabolites_ak.R` & `06_reporter_metabolites_ak.py`)
- **What it does:** Patil & Nielsen algorithm (Patil & Nielsen 2005, *PNAS*) implemented in pure R and Python using the full **Human-GEM 2.0 SBML model (`Human-GEM.xml`)**. Converts gene $p$-values to $Z$-scores ($Z_g = \Phi^{-1}(1 - p_g/2)$), computes raw metabolite $Z$-scores ($Z_{\text{raw}} = \frac{1}{\sqrt{k}} \sum Z_g$), performs 1,000 Monte Carlo background permutations, and calculates corrected Reporter $Z$-scores ($Z_{\text{corr}} = \frac{Z_{\text{raw}} - \mu_k}{\sigma_k}$).
- **Run Commands:**
  ```bash
  # R implementation
  Rscript improvements/06_reporter_metabolites_ak.R

  # Python implementation
  python3 improvements/06_reporter_metabolites_ak.py --gem improvements/data/Human-GEM.xml
  ```
- **Outputs:** `improvements/results_06_reporter_metabolites/` (`gem_reporter_metabolites_summary.csv`, `gem_reporter_metabolites_zscores.png`).

---

## 4. Automated Verification & Testing

You can run the entire test suite at any time from your repository root:

```bash
Rscript improvements/test_ak_pipelines.R
```

### Expected Output:
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

## 5. Next Action Items for Your Research Project

1. **Cross-Disease Comparison:** Apply `01_deg_overlap_jaccard_ak.R` and `04_cross_cohort_concordance_ak.R` to compare your Obesity Adipose dataset (Week 5) against your NAFLD datasets (Week 4).
2. **GSEA Desktop Enrichment:** Use `improvements/results_03_enrichment_gsea/gsea_ranked_genes.rnk` for GSEA Desktop or FGSEA pathway analysis.
3. **Spatial Omics Progression:** Advance to **Module 3 (Spatial Omics Analysis)** in `Workshops/Module_3_Spatial_Omics_Analysis/Code/` to map the *TREM2+ / SPP1+* macrophage niche spatially in liver tissue using Visium and Xenium data.

Keep up the fantastic work!

Best regards,  
**Dr. Ali Kaynar**  
*King's College London*
