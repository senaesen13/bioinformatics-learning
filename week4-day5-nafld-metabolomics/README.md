# Week 4 Day 5 — NAFLD Plasma Metabolomics (MDPI/PMC9312563)

This folder contains plasma metabolomics data from the study:  
**"Plasma Metabolomics and Machine Learning-Driven Novel Diagnostic Signature for Non-Alcoholic Steatohepatitis"** (Biomedicines 2022, PMID: [35884973](https://pubmed.ncbi.nlm.nih.gov/35884973/), PMCID: [PMC9312563](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9312563/)).

---

## 📁 Files Included

1. **`biomedicines-10-01669.pdf`** — Full text of the published article.
2. **`biomedicines-10-01669-s001/biomedicines-1795196-supplementary.pdf`** — Official supplementary material containing figures and Table S1.
3. **`Table_S1_Metabolite_Levels.csv`** — Pure R-parsed structured table containing all 79 detected metabolites, their mean concentrations (± SD) in Healthy Controls, NAFL, and NASH groups, fold-changes, and statistical significance values (P-values and FDR-corrected Q-values).

---

## 📊 Dataset Overview (Table S1 Structure)

The parsed CSV contains the following columns:
* `No`: Original index.
* `Category`: Class of metabolite (`Amino acids`, `Kynurenine pathway metabolites and nucleosides`, `Organic acids`, `Fatty acids`).
* `Metabolite`: Common name of the metabolite (e.g. `Glutamic acid`, `Isocitric acid`, `Palmitoleic acid`).
* `Control_Mean` / `Control_SD`: Plasma concentration (ng/μL) in healthy controls.
* `NAFL_Mean` / `NAFL_SD`: Plasma concentration (ng/μL) in non-alcoholic fatty liver patients.
* `NASH_Mean` / `NASH_SD`: Plasma concentration (ng/μL) in non-alcoholic steatohepatitis patients.
* `Normalized_NAFL` / `Normalized_NASH`: Fold change relative to the control group.
* `Pval_Control_vs_NAFL` / `Pval_Control_vs_NASH` / `Pval_NAFL_vs_NASH`: Wilcoxon rank-sum test p-values.
* `Pval_Kruskal_Wallis`: Non-parametric multi-group comparison p-value.
* `Qval_FDR`: Benjamini-Hochberg false discovery rate adjusted q-value.

---

## 💡 Key Biological Signatures to Cross-Reference

The authors formulated the **"MetaNASH Score"** to diagnose NASH using three key plasma biomarkers that change dynamically:
1. **Glutamic acid** (Gradually increases: Control `6.52` → NAFL `12.74` → NASH `15.88` ng/μL, q < 0.001).
2. **Isocitric acid** (Significantly elevated in NASH vs NAFL: `0.39` vs `0.31` ng/μL, p = 0.011).
3. **Aspartic acid** (Upregulated in NASH: Control `18.57` → NASH `43.18` ng/μL, p = 0.007).

### 🧬 Connecting to your Bulk/Single-Cell Transcriptomics Pipelines

To validate these findings mechanistically, cross-reference the expression of the enzymes and transporters responsible for these metabolite shifts in your **bulk RNA-seq results** (Days 1–3) and **scRNA-seq datasets** (Week 6):

* **Glutamate / Aspartate Metabolism:**
  * Check expression of *GLS* / *GLS2* (Glutaminase), *GLUL* (Glutamine synthetase), and *GOT1* / *GOT2* (Aspartate aminotransferase — converts glutamate + oxaloacetate to aspartate + alpha-ketoglutarate).
  * Look up glutamate/aspartate transporters: *SLC1A1*, *SLC1A3*, *SLC25A12*, *SLC25A13*.
* **TCA Cycle & Organic Acids (Isocitrate, alpha-Ketoglutarate, Malate):**
  * Check *IDH1* / *IDH2* (Isocitrate dehydrogenase — converts isocitrate to alpha-ketoglutarate).
  * Check *MDH1* / *MDH2* (Malate dehydrogenase) and *ME1* / *ME2* (Malic enzyme).
* **Desaturation (Myristoleic, Palmitoleic Acid):**
  * Check *SCD* (Stearoyl-CoA desaturase — responsible for converting palmitate to palmitoleic acid). Palmitoleic acid is significantly upregulated in NAFL (1.69x) and NASH (2.20x).
