# NAFLD Cross-Cohort Pairwise Comparison: GSE130970 vs GSE162694 and GSE135251

Generated: 2026-07-14 | Significance threshold: padj < 0.05 AND |MLE log2FC| > 1  
Gene matching: by **gene symbol** (cohorts use different ID systems: Ensembl vs Entrez)

---

## Cohort Overview

| Cohort | Folder | Samples (Normal / NAFLD) | Control definition | Sig genes |
|---|---|---|---|---|
| GSE162694 (Day 1, Discovery) | `week4-day1-nafld-bulk-rnaseq/` | 31 / 112 | Histologically normal liver | 485 |
| GSE135251 (Day 2, Validation) | `week4-day2-nafld-validation-rnaseq/` | 10 / 206 | Disease: Control | 1058 |
| GSE130970 (Day 3, this cohort) | `week4-day3-nafld-gse130970-rnaseq/` | 4 / 74 | NAS = 0 | 180 |

---

## Comparison A: GSE130970 vs GSE162694 (Day 1 Discovery)

### Overlap

| Metric | Value |
|---|---|
| GSE130970 significant genes | 180 |
| GSE162694 significant genes | 485 |
| Shared significant genes (by symbol) | **38** |
| Overlap as % of GSE130970 | 21.1% |
| Universe (shared tested genes by symbol) | 15,036 |

### Statistical Agreement

| Metric | Value |
|---|---|
| Fisher's exact OR | **8.68** |
| Fisher's exact p-value | **8.7e-21** |
| Direction concordance | **94.7%** (36/38 concordant) |
| Pearson r (sig overlap, n=38) | **0.656** (p = 7.7e-06) |
| Spearman ρ (sig overlap, n=38) | **0.411** (p = 0.011) |
| Pearson r (genome-wide, n=15,036) | **0.461** |

### TREM2 / SPP1 / GPNMB

| Gene | GSE130970 log2FC | GSE130970 padj | GSE130970 sig | GSE162694 log2FC | GSE162694 padj | GSE162694 sig | Match |
|---|---|---|---|---|---|---|---|
| TREM2 | +1.51 | 1.40e-01 | NO | +2.48 | 5.17e-13 | YES | Concordant direction; not sig in Day 3 (n=4 controls) |
| SPP1 | +0.18 | 8.97e-01 | NO | +1.87 | 1.10e-07 | YES | Same direction; below LFC threshold in Day 3 |
| GPNMB | +0.30 | 6.05e-01 | NO | +1.19 | 2.25e-09 | YES | Same direction; below threshold in Day 3 |

### LFC Correlation Plot

See `results/lfc_corr_d3_vs_d1.png`

---

## Comparison B: GSE130970 vs GSE135251 (Day 2 Validation)

### Overlap

| Metric | Value |
|---|---|
| GSE130970 significant genes | 180 |
| GSE135251 significant genes | 1058 |
| Shared significant genes (by symbol) | **52** |
| Overlap as % of GSE130970 | 28.9% |
| Universe (shared tested genes by symbol) | 13,016 |

### Statistical Agreement

| Metric | Value |
|---|---|
| Fisher's exact OR | **4.78** |
| Fisher's exact p-value | **1.6e-16** |
| Direction concordance | **84.6%** (44/52 concordant) |
| Pearson r (sig overlap, n=52) | **0.510** (p = 1.1e-04) |
| Spearman ρ (sig overlap, n=52) | **0.353** (p = 0.011) |
| Pearson r (genome-wide, n=13,016) | **0.219** |

### TREM2 / SPP1 / GPNMB

| Gene | GSE130970 log2FC | GSE130970 padj | GSE130970 sig | GSE135251 log2FC | GSE135251 padj | GSE135251 sig | Match |
|---|---|---|---|---|---|---|---|
| TREM2 | +1.51 | 1.40e-01 | NO | +2.52 | 2.22e-07 | YES | Concordant direction; not sig in Day 3 |
| SPP1 | +0.18 | 8.97e-01 | NO | +1.41 | 2.76e-03 | YES | Concordant direction; below threshold in Day 3 |
| GPNMB | +0.30 | 6.05e-01 | NO | +0.91 | 6.00e-03 | NO | Concordant direction; below threshold in both |

### LFC Correlation Plot

See `results/lfc_corr_d3_vs_d2.png`

---

## Side-by-Side Summary

| Metric | Comparison A (D3 vs D1) | Comparison B (D3 vs D2) |
|---|---|---|
| Overlap genes | 38 | 52 |
| Overlap % of D3 | 21.1% | 28.9% |
| Fisher's OR | **8.68** | 4.78 |
| Fisher's p | 8.7e-21 | 1.6e-16 |
| Direction concordance | **94.7%** | 84.6% |
| Pearson r (sig overlap) | **0.656** | 0.510 |
| Spearman ρ (sig overlap) | **0.411** | 0.353 |
| Genome-wide Pearson r | **0.461** | 0.219 |

---

## Conclusion

**GSE130970 shows stronger agreement with GSE162694 (Day 1 Discovery, Comparison A)
than with GSE135251 (Day 2 Validation, Comparison B)**, on every metric:

- Higher Fisher's OR (8.68 vs 4.78): the overlap is more enriched above chance in Comparison A
- Higher direction concordance (94.7% vs 84.6%): 3 out of 8 discordant genes fewer in Comparison A
- Higher Pearson r on shared significant genes (0.656 vs 0.510): stronger quantitative agreement
- Higher genome-wide Pearson r (0.461 vs 0.219): more coherent transcriptome-wide LFC patterns

**Why Comparison A is stronger:**

1. **Fibrosis stage distribution match.** GSE162694 (Day 1) spans all fibrosis stages
   F0–F4 (median ~F1) with 31 controls, similar to GSE130970's F0-dominated distribution
   (F0+F1 = 53/78 samples). In contrast, GSE135251 (Day 2) is heavily enriched for
   advanced fibrosis (NASH F3/F4 = 68/216 samples), producing a transcriptional profile
   dominated by hepatocyte loss and ECM replacement that is less similar to earlier-stage
   disease.

2. **Control group quality.** GSE162694 has 31 normal controls (matched histology), while
   GSE135251 has 10 and GSE130970 has only 4. The 31-control cohort produces more robust
   reference-group variance estimates, making the effect size estimates more comparable
   to GSE130970's (also representing a mixed-stage NAFLD population).

3. **Gene ID system.** Both Day 1 (Ensembl → symbol) and Day 3 (Entrez → symbol) mappings
   introduce some noise through org.Hs.eg.db, but the Day 2 dataset has an additional
   layer of cross-cohort complexity due to GTF version differences in the two mapping chains.

**Caveat:** GSE130970's n=4 Normal group severely limits the discovery list to 180 genes —
only the most robust signals survive. This is why Comparison B (Day 3 vs Day 2) shows
more raw overlap (52 genes) despite weaker statistical agreement: the larger Day 2
gene list is more likely to overlap with any given gene set by chance. The Fisher's OR
controls for this and correctly identifies Day 1 as the stronger partner.

**Marker summary:** TREM2, SPP1, and GPNMB show the same positive direction (upregulated
in NAFLD) in GSE130970 as in the other two cohorts, but do not reach individual gene
significance due to the small control group. The concordant direction across all three
cohorts strengthens confidence in these markers as genuine NAFLD signals despite the
GSE130970 power limitation.
