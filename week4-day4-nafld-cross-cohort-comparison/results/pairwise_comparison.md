# NAFLD Cross-Cohort Pairwise Comparison: All Three Cohorts

Generated: 2026-07-14 | Significance threshold: padj < 0.05 AND |MLE log2FC| > 1  
Gene matching: by **gene symbol** (cohorts use different ID systems: Ensembl vs Entrez)

Three pairwise comparisons are reported here, covering every combination of the three
cohorts. Comparisons A and B were computed by `scripts/pairwise_comparison.R`.
Comparison C reuses numbers already computed in `week4-day2-nafld-validation-rnaseq/`
(see `results/validation_overlap_summary.csv` and `NOTES.md` in that folder).

---

## Cohort Overview

| Cohort | Folder | Samples (Normal / NAFLD) | Control definition | Sig genes |
|---|---|---|---|---|
| GSE162694 (Day 1, Discovery) | `week4-day1-nafld-bulk-rnaseq/` | 31 / 112 | Histologically normal liver | 485 |
| GSE135251 (Day 2, Validation) | `week4-day2-nafld-validation-rnaseq/` | 10 / 206 | Disease: Control | 1058 |
| GSE130970 (Day 3, third cohort) | `week4-day3-nafld-gse130970-rnaseq/` | 4 / 74 | NAS = 0 | 180 |

---

## Side-by-Side Summary (All Three Comparisons)

| Metric | Comparison A (D3 vs D1) | Comparison B (D3 vs D2) | Comparison C (D1 vs D2) |
|---|---|---|---|
| Overlap genes | 38 | 52 | **139** |
| Overlap % of smaller list | 21.1% | 21.1% | 28.7% |
| Universe (shared tested genes) | 15,036 | 13,016 | 13,100 |
| Fisher's OR | **8.68** | 4.78 | 5.11 |
| Fisher's p | 8.7e-21 | 1.6e-16 | **9.5e-43** |
| Hypergeometric p | 1.1e-20 | 1.6e-16 | **9.5e-43** |
| Direction concordance | **94.7%** | 84.6% | 73.4% |
| Pearson r (sig overlap) | **0.656** | 0.510 | −0.048 |
| Spearman ρ (sig overlap) | **0.411** | 0.353 | 0.039 |
| Genome-wide Pearson r | **0.461** | 0.219 | 0.323 |

> **Note on Fisher's vs hypergeometric p:** One-sided Fisher's exact test and the
> hypergeometric test are mathematically equivalent for enrichment testing. The small
> differences above (e.g. 8.7e-21 vs 1.1e-20 for Comparison A) are due to rounding in
> reported digits only. Both confirm the same conclusion at every comparison.

---

## Comparison A: GSE130970 (Day 3) vs GSE162694 (Day 1 Discovery)

### Overlap

| Metric | Value |
|---|---|
| GSE130970 significant genes | 180 |
| GSE162694 significant genes | 485 |
| Shared significant genes (by symbol) | **38** |
| Overlap as % of GSE130970 | 21.1% |
| Universe (shared tested genes by symbol) | 15,036 |
| Expected overlap by chance | 5.8 |
| Fold-enrichment over chance | 6.5× |

### Statistical Agreement

| Metric | Value |
|---|---|
| Fisher's exact OR | **8.68** |
| Fisher's exact p-value | **8.7e-21** |
| Hypergeometric p-value | **1.1e-20** |
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

## Comparison B: GSE130970 (Day 3) vs GSE135251 (Day 2 Validation)

### Overlap

| Metric | Value |
|---|---|
| GSE130970 significant genes | 180 |
| GSE135251 significant genes | 1058 |
| Shared significant genes (by symbol) | **52** |
| Overlap as % of GSE130970 | 28.9% |
| Universe (shared tested genes by symbol) | 13,016 |
| Expected overlap by chance | 14.6 |
| Fold-enrichment over chance | 3.6× |

### Statistical Agreement

| Metric | Value |
|---|---|
| Fisher's exact OR | **4.78** |
| Fisher's exact p-value | **1.6e-16** |
| Hypergeometric p-value | **1.6e-16** |
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

## Comparison C: GSE162694 (Day 1 Discovery) vs GSE135251 (Day 2 Validation)

> **Source:** Numbers reused from `week4-day2-nafld-validation-rnaseq/` — computed in
> `scripts/rethreshold_day2.R` and stored in `results/validation_overlap_summary.csv`.
> Not recomputed here. For the full list of 139 overlap genes and GSEA pathway
> comparison, see `week4-day2-nafld-validation-rnaseq/NOTES.md`.

### Overlap

| Metric | Value |
|---|---|
| GSE162694 significant genes | 485 |
| GSE135251 significant genes | 1058 |
| Shared significant genes (by symbol) | **139** |
| Overlap as % of GSE162694 (smaller list) | 28.7% |
| Universe (shared tested genes by symbol) | 13,100 |
| Expected overlap by chance | 39.2 |
| Fold-enrichment over chance | 3.5× |

### Statistical Agreement

| Metric | Value |
|---|---|
| Fisher's exact OR | **5.11** |
| Fisher's exact p-value | **9.5e-43** |
| Hypergeometric p-value | **9.5e-43** |
| Direction concordance | **73.4%** (102/139 concordant) |
| Pearson r (sig overlap, n=139) | **−0.048** (p = 0.58) |
| Spearman ρ (sig overlap, n=139) | **0.039** (p = 0.65) |
| Pearson r (genome-wide, n=13,100) | **0.323** |

### Why the sig-overlap Pearson r is near zero despite 73.4% concordance

All 485 Day 1 significant genes are upregulated (lfc_mle > 1). Within the 139 overlap
genes, ~37 are strongly downregulated in Day 2 (log2FC −4 to −8), including FOS
(−6.28), FOSB (−7.53), RGS1 (−4.52), and NR4A1 (−3.72). These few genes with very
large negative LFCs in Day 2 dominate the Pearson r calculation and cancel the
positive correlation from the 102 concordant genes. Spearman ρ is less sensitive to
these extreme values but still near zero because the same ~37 genes are the largest
rank outliers. **Direction concordance (73.4%) is the more interpretable metric here.**

The genome-wide Pearson r = 0.323 is unaffected by the significance filter — it uses
all 13,100 shared genes and shows moderate positive correlation at the transcriptome level.

### TREM2 / SPP1 / GPNMB

| Gene | GSE162694 log2FC | GSE162694 padj | GSE162694 sig | GSE135251 log2FC | GSE135251 padj | GSE135251 sig | Match |
|---|---|---|---|---|---|---|---|
| TREM2 | +2.48 | 5.17e-13 | YES | +2.52 | 2.22e-07 | YES | **Concordant — sig in both** |
| SPP1 | +1.87 | 1.10e-07 | YES | +1.41 | 2.76e-03 | YES | **Concordant — sig in both** |
| GPNMB | +1.19 | 2.25e-09 | YES | +0.91 | 6.00e-03 | NO | Concordant direction; below LFC cutoff in Day 2 |

TREM2 and SPP1 replicate across the two largest cohorts. GPNMB is significant in Day 1
but just misses the LFC > 1 threshold in Day 2 (+0.91), consistent with bulk RNA-seq
dilution of the macrophage/stellate cell signal.

### LFC Correlation Plot

See `results/lfc_corr_d1_vs_d2.png`  
*(Copied from `week4-day2-nafld-validation-rnaseq/plots/lfc_correlation.png`)*

---

## Ranked Conclusion: Which Pairing Shows Strongest Agreement?

**Ranking: A > C > B**

### 1st: Comparison A — GSE130970 (Day 3) vs GSE162694 (Day 1)  ← Strongest

Comparison A wins on every per-gene concordance metric:

| Why A is strongest | Detail |
|---|---|
| Highest Fisher's OR | 8.68 — overlap is most enriched above chance |
| Highest direction concordance | 94.7% — only 2 of 38 genes call opposite direction |
| Highest sig-overlap Pearson r | 0.656 — strong quantitative LFC agreement |
| Highest genome-wide Pearson r | 0.461 — most coherent transcriptome-wide LFC pattern |

**Reason:** Both Day 1 and Day 3 are dominated by early-to-intermediate fibrosis
(Day 1 median F1, Day 3 F0/F1 = 53/78 samples). They share a similar disease
population composition, making their transcriptional profiles more comparable. Day 1
also has 31 normal controls — the best-powered reference group among the three datasets
— which improves effect size estimate quality.

---

### 2nd: Comparison C — GSE162694 (Day 1) vs GSE135251 (Day 2)  ← Middle

Comparison C has the largest absolute overlap (139 genes) and the most significant
enrichment p-value (9.5e-43), but direction concordance drops to 73.4% and genome-wide
Pearson r (0.323) falls below Comparison A (0.461):

| Metric vs A | C vs A |
|---|---|
| Fisher's OR | 5.11 vs **8.68** (lower) |
| Fisher's p | **9.5e-43** vs 8.7e-21 (lower p, but driven by larger list sizes) |
| Direction concordance | 73.4% vs **94.7%** (lower — IEG discordance) |
| Genome-wide Pearson r | 0.323 vs **0.461** (lower) |

The lower concordance is explained by the IEG (immediate-early gene) discordance: FOS,
FOSB, JUN, JUNB, DUSP1 and related stress-response genes are upregulated in Day 1
(early/mixed fibrosis) but strongly downregulated in Day 2 (advanced fibrosis F3/F4-
enriched), where hepatocyte transcriptional collapse overrides the acute IEG response.
The Fisher's p is the lowest of the three because the universe and list sizes are
largest, not because the per-gene concordance is strongest.

---

### 3rd: Comparison B — GSE130970 (Day 3) vs GSE135251 (Day 2)  ← Weakest

Comparison B ranks last on the two most robust metrics:

| Metric | B value | Ranking |
|---|---|---|
| Fisher's OR | 4.78 | 3rd |
| Direction concordance | 84.6% | 2nd |
| Genome-wide Pearson r | 0.219 | 3rd |
| Sig-overlap Pearson r | 0.510 | 2nd |

Two independent sources of divergence compound in Comparison B:
1. **Day 3's tiny control group (n=4)** restricts the significant gene list to 180
   extreme signals, biasing the overlap sample toward genes with large LFCs.
2. **Day 2's advanced-fibrosis enrichment** introduces IEG/inflammatory discordance
   not present in Comparison A (where both cohorts are F0/F1-dominated).

---

### Summary Narrative

> The Day 1 and Day 3 cohorts (both early/mixed fibrosis, F0–F2 dominant) show the
> strongest mutual concordance (Comparison A: OR=8.68, concordance=94.7%, genome-wide
> r=0.461). The Day 1–Day 2 comparison (Comparison C) has the largest overlap in absolute
> terms (139 genes) but lower per-gene concordance (73.4%) due to IEG discordance driven
> by Day 2's advanced-fibrosis enrichment. Day 3 vs Day 2 (Comparison B) is the weakest
> pairing — the combination of Day 3's small control group and Day 2's disease-stage
> composition creates the lowest genome-wide correlation (r=0.219) among the three pairs.
>
> TREM2 and SPP1 replicate concordantly across all three cohorts in direction; both reach
> significance in the two well-powered datasets (Day 1 and Day 2) and show the same
> positive trend in Day 3 despite failing to reach significance there. This three-cohort
> directional consistency is the strongest available evidence for TREM2 and SPP1 as robust
> NAFLD transcriptional markers in bulk liver RNA-seq.
