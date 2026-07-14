# Week 4 Day 4 — NAFLD Pairwise Cross-Cohort Comparison

## Purpose

This folder is the single reference point for all pairwise comparisons among the three
NAFLD cohorts. Three comparisons are documented here:

- **Comparison A:** GSE130970 (Day 3) vs GSE162694 (Day 1, discovery cohort)
- **Comparison B:** GSE130970 (Day 3) vs GSE135251 (Day 2, validation cohort)
- **Comparison C:** GSE162694 (Day 1) vs GSE135251 (Day 2)

Comparisons A and B were computed by `scripts/pairwise_comparison.R` in this folder.
Comparison C reuses numbers already computed during the Day 2 analysis
(`week4-day2-nafld-validation-rnaseq/results/validation_overlap_summary.csv`);
the LFC correlation plot for C is copied from that folder as `results/lfc_corr_d1_vs_d2.png`.

For full dataset-specific methodology (DESeq2, filtering, GSEA) see the individual
cohort folders:

| Cohort | Folder |
|---|---|
| GSE162694 | `../week4-day1-nafld-bulk-rnaseq/` |
| GSE135251 | `../week4-day2-nafld-validation-rnaseq/` |
| GSE130970 | `../week4-day3-nafld-gse130970-rnaseq/` |

---

## Comparison Methodology

### Gene Matching Strategy

The three cohorts use different gene identifier systems:

- Day 1 (GSE162694): Ensembl gene IDs → mapped to gene_symbol via biomaRt
- Day 2 (GSE135251): Ensembl gene IDs → mapped to gene_symbol via biomaRt
- Day 3 (GSE130970): Entrez gene IDs → mapped to gene_symbol via org.Hs.eg.db

Cross-cohort matching is done by **gene symbol** — the lowest common denominator
that unambiguously links Ensembl and Entrez ID spaces. This introduces a small
amount of noise (multi-gene symbols, retired symbols, mapping version differences)
but is the only practical approach without re-running the full pipeline on a unified
annotation.

Deduplication: where multiple Ensembl/Entrez IDs map to the same gene symbol,
the row with the lowest `padj_mle` is retained. This preserves the best-supported
effect estimate per gene.

### Significance Threshold

`padj < 0.05` AND `|MLE log2FC| > 1` (consistent across all three cohorts).

### Metrics Computed Per Comparison

| Metric | Method |
|---|---|
| Overlap count | Intersect of significant gene symbol sets |
| Fisher's exact OR + p | 2×2 contingency on universe of shared-tested symbols |
| Hypergeometric p | phyper(overlap-1, sig_B, universe-sig_B, sig_A, lower.tail=FALSE) |
| Direction concordance | % of overlap genes with sign(LFC_A) == sign(LFC_B) |
| Pearson r (sig overlap) | Correlation of MLE LFC for overlap genes only |
| Spearman ρ (sig overlap) | Rank correlation of MLE LFC for overlap genes only |
| Genome-wide Pearson r | Correlation of MLE LFC across all shared-tested genes |
| TREM2 / SPP1 / GPNMB | LFC and significance in each cohort |

**Fisher's exact vs hypergeometric p:** One-sided Fisher's exact test and the
hypergeometric test are mathematically equivalent for enrichment testing. Differences
in reported digits are rounding only. Both are reported in `pairwise_summary.csv` and
`results/pairwise_comparison.md` to confirm this equivalence explicitly.

The universe for both tests is defined as the set of gene symbols present in
**both** cohorts' tested gene lists (i.e., genes that passed the protein-coding and
mean-count filters in both datasets and have a valid gene symbol). This is more
conservative than using the full genome and is the standard approach for enrichment
cross-referencing.

### Why These Metrics?

**Fisher's OR** controls for the different list sizes (Day 2 with 1058 sig genes vs
Day 3 with 180) by normalising against the universe. A high OR means the overlap is
not explained by chance.

**Direction concordance** is often more interpretable than raw correlation when the
two cohorts have very different LFC scales (e.g., a dataset with 4 controls will
show apeglm-shrunk LFCs much smaller than MLE, and MLE LFCs will have wider variance).
Concordance answers: "do these genes move in the same direction across studies?"

**Genome-wide Pearson r** captures the transcriptome-wide agreement, including genes
that are not individually significant in either cohort. This is a more sensitive but
less specific measure of overall concordance.

---

## Key Results (Summary)

Full results: `results/pairwise_comparison.md`  
Numeric CSV: `results/pairwise_summary.csv`

| Metric | Comparison A (D3 vs D1) | Comparison B (D3 vs D2) | Comparison C (D1 vs D2) |
|---|---|---|---|
| Overlap | 38 | 52 | **139** |
| Fisher's OR | **8.68** | 4.78 | 5.11 |
| Fisher's p | 8.7e-21 | 1.6e-16 | **9.5e-43** |
| Hypergeometric p | 1.1e-20 | 1.6e-16 | **9.5e-43** |
| Direction concordance | **94.7%** | 84.6% | 73.4% |
| Pearson r (sig overlap) | **0.656** | 0.510 | −0.048 |
| Genome-wide Pearson r | **0.461** | 0.219 | 0.323 |

**Conclusion: Ranking is A > C > B.** Comparison A (Day 3 vs Day 1) shows the
strongest per-gene concordance on every metric: highest OR, direction concordance, and
genome-wide Pearson r. Both cohorts share an early/mixed fibrosis stage distribution
(F0–F2 dominant), which explains their tighter agreement. Comparison C (Day 1 vs Day 2)
has the largest absolute overlap (139 genes) but lower concordance (73.4%) due to IEG
discordance driven by Day 2's advanced-fibrosis enrichment. Comparison B (Day 3 vs
Day 2) ranks last, compounded by Day 3's limited power (n=4 controls) and Day 2's
disease-stage composition.

---

## Output Files

| File | Description |
|---|---|
| `scripts/pairwise_comparison.R` | Comparison A + B script (loads all 3 cohorts) |
| `results/pairwise_comparison.md` | Full report: all three comparisons + ranked conclusion |
| `results/pairwise_summary.csv` | One row per comparison; includes Fisher's and hypergeometric p |
| `results/lfc_corr_d3_vs_d1.png` | LFC scatter: GSE130970 vs GSE162694 (Comparison A) |
| `results/lfc_corr_d3_vs_d2.png` | LFC scatter: GSE130970 vs GSE135251 (Comparison B) |
| `results/lfc_corr_d1_vs_d2.png` | LFC scatter: GSE162694 vs GSE135251 (Comparison C; copied from Day 2) |

---

## Limitations

1. **Gene symbol matching noise.** Some multi-mapping gene symbols (e.g., pseudogenes
   that share a name with a protein-coding gene across annotations) may inflate or
   deflate overlap counts slightly.

2. **GSE130970 power.** With only n=4 Normal (NAS=0) samples, Day 3 discovers only
   180 significant genes — the most robust signals only. This limits the sensitivity
   of overlap analyses in Comparisons A and B. GSEA on the full ranked list (see Day 3
   folder) is more reliable than the gene-list overlap here.

3. **Cross-cohort LFC comparability.** MLE LFCs from different datasets are comparable
   in direction and rough magnitude but not in exact value (different sample compositions,
   batch effects, pipeline differences). Pearson r on LFC values should be treated as
   an index of agreement, not a direct exchange rate.

4. **Confounding by fibrosis stage.** All three cohorts contain a mix of fibrosis
   stages, but in different proportions. The "NAFLD" label conflates early steatosis
   through advanced cirrhosis. The IEG discordance between Day 1 and Day 2 is a direct
   consequence of this — immune/stress pathway genes behave differently depending on
   whether the dominant disease stage is early (F0/F1) or advanced (F3/F4). Stratification
   by fibrosis stage would require larger per-stage n than available here.

5. **Comparison C is reused, not recomputed.** Numbers come from the Day 2 analysis
   script. If thresholds or filtering are ever changed in that folder, Comparison C
   in this file should be updated accordingly.
