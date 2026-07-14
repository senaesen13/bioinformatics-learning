# Week 4 Day 4 — NAFLD Pairwise Cross-Cohort Comparison

## Purpose

Compare GSE130970 (a third independent NAFLD cohort, analysed in Day 3) separately
against each of the two existing cohorts to quantify transcriptional concordance. Two
pairwise comparisons are run:

- **Comparison A:** GSE130970 vs GSE162694 (Day 1, discovery cohort)
- **Comparison B:** GSE130970 vs GSE135251 (Day 2, validation cohort)

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
| Direction concordance | % of overlap genes with sign(LFC_A) == sign(LFC_B) |
| Pearson r (sig overlap) | Correlation of MLE LFC for overlap genes only |
| Spearman ρ (sig overlap) | Rank correlation of MLE LFC for overlap genes only |
| Genome-wide Pearson r | Correlation of MLE LFC across all shared-tested genes |
| TREM2 / SPP1 / GPNMB | LFC and significance in each cohort |

The universe for Fisher's test is defined as the set of gene symbols present in
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

| Metric | Comparison A (D3 vs D1) | Comparison B (D3 vs D2) |
|---|---|---|
| Overlap | 38 | 52 |
| Fisher's OR | **8.68** | 4.78 |
| Fisher's p | 8.7e-21 | 1.6e-16 |
| Direction concordance | **94.7%** | 84.6% |
| Pearson r (sig overlap) | **0.656** | 0.510 |
| Genome-wide Pearson r | **0.461** | 0.219 |

**Conclusion:** GSE130970 (Day 3) is more concordant with GSE162694 (Day 1) than with
GSE135251 (Day 2) on every metric. The most likely explanation is that Day 1 and Day 3
both contain predominantly early-to-intermediate fibrosis (F0/F1/F2) while Day 2
is enriched for advanced fibrosis (F3/F4), producing a different transcriptional
programme dominated by ECM remodelling and hepatocyte loss.

---

## Output Files

| File | Description |
|---|---|
| `scripts/pairwise_comparison.R` | Full comparison script (loads all 3 cohorts, runs both comparisons) |
| `results/pairwise_comparison.md` | Full markdown report with tables and conclusion |
| `results/pairwise_summary.csv` | One-row-per-comparison numeric summary |
| `results/lfc_corr_d3_vs_d1.png` | LFC scatter: GSE130970 vs GSE162694 (Comparison A) |
| `results/lfc_corr_d3_vs_d2.png` | LFC scatter: GSE130970 vs GSE135251 (Comparison B) |

---

## Limitations

1. **Gene symbol matching noise.** Some multi-mapping gene symbols (e.g., pseudogenes
   that share a name with a protein-coding gene across annotations) may inflate or
   deflate overlap counts slightly.

2. **GSE130970 power.** With only n=4 Normal (NAS=0) samples, Day 3 discovers only
   180 significant genes — the most robust signals only. This limits the sensitivity
   of overlap analyses. GSEA on the full ranked list (see Day 3 folder) is more
   reliable than the gene-list overlap here.

3. **Cross-cohort LFC comparability.** MLE LFCs from different datasets are comparable
   in direction and rough magnitude but not in exact value (different sample compositions,
   batch effects, pipeline differences). Pearson r on LFC values should be treated as
   an index of agreement, not a direct exchange rate.

4. **Confounding by fibrosis stage.** All three cohorts contain a mix of fibrosis
   stages, but in different proportions. The "NAFLD" label conflates early steatosis
   through advanced cirrhosis. Stratification by fibrosis stage would require larger
   per-stage n than available here.
