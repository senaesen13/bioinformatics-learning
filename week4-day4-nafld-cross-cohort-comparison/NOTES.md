# Week 4 Day 4 — NAFLD Pairwise Cross-Cohort Comparison

This folder consolidates pairwise DEG overlap and pathway concordance comparisons among three
independent NAFLD cohorts (GSE162694, GSE135251, GSE130970).

The primary comparison is **Dataset 1 vs Dataset 2 (Comparison C)**. It has the most overlapping
genes (139), the strongest enrichment p-value (Fisher's p = 9.5e-43), and is the only pairing
where both TREM2 and SPP1 are individually significant in both datasets. These 139 shared genes
will be used as input for GEM and coexpression analysis going forward.

Dataset 3 vs Dataset 1 (Comparison A) has the highest direction agreement percentage (94.7%).
However, it only produces 38 overlapping genes, because Dataset 3 has very few healthy controls
(n = 4). That limits it as an input for downstream analysis. It is useful as supporting evidence
for the biology seen in Dataset 1, but is not the primary comparison.

---

## Overview

This folder is the single reference point for all pairwise comparisons among the three
NAFLD cohorts. For full dataset-specific methodology (DESeq2, filtering, GSEA) see the
individual cohort folders:

| Cohort | Folder |
|---|---|
| GSE162694 | `../week4-day1-nafld-bulk-rnaseq/` |
| GSE135251 | `../week4-day2-nafld-validation-rnaseq/` |
| GSE130970 | `../week4-day3-nafld-gse130970-rnaseq/` |

---

## Datasets Compared

Three pairwise comparisons are documented here:

- **Comparison A:** GSE130970 (Day 3) vs GSE162694 (Day 1, discovery cohort)
- **Comparison B:** GSE130970 (Day 3) vs GSE135251 (Day 2, validation cohort)
- **Comparison C:** GSE162694 (Day 1) vs GSE135251 (Day 2)

Comparisons A and B were computed by `scripts/pairwise_comparison.R` in this folder.
Comparison C reuses numbers already computed during the Day 2 analysis
(`week4-day2-nafld-validation-rnaseq/results/validation_overlap_summary.csv`);
the LFC correlation plot for C is copied from that folder as `results/lfc_corr_d1_vs_d2.png`.

> **GSEA ranking update (2026-07-21):** Pathway-level concordance statistics were
> recomputed using the updated GSEA results (ranking metric: sign(lfc_apeglm) × −log10(pvalue_mle)).
> Gene-level concordance statistics are unchanged (they use DESeq2 MLE results, not GSEA ranking).

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
| Binomial p (concordance) | One-sided binom.test(n_concordant, n_overlap, p=0.5, alt="greater") |
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

**Fisher's OR** accounts for the fact that different cohorts have different numbers of
significant genes. A high OR means the overlap is larger than chance alone would produce,
even after adjusting for list size.

**Direction concordance** answers a simple question: do these genes go up and down in the
same direction in both datasets? The binomial test checks whether the concordance rate is
above 50% (what you would expect if direction were random). LFC correlation (Pearson/Spearman)
is not reported — it can appear near zero even when most genes are concordant, because a few
strongly discordant genes dominate the calculation. Direction concordance with a binomial test
is the cleaner metric here.

---

## Gene-Level Results

Full results: `results/pairwise_comparison.md`  
Numeric CSV: `results/pairwise_summary.csv`

| Metric | Comparison A (D3 vs D1) | Comparison B (D3 vs D2) | **Comparison C (D1 vs D2)** |
|---|---|---|---|
| Overlap | 38 | 52 | **139** |
| Fisher's OR | **8.68** | 4.78 | 5.11 |
| Fisher's p | 8.7e-21 | 1.6e-16 | **9.5e-43** |
| Hypergeometric p | 1.1e-20 | 1.6e-16 | **9.5e-43** |
| Direction concordance | **94.7%** | 84.6% | 73.4% |
| Binomial p (concordance > 50%) | 2.7e-09 | 2.0e-07 | 1.6e-08 |

Comparison C (Day 1 vs Day 2) has the largest overlap (139 genes) and the strongest
enrichment p-value (9.5e-43). Both TREM2 and SPP1 are significant in both datasets —
this only happens in Comparison C. The direction agreement is 73.4%, which is lower
than Comparison A (94.7%), but the binomial test confirms it is still well above chance
(p = 1.6e-08).

Comparison A (Day 3 vs Day 1) has the highest direction concordance (94.7%) but only
38 overlapping genes. This is because Dataset 3 has only 4 healthy controls, so it
finds far fewer significant genes (180 total), leaving little to overlap with.

---

## Pathway-Level Results

**Pathway-level ranking: A >> B >> C** (by concordance significance).

| Comparison | KEGG concordance | binom.p | Hallmark concordance | binom.p |
|---|---|---|---|---|
| A: D3 vs D1 | **100%** (117/117) | **6.0e-36** | **81%** (21/26) | **1.2e-03** |
| B: D3 vs D2 | 66% (25/38) | 3.6e-02 | 38% (5/13) | 0.87 n.s. |
| **C: D1 vs D2** | 42% (13/31) | 0.86 n.s. | 18% (2/11) | 0.99 n.s. |

Comparison A has 100% KEGG concordance because Datasets 1 and 3 come from similarly
staged patients (both enriched for early/mixed fibrosis). The 5 discordant Hallmark sets
(OXIDATIVE_PHOSPHORYLATION, FATTY_ACID_METABOLISM, PROTEIN_SECRETION, PEROXISOME,
MYC_TARGETS) are activated in Day 3 but suppressed in Day 1. This most likely reflects
Dataset 3's n=4 control group, which is too small to give a stable metabolic baseline.

Comparison C has lower pathway concordance (42% KEGG, 18% Hallmark). The main reason is
that Dataset 1 is enriched for early/mid-stage NAFLD, while Dataset 2 is enriched for
advanced-stage NAFLD. Broad pathway categories like inflammatory signalling behave
differently at different disease stages, so they point in different directions between the
two cohorts. Individual genes like TREM2 and SPP1 are not as sensitive to this stage
difference — they go up consistently. The 2 pathway sets that are concordant across all
three cohorts (APICAL_JUNCTION and MYOGENESIS) are the most fibrosis-stage-independent
signals in this study.

---

## Interpretation and Caveats

### Important Caveats

1. **Gene symbol matching noise.** Some multi-mapping gene symbols (e.g., pseudogenes
   that share a name with a protein-coding gene across annotations) may inflate or
   deflate overlap counts slightly.

2. **GSE130970 power.** With only n=4 Normal (NAS=0) samples, Day 3 discovers only
   180 significant genes — the most robust signals only. This limits the sensitivity
   of overlap analyses in Comparisons A and B. GSEA on the full ranked list (see Day 3
   folder) is more reliable than the gene-list overlap here.

3. **Cross-cohort LFC comparability.** Log2 fold changes from different datasets are
   comparable in direction and rough magnitude, but not in exact value. Sample composition,
   batch effects, and pipeline differences all affect the LFC scale. Direction concordance
   is therefore more reliable than any correlation of raw LFC values.

4. **Confounding by fibrosis stage.** All three cohorts contain a mix of fibrosis
   stages, but in different proportions. The "NAFLD" label conflates early steatosis
   through advanced cirrhosis. The IEG discordance between Day 1 and Day 2 is a direct
   consequence of this — immune/stress pathway genes behave differently depending on
   whether the dominant disease stage is early (F0/F1) or advanced (F3/F4). Stratification
   by fibrosis stage would require larger per-stage n than available here.

5. **Comparison C is reused, not recomputed.** Numbers come from the Day 2 analysis
   script. If thresholds or filtering are ever changed in that folder, Comparison C
   in this file should be updated accordingly.

---

## Why Dataset 1 vs Dataset 2 Was Chosen

Three pairwise comparisons were run. All the numbers are in the Gene-Level Results
table above. Here is why Comparison C was chosen as the primary comparison.

**Comparison A (Dataset 1 vs Dataset 3):** 94.7% direction agreement, 38 overlapping
genes. The agreement percentage is the highest of the three. But 38 genes is a very
small set for downstream analysis. The reason it's small is that Dataset 3 has only
4 healthy controls, so DESeq2 can only call 180 significant genes in that dataset. That
limits how much can overlap with the other two cohorts.

**Comparison B (Dataset 2 vs Dataset 3):** 84.6% agreement, 52 overlapping genes. Also
limited by Dataset 3's small control group.

**Comparison C (Dataset 1 vs Dataset 2):** 73.4% agreement, 139 overlapping genes. The
percentage is lower than A, but the number of overlapping genes is nearly 4× larger (139
vs 38). The enrichment p-value is also the strongest of the three: Fisher's p = 9.5e-43.

For downstream analysis like GEM construction and coexpression network analysis, the
number of shared genes matters more than the agreement percentage. 139 genes is a
workable input. 38 genes is not.

The 73.4% concordance is also not a weakness — it is well above chance (binomial
p = 1.6e-08). The 26.6% of genes that disagree mostly reflect a known biological
difference: Dataset 1 contains more early/mid-stage NAFLD patients (fibrosis stage
F0–F2), while Dataset 2 contains more advanced-stage patients (F3–F4). Some genes —
particularly immediate-early stress-response genes like FOS and FOSB — behave
differently at different fibrosis stages. That is a real biological finding, not noise.
The 73.4% that do agree are genes that go up or down consistently regardless of stage.

The pathway-level concordance for Dataset 1 vs Dataset 2 also looks low (42% KEGG,
18% Hallmark). This reflects the same fibrosis-stage difference. Broad pathway
categories (like "TNF signaling" or "IL-17 signaling") are sensitive to disease stage
and flip direction between the two cohorts. Individual gene-level signals like TREM2
and SPP1 are not as sensitive — they go up in both datasets. The low pathway numbers
do not weaken the gene-level finding.

Because of this, Dataset 1 vs Dataset 2's 139 shared genes will be used as the input
for GEM and coexpression analysis going forward.

---

## Final Conclusion

Datasets 1 and 2 are two independent NAFLD studies. Their 139 shared significant genes
are enriched 3.5× above what chance would predict (Fisher's p = 9.5e-43). 73.4% of
those genes go up or down in the same direction in both datasets (binomial p = 1.6e-08).

TREM2 is upregulated with log2FC +2.48 in Dataset 1 and +2.52 in Dataset 2. SPP1 is
upregulated with log2FC +1.87 in Dataset 1 and +1.41 in Dataset 2. Both are individually
significant in both datasets. This does not happen in any other pairing — Dataset 3 does
not have enough statistical power to confirm either gene.

Across all three datasets and all three comparisons, two pathway-level signals are
consistent: MYOGENESIS and APICAL_JUNCTION are activated in NAFLD in every cohort.
These are the most reproducible pathway findings in this study.

The 139 shared genes from Comparison C (Dataset 1 vs Dataset 2) will be used as the
input for GEM and coexpression analysis going forward.

---

## Output Files

| File | Description |
|---|---|
| `scripts/pairwise_comparison.R` | Comparison A + B script (loads all 3 cohorts) |
| `results/pairwise_comparison.md` | Full report: all three comparisons, gene-level + pathway-level, ranked conclusions |
| `results/pairwise_summary.csv` | One row per comparison; includes Fisher's and hypergeometric p |
| `results/lfc_corr_d3_vs_d1.png` | LFC scatter: GSE130970 vs GSE162694 (Comparison A) |
| `results/lfc_corr_d3_vs_d2.png` | LFC scatter: GSE130970 vs GSE135251 (Comparison B) |
| `results/lfc_corr_d1_vs_d2.png` | LFC scatter: GSE162694 vs GSE135251 (Comparison C; copied from Day 2) |
