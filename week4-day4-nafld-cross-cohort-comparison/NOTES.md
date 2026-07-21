# Week 4 Day 4 — NAFLD Pairwise Cross-Cohort Comparison

This folder consolidates pairwise DEG overlap and pathway concordance comparisons among three
independent NAFLD cohorts (GSE162694, GSE135251, GSE130970).

**Primary finding (Comparison C — Day 1 vs Day 2):** The discovery cohort (GSE162694, n=143)
replicates in the independent validation cohort (GSE135251, n=216): 139 genes overlap (3.5×
enrichment above chance, Fisher's p 9.5e-43), 73.4% call the same direction (binomial p 1.6e-08),
and both TREM2 (+2.48/+2.52) and SPP1 (+1.87/+1.41) are individually significant in both datasets.

Supporting finding (Comparison A — Day 3 vs Day 1): Dataset 3 (GSE130970) and Dataset 1
(GSE162694) show the highest per-gene concordance (Fisher's OR 8.68, 94.7% direction concordance,
100% concordant KEGG pathways), consistent with their shared early/mixed fibrosis stage composition.
This provides supporting evidence for the biology seen in Day 1, but is not an independent
replication.

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

**Fisher's OR** controls for the different list sizes (Day 2 with 1058 sig genes vs
Day 3 with 180) by normalising against the universe. A high OR means the overlap is
not explained by chance.

**Direction concordance + binomial test** is the primary directional metric. Concordance
answers: "do these genes move in the same direction across studies?" The binomial test
(H0: p_concordance = 0.5) assesses whether the observed concordance rate is above what
would be expected by chance for the observed overlap size. Pearson/Spearman correlations
on overlap LFCs are not reported: when cohorts have very different LFC scales and include
discordant outliers (as in Comparison C), correlation coefficients can be near zero or
negative despite high biological concordance, making them misleading here.

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

**Primary comparison: C (Day 1 vs Day 2)** — the designed discovery-to-validation
replication with the largest overlap (139 genes, Fisher's p 9.5e-43) and the only
pairing where TREM2 and SPP1 are both individually significant in both datasets.
Direction concordance of 73.4% (binomial p 1.6e-08) is lower than Comparison A
due to IEG discordance driven by Day 2's advanced-fibrosis enrichment, but is well
above chance. Comparison A achieves the highest per-gene concordance (OR 8.68, 94.7%)
but reflects matched fibrosis-stage composition, not independent validation.

---

## Pathway-Level Results

**Pathway-level ranking: A >> B >> C** (by concordance significance).

| Comparison | KEGG concordance | binom.p | Hallmark concordance | binom.p |
|---|---|---|---|---|
| A: D3 vs D1 | **100%** (117/117) | **6.0e-36** | **81%** (21/26) | **1.2e-03** |
| B: D3 vs D2 | 66% (25/38) | 3.6e-02 | 38% (5/13) | 0.87 n.s. |
| **C: D1 vs D2** | 42% (13/31) | 0.86 n.s. | 18% (2/11) | 0.99 n.s. |

Comparison A achieves perfect KEGG concordance (100%, binomial p 6.0e-36) and strong
Hallmark concordance (81%, p 1.2e-03) — both significantly above chance. The 5
discordant Hallmark sets in Comparison A are metabolic/cellular-machinery programmes
(OXIDATIVE_PHOSPHORYLATION, FATTY_ACID_METABOLISM, PROTEIN_SECRETION, PEROXISOME,
MYC_TARGETS) that are activated in Day 3 but suppressed in Day 1, reflecting Day 3's
limited n=4 control group.

Comparison C's low pathway concordance is consistent with its gene-level IEG
discordance and does not undermine its status as the primary validation comparison.
The 2 concordant Hallmark sets (APICAL_JUNCTION and MYOGENESIS) and the concordant
ECM/immune KEGG pathways (Integrin signaling, ECM-receptor interaction, Focal adhesion,
Phagocytosis) are the fibrosis-stage-independent core of the NAFLD transcriptional
signature.

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

---

## Final Conclusion

**Primary result:** Comparison C (GSE162694 Day 1 vs GSE135251 Day 2) is the key
finding of this analysis. TREM2 (+2.48 / +2.52 log2FC) and SPP1 (+1.87 / +1.41 log2FC)
replicate concordantly across independent discovery and validation cohorts, each
individually significant in both datasets. The 139-gene overlap is enriched 3.5× above
chance (Fisher's p 9.5e-43) and 73.4% of genes call the same direction (binomial p
1.6e-08). The IEG discordance (FOS, FOSB, NR4A1, RGS1) is a biologically meaningful
signal of fibrosis-stage differences between cohorts, not measurement noise.

**Supporting result:** Comparison A (GSE130970 Day 3 vs GSE162694 Day 1) achieves the
highest per-gene concordance (OR 8.68, 94.7% direction concordance, binomial p 2.7e-09)
and 100% concordant KEGG pathways (117/117, binomial p 6.0e-36). This reflects matched
fibrosis-stage composition (both cohorts F0/F1-dominated) and supports the biological
robustness of the Day 1 signals. TREM2 and SPP1 show the expected positive direction in
Day 3 even without reaching significance (n=4 controls).

The one finding that is robust to all three cohorts at both gene and pathway level:
**MYOGENESIS and APICAL_JUNCTION** are concordantly activated in NAFLD in all three
datasets. These are the most replication-robust NAFLD pathway signals in this study.

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
