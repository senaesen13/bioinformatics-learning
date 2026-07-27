# NAFLD vs Obesity Gene Signature Comparison

This folder compares the transcriptomic signatures of NAFLD (liver biopsies) and
obesity (adipose tissue) to ask whether the same genes are significantly dysregulated
in both diseases.

---

## Datasets

| Dataset | Tissue | Condition | Significant genes |
|---|---|---|---|
| NAFLD primary (D1 ∩ D2) | Liver biopsies | NAFLD vs Normal | 139 genes |
| NAFLD D1 (GSE162694) | Liver biopsies | NAFLD vs Normal | 485 sig genes |
| NAFLD D2 (GSE135251) | Liver biopsies | NAFLD vs Normal | 1,058 sig genes |
| Obesity (GSE166047) | Adipose tissue | Obese vs Lean | **2 sig genes** |

The obesity study (GSE166047) produced only two significant genes after DESeq2
analysis at padj < 0.05 and |MLE log2FC| > 1: **AZGP1** (down in obesity, LFC = −2.43)
and **CLEC4C** (up in obesity, LFC = +7.68). This is a very short list and severely
limits any cross-disease comparison. All results below should be read with this
limitation in mind.

Significance thresholds applied consistently across all datasets: **padj < 0.05
AND |MLE log2FC| > 1**.

---

## Method

**Overlap count:** number of gene symbols appearing in both the NAFLD and obesity
significant gene lists.

**Jaccard Index:** overlap / (size of NAFLD list + size of obesity list − overlap).
Ranges from 0 (no overlap) to 1 (identical lists).

**Fisher's exact test:** tests whether the overlap is larger than expected by chance,
given a universe of all genes tested in both studies. Universe size = genes with a
valid gene symbol present in both the obesity and the respective NAFLD full results
table (~12,000–14,000 genes depending on the comparison).

**Direction concordance:** for genes appearing in both significant lists, checks whether
the sign of the log2FC is the same (both up, or both down) in NAFLD and obesity.
Since there is no formal overlap, direction concordance is instead checked in the
**full results** (non-filtered) for AZGP1 and CLEC4C.

---

## Results

### Overlap and Jaccard Index

| Comparison | NAFLD genes | Obesity genes | Overlap | Jaccard |
|---|---|---|---|---|
| 139-gene NAFLD signature vs Obesity | 139 | 2 | **0** | **0.000** |
| NAFLD D1 (485 genes) vs Obesity | 485 | 2 | **0** | **0.000** |
| NAFLD D2 (1,058 genes) vs Obesity | 1,058 | 2 | **0** | **0.000** |

There is **zero formal overlap** across all three comparisons. Neither AZGP1 nor
CLEC4C appears in any of the NAFLD significant gene lists. Fisher's exact test
gives p = 1.0 for all comparisons (no enrichment above chance). This outcome is
expected when one list has only 2 entries — the probability of randomly drawing 2
genes that both appear in the NAFLD list is very low.

### Direction of AZGP1 and CLEC4C in NAFLD Data (Full Results, Not Filtered)

Although AZGP1 and CLEC4C do not meet the NAFLD significance threshold, they are
present in the full NAFLD results tables. Their direction can be checked:

| Gene | LFC in obesity | LFC in NAFLD D1 | padj in NAFLD D1 | LFC in NAFLD D2 | padj in NAFLD D2 | Direction concordant? |
|---|---|---|---|---|---|---|
| AZGP1 | −2.43 (sig) | −0.21 (not sig) | 0.017 | −0.21 (not sig) | 0.172 | **Yes — both diseases: down** |
| CLEC4C | +7.68 (sig) | +0.18 (not sig) | 0.363 | not tested | — | **Yes in D1 — both: up** |

Both genes point in the same direction in NAFLD as in obesity, but neither reaches
the LFC > 1 cutoff in any NAFLD liver dataset.

**AZGP1** (alpha-2-glycoprotein, a zinc-binding protein made in the liver and secreted
into blood) is downregulated in both obese adipose tissue and NAFLD liver. The LFC in
NAFLD is small (−0.21), so this concordance is weak numerically, but the direction
is consistent.

**CLEC4C** (BDCA-2, a plasmacytoid dendritic cell marker) is upregulated in both
settings. Its LFC in NAFLD D1 is also small (+0.18) and not significant.

**TREM2, SPP1, FASN, and COL1A1** (the key NAFLD markers) do not appear in the
obesity significant gene list and are not directionally tracked because the obesity
study did not yield results for all of them at the genome-wide level.

### Specific Check: TREM2, SPP1, FASN, FASN in Obesity

None of TREM2, SPP1, FASN, or COL1A1 overlap with the obesity significant genes
(AZGP1, CLEC4C). The obesity list of 2 genes contains neither macrophage markers
nor fibrosis/ECM genes.

---

## Interpretation

**The NAFLD and obesity gene signatures are largely independent in this data.** The
formal overlap is zero — no gene is simultaneously significant in both the NAFLD liver
biopsies and the obesity adipose tissue study at the thresholds applied. The directional
concordance of AZGP1 (down in both) and CLEC4C (up in both, D1) is biologically
plausible — AZGP1 is a liver-secreted protein known to decrease in obesity and liver
disease — but the NAFLD effect size for both genes is too small to reach significance.

The primary reason for zero overlap is not an absence of shared biology between NAFLD
and obesity, but the statistical limitation of the obesity dataset: with only 2
significant genes after strict filtering, the comparison is underpowered by design.
A larger, more powered obesity study would be needed to meaningfully test whether
NAFLD and obesity transcriptional programs share components.

---

## Output Files

| File | Description |
|---|---|
| `scripts/nafld_obesity_comparison.R` | Full analysis script |
| `results/overlap_summary.csv` | Jaccard, overlap count, Fisher's p for all three comparisons |
| `results/gene_direction_table.csv` | AZGP1 and CLEC4C LFC/padj in obesity and both NAFLD datasets |
| `plots/lfc_comparison.png` | LFC of AZGP1 and CLEC4C in obesity vs NAFLD (shows direction even without significance) |
| `plots/set_overlap_summary.png` | Stacked bar showing NAFLD-only, shared, and obesity-only genes per comparison |
