# Week 4 Day 2 — NAFLD RNA-seq Validation Cohort

## Dataset: GSE135251

| Field | Value |
|---|---|
| GEO Accession | GSE135251 |
| Publication | Govaere et al. (2020), *Science Translational Medicine* |
| Platform | Illumina HiSeq (GPL18573) |
| Tissue | Liver biopsy (human) |
| Total samples | 216 |
| Normal (healthy controls) | 10 |
| NAFLD (all stages NAFL + NASH F0–F4) | 206 |
| Alignment | GRCh38; raw HTSeq counts, one file per sample |
| Gene identifiers | Ensembl IDs (no version suffixes) |

**Sample breakdown (as in paper):**
control=10 · NAFL=51 · NASH_F0-F1=34 · NASH_F2=53 · NASH_F3=54 · NASH_F4=14

This dataset is the independent validation cohort for the discovery analysis in
`week4-day1-nafld-bulk-rnaseq/` (GSE162694, Suppli et al. 2021).

---

## Grouping Strategy

All NAFL and NASH samples (regardless of fibrosis stage F0–F4) are collapsed into
a single **"NAFLD"** group and compared against healthy controls (**"Normal"**).
This matches the binary Normal vs NAFLD design used in the discovery cohort.

The disease status column in GEO metadata is `disease:ch1`:
- `"disease: Control"` → Normal (n = 10)
- `"disease: NAFLD"` → NAFLD (n = 206)

---

## Methodology

Identical pipeline to the discovery cohort:

1. **Protein-coding filter** — applied before count filtering. Same gene list derived
   from biomaRt in the discovery cohort (`protein_coding_gene_list.rds`).
   Keeps only `gene_biotype == protein_coding` Ensembl IDs.

2. **Mean-count filter** — `rowMeans(counts) >= 10`.

3. **DESeq2** — design `~ condition` (reference = Normal). Both MLE (un-shrunk) and
   apeglm-shrunk LFC are computed and stored in the results table.

4. **Significance threshold** — `padj < 0.01` AND `|log2FC_MLE| > 2`.
   MLE LFC is used for threshold filtering; apeglm LFC is reported in the output table
   (consistent with the discovery cohort approach).

### Filtering Counts

| Stage | Genes |
|---|---|
| Starting count matrix (HTSeq GRCh38) | 64,253 |
| After protein-coding filter | 16,390 |
| After mean-count filter (rowMeans ≥ 10) | **13,166** |

---

## Significant Gene Results

**Threshold:** padj < 0.01 AND |MLE log2FC| > 2

| Direction | Count |
|---|---|
| Up in NAFLD | 63 |
| Down in NAFLD | 40 |
| **Total** | **103** |

The validation cohort yields more significant genes (103 vs 24 in discovery) despite
having only 10 controls, because the NAFLD group here is enriched for advanced fibrosis
(F3+F4 = 68/206 samples), which increases average effect size relative to controls.

---

## TREM2 / SPP1 / GPNMB

| Gene | MLE log2FC | padj | In sig list? |
|---|---|---|---|
| **TREM2** | +2.52 | 2.2e-07 | **YES** |
| SPP1 | +1.41 | 2.8e-03 | NO — LFC below 2× |
| GPNMB | +0.91 | 6.0e-03 | NO — LFC below 2× |

**TREM2** replicates as significant with concordant direction (upregulated in NAFLD
in both cohorts, similar magnitude: +2.48 discovery, +2.52 validation).

**SPP1 and GPNMB** are statistically significant (padj < 0.01 in discovery; padj < 0.05
in validation) and consistently upregulated, but their fold changes stay below 2× in
both cohorts. This is expected for bulk RNA-seq across mixed fibrosis stages — both
genes are expressed in macrophage/stellate cell subpopulations that are diluted within
whole-liver biopsies. They are genuine NAFLD markers but the bulk RNA-seq fold change
does not reach the stringent |LFC| > 2 threshold.

---

## Cross-Cohort Validation Summary

**7 genes significant in both cohorts** (Ensembl ID match):

| Gene | disc log2FC | disc padj | val log2FC | val padj | Concordant |
|---|---|---|---|---|---|
| TREM2 | +2.48 | 5.2e-13 | +2.52 | 2.2e-07 | YES |
| CHI3L1 | +2.23 | 1.3e-09 | +2.19 | 1.1e-04 | YES |
| STMN2 | +2.68 | 4.0e-08 | +2.57 | 3.6e-03 | YES |
| PDZK1IP1 | +2.20 | 1.4e-07 | +2.23 | 2.2e-04 | YES |
| PRAMEF10 | +3.52 | 2.5e-10 | +3.23 | 3.2e-05 | YES |
| NR4A1 | +2.04 | 8.5e-08 | −3.72 | 6.4e-35 | NO |
| RGS1 | +2.04 | 1.2e-06 | −4.52 | 1.6e-49 | NO |

NR4A1 and RGS1 show opposite directions between cohorts (up in discovery, strongly
down in validation). Both regulate immune cell function; the discordance likely reflects
cohort-level differences in inflammatory cell composition. The large and significant
downregulation in the validation cohort is noteworthy and may reflect hepatocyte
dedifferentiation effects in advanced fibrosis that dominate in GSE135251.

| Metric | Value |
|---|---|
| Discovery sig genes | 24 |
| Validation sig genes | 103 |
| Overlap | 7 genes |
| Overlap % of discovery | 29.2% |
| Fisher's exact OR | 55.52 |
| Fisher's exact p-value | 4.7e-10 |
| Direction concordance | 71.4% (5/7) |
| Spearman ρ (sig overlap, n=7) | 0.929 (p = 0.007) |
| Pearson r (sig overlap, n=7) | 0.624 (p = 0.134 — underpowered at n=7) |
| Pearson r (genome-wide, 13,100 genes) | 0.323 |

The overlap is highly enriched above chance (OR = 55.5, p = 4.7e-10). The low raw
overlap count (7 genes) is a consequence of the stringent |LFC| > 2 threshold
producing small gene lists. At the genome-wide level (13,100 shared tested genes),
the two cohorts show moderate positive correlation (Pearson r = 0.32), consistent
with independent patient populations studying the same disease but with different
fibrosis stage distributions.

---

## Output Files

| File | Description |
|---|---|
| `scripts/validation_gse135251.R` | Full pipeline: download → filter → DESeq2 → cross-cohort |
| `plots/validation_overlap.png` | Bar chart: discovery-only / both / validation-only gene counts |
| `plots/lfc_correlation.png` | LFC scatter (13,100 shared genes), purple = both sig, yellow = markers |
| `results/gse135251_results.csv` | Full DESeq2 results (MLE + apeglm LFC, gene symbol) |
| `results/validation_overlap_summary.csv` | Overlap stats, Fisher's test, concordance, correlation |
| `results/protein_coding_gene_list.rds` | Protein-coding Ensembl IDs (biomaRt, for pipeline reproducibility) |
| `data/GSE135251/` | Raw HTSeq count files — gitignored, re-downloadable via script |

---

## Tools and Citations

- **GEOquery** (Davis & Meltzer, 2007) — GEO data download
- **DESeq2** (Love, Huber & Anders, 2014) — differential expression
- **apeglm** (Zhu, Ibrahim & Love, 2018) — LFC shrinkage
- **org.Hs.eg.db** — Ensembl → gene symbol mapping
- **ggplot2** + **ggrepel** — visualisation
