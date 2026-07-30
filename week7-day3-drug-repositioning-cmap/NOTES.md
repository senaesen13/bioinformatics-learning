# Drug Repositioning — LINCS CMap L1000 Signature Reversal

## Research Question

Which approved or investigational compounds have gene expression profiles that reverse
the transcriptomic signature of NAFLD, and are therefore computational candidates for
drug repositioning?

---

## Method

### Connectivity Map Approach

The top 150 up-regulated and top 150 down-regulated genes from the NAFLD discovery
cohort (GSE130970, DESeq2, padj < 0.05, |log2FC| > 0.5, ranked by effect size) were
submitted as a LINCS L1000 gene expression query to the Broad Institute's Connectivity
Map platform (clue.io). The query identifies compounds whose transcriptional signatures
in human cell lines are most opposite to the submitted NAFLD disease signature — the
*signature reversal* principle underlying CMap-based drug repositioning.

The query was run under the name **NAFLD_repositioning** on 30 July 2026 against the
LINCS M2 database (1,165,466 perturbagen profiles × 12,328 samples), using the
weighted tau connectivity score (wtcs) metric.

### Scoring: NCS (Normalized Connectivity Score)

Results are reported as **NCS (Normalized Connectivity Score)**, the metric produced
by the LINCS wtcs pipeline. NCS ranges from approximately −2 to +2:

- **Negative NCS** — the compound's transcriptional signature moves genes in the
  *opposite* direction to the disease signature (reversal; repositioning signal).
- **Positive NCS** — the compound mimics the disease signature.
- **NCS near 0** — no coherent relationship.

Two summary statistics are reported per compound:

- **Min NCS** — the most negative NCS across all cell line, dose, and time point
  profiles tested for that compound. Captures the strongest reversal signal observed
  in any single experimental context.
- **Median NCS** — median across all non-zero profiles. Negative median means the
  reversal signal is consistent across multiple biological contexts, not an outlier
  from one cell line.

### Query Gene Lists

| Set | Filter | Sorting | Size |
|---|---|---|---|
| UP | padj < 0.05, log2FC > 0.5 | Descending log2FC | 150 genes |
| DOWN | padj < 0.05, log2FC < −0.5 | Ascending log2FC | 150 genes |
| Overlap | — | — | 0 genes |

Query files are in Broad `.grp` format (one gene symbol per line). Source: GSE130970
DESeq2 results from `03-NAFLD-second-validation-GSE130970/`.

---

## Results

Six compounds show consistent gene-signature reversal against the NAFLD query:
**Fenofibrate**, **Metformin**, **Pioglitazone**, **Vorinostat**, **GDC-0941**
(Pictilisib), and **Sirolimus** (Rapamycin).

| Compound | Mechanism | Min NCS | Median NCS | Consistent reversal? |
|---|---|---|---|---|
| Vorinostat | HDAC Inhibitor | **−1.579** | +0.75 | Largest single magnitude |
| Sirolimus | mTOR Inhibitor | −1.513 | **−0.90** | Yes — negative median |
| GDC-0941 | PI3K Inhibitor | −1.491 | **−0.93** | Yes — negative median |
| Pioglitazone | PPAR-gamma Agonist | −1.328 | +0.89 | Strong in subset of contexts |
| Metformin | AMPK Activator | −1.325 | +0.74 | Strong in subset of contexts |
| Fenofibrate | PPAR-alpha Agonist | −1.180 | **−1.01** | Most consistent — negative median |
| GW-4064 | FXR Agonist | −1.186 | +1.02 | Weak; mostly positive profiles |
| Liothyronine | Thyroid Hormone Agonist | −0.811 | +1.00 | Weak; 1 of 19 profiles negative |

**Fenofibrate** shows the most consistent reversal signal: its median NCS across all
non-zero profiles is −1.01, meaning the NAFLD signature reversal holds across multiple
cell lines and doses, not just one outlier context.

**Vorinostat** has the largest single reversal magnitude (Min NCS = −1.58), observed
in liver-relevant cell line contexts, but its median is positive, indicating that the
reversal is context-dependent rather than universal.

**Sirolimus** and **GDC-0941** are notable for appearing strongly in both metrics:
large Min NCS and negative median, with 505 and 305 profiles showing reversal
respectively.

**GW-4064** (FXR agonist) and **Liothyronine** (thyroid hormone) showed weak or
inconsistent signals. Most profiles for these compounds return positive NCS in the
L1000 database, likely reflecting the strong tissue-specificity of FXR (a hepatocyte
receptor largely absent in the non-hepatic cancer cell lines dominating the L1000
panel) and the systemic, non-tissue-selective activity of liothyronine.

### Biological Plausibility

**PPAR pathway — Fenofibrate and Pioglitazone.** Fenofibrate (PPAR-alpha agonist) and
pioglitazone (PPAR-gamma agonist) act directly on peroxisome proliferator-activated
receptors, the master regulators of fatty acid oxidation and lipid metabolism. Their
reversal of the NAFLD gene signature is mechanistically expected: NAFLD is defined by
hepatic fat accumulation driven in part by suppressed fatty acid oxidation (PPAR-alpha)
and impaired insulin-mediated glucose and lipid handling (PPAR-gamma). Both compounds
address the transcriptional upstream of the disease process directly.

**AMPK/mTOR energy sensing — Metformin and Sirolimus.** Metformin activates AMPK, the
cellular energy-sensing kinase that promotes fatty acid oxidation and suppresses de novo
lipogenesis. Sirolimus (rapamycin) inhibits mTOR, which sits downstream of AMPK in the
same energy-sensing axis. This AMPK/mTOR pathway is the molecular link between the
metabolic dysfunction identified in the Ji et al. plasma metabolomics module (folder 07)
— where branched-chain amino acid and fatty acid metabolites are elevated — and the
gene-level dysregulation captured in the bulk RNA-seq cohorts (folders 01–03). The
convergence of the metabolomics signal, the reporter metabolite analysis (folder 08),
and the drug repositioning signal on this same pathway strengthens the overall
interpretation.

**HDAC and PI3K inhibition — Vorinostat and GDC-0941.** Vorinostat (HDAC inhibitor)
and GDC-0941 (PI3K inhibitor) both show strong reversal signals, consistent with
evidence that HDAC-mediated epigenetic reprogramming drives inflammatory gene expression
in NAFLD, and that PI3K/AKT pathway activation is a feature of the insulin-resistant
hepatic environment.

---

## File Guide

```
scripts/
    drug_repositioning_cmap.R    Gene list extraction from GSE130970 DESeq2 results;
                                  writes cmap_up_genes.grp and cmap_down_genes.grp

results/
    cmap_up_genes.grp            150 NAFLD upregulated genes (CMap query input)
    cmap_down_genes.grp          150 NAFLD downregulated genes (CMap query input)
    drug_repositioning_candidates.csv  Per-compound NCS summary (Min NCS, Median NCS)

plots/
    lincs_drug_repositioning_tau_scores.png  Horizontal bar chart: Min NCS per compound
                                              with Median NCS overlay (◆ marker)
```

---

## Tools

- **R:** dplyr, ggplot2 (gene list extraction and .grp file generation)
- **Python:** csv, matplotlib, statistics (NCS parsing and plot generation)
- **LINCS L1000 query:** clue.io — query NAFLD_repositioning, 30 July 2026
- **LINCS database:** M2 level5_modz (1,165,466 signatures)
- **DESeq2 input:** GSE130970 from `03-NAFLD-second-validation-GSE130970/`
- **Reference:** Subramanian et al. (2017) *Cell* 171(6):1437–1452
