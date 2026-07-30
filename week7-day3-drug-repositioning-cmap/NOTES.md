# Drug Repositioning — LINCS CMap L1000 Signature Reversal

## Research Question

Given the transcriptomic signature of NAFLD from bulk RNA-seq (GSE130970), which
approved or investigational compounds have gene expression profiles that most strongly
reverse that signature — and are therefore candidates for drug repositioning?

---

## Method

### Connectivity Map (CMap) and LINCS L1000

The **Broad Institute Connectivity Map (CMap)** is a large compendium of gene expression
profiles measured after treating human cell lines with thousands of small molecules
(drugs and tool compounds). Each perturbation is profiled using the **LINCS L1000**
assay, which measures the expression of ~1,000 landmark genes.

The central idea of the approach is **signature reversal**: a drug is a repositioning
candidate if its transcriptional signature in cell lines is the *opposite* of the disease
signature. If NAFLD upregulates a set of genes and the drug downregulates those same
genes (and vice versa for downregulated genes), the drug could potentially correct the
disease-associated dysregulation.

The degree of reversal is quantified by the **NCS (Normalized Connectivity Score)**
from the weighted tau connectivity scoring (wtcs) metric. More negative NCS = stronger
reversal of the disease signature. The NCS ranges from approximately −2 to +2 in the
LINCS M2 database used here.

**Note on NCS vs Tau:** The Tau score (range −100 to +100) shown on the clue.io web
interface is derived by percentile-ranking an NCS against the full background distribution
of all perturbagen profiles. The raw downloaded GCT files contain NCS values, not Tau
scores. The Min NCS values reported here are the most negative NCS observed across all
cell line, dose, and time point combinations for each compound — they are directionally
equivalent to Tau (negative = reversal) but on a different scale.

### Query

- **Query name:** NAFLD_repositioning
- **Query date:** 30 July 2026
- **LINCS database:** level5_modz_n1165466x12328 (1,165,466 signatures × 12,328 samples)
- **Scoring metric:** wtcs (weighted tau connectivity score)
- **Reference:** Subramanian et al. (2017) *Cell* 171(6):1437–1452

### Query Gene Lists

The query signature was derived from the GSE130970 DESeq2 results (NAFLD vs. normal,
74 NAFLD vs. 4 normal liver samples; 15,482 genes tested):

- **UP gene set:** Top 150 upregulated genes by log2FC, filtered by padj < 0.05 and
  lfc > 0.5, sorted descending.
- **DOWN gene set:** Top 150 downregulated genes, filtered by padj < 0.05 and lfc < −0.5,
  sorted ascending.
- **Overlap between lists: 0 genes** (disjoint, as required by CMap).

Both query files are in `.grp` format (one gene symbol per line, no header, no quotes).

---

## Results

### Notes on Database Coverage

Two compounds from the original mechanistic selection are **not present** in the LINCS
L1000 database used here:

- **Obeticholic Acid** — not profiled in the M2 L1000 dataset. Replaced by **GW-4064**,
  the most widely used FXR agonist tool compound, which is present in the database.
- **Resmetirom** — FDA approved for NASH in 2024, too recent for this L1000 dataset.
  Replaced by **Liothyronine** (T3, the natural thyroid hormone), the only thyroid
  hormone receptor agonist present in the database.

Two compounds appear under different names:

- **Rapamycin** (original) = **Sirolimus** in the database (alternative INN name).
- **Pictilisib** (original) = **GDC-0941** in the database (development code used).

### Key Results

Results are ranked by Min NCS (most negative = strongest reversal in any cell line):

| Rank | Compound | Mechanism | Min NCS | Profiles (neg/total) |
|---|---|---|---|---|
| 1 | **Vorinostat** | HDAC Inhibitor | **−1.5791** | 500 / 3,286 |
| 2 | Sirolimus | mTOR Inhibitor | −1.5134 | 505 / 1,951 |
| 3 | GDC-0941 | PI3K Inhibitor | −1.4906 | 305 / 1,045 |
| 4 | Pioglitazone | PPAR-gamma Agonist | −1.3282 | 43 / 264 |
| 5 | Metformin | AMPK/Insulin Sensitizer | −1.3248 | 44 / 237 |
| 6 | GW-4064 | FXR Agonist | −1.1855 | 5 / 46 |
| 7 | Fenofibrate | PPAR-alpha Agonist | −1.1804 | 6 / 24 |
| 8 | Liothyronine | Thyroid Hormone Agonist | −0.8112 | 1 / 19 |

**"Profiles (neg/total)"** = number of cell line / dose / time point combinations
showing a negative NCS out of all profiles tested for that compound. A higher fraction
means the reversal signal is more consistent across biological contexts.

### Interpretation

**Vorinostat** (HDAC inhibitor) shows the most negative Min NCS (−1.5791), placing it
as the top computational repositioning candidate against the NAFLD transcriptomic
signature. This aligns with published evidence that HDAC inhibitors reduce hepatic
inflammation and fibrosis in NAFLD animal models, with the proposed mechanism involving
histone deacetylase-mediated epigenetic reprogramming of inflammatory gene expression.

**Sirolimus** (mTOR inhibitor, = rapamycin) and **GDC-0941** (PI3K inhibitor) both
show strong reversal signals (NCS ≈ −1.51 and −1.49 respectively) and are consistent
across a substantial fraction of profiles. mTOR/PI3K signalling is upregulated in NAFLD
and both compounds are known to reduce hepatic lipid accumulation in preclinical models.

**Fenofibrate** (PPAR-alpha agonist) is notable for having a **negative median NCS
(−1.01)** across its non-zero profiles, meaning it consistently reverses the NAFLD
signature rather than showing reversal in only a subset of contexts. Fenofibrate is
already in clinical use for hypertriglyceridaemia and has been investigated in NAFLD
trials.

**Liothyronine** (thyroid hormone) shows the weakest reversal (NCS = −0.81, only 1 of
19 profiles negative), consistent with the fact that non-selective thyroid hormone
agonists cause cardiac side effects at the doses needed for hepatic effect — thyroid
hormone receptor beta selectivity (the mechanism of resmetirom) was specifically
developed to avoid this.

**GW-4064** (FXR agonist) shows reversal in only 5 of 46 profiles (Min NCS = −1.19),
with the majority of profiles showing positive NCS. This likely reflects the strong
cell-line specificity of FXR expression, since FXR (NR1H4) is primarily a hepatocyte
receptor and most L1000 cell lines are non-hepatic.

### Comparison with Original Simulated Results

| Original (Simulated) | Simulated Tau | Real Compound | Real Min NCS |
|---|---|---|---|
| Obeticholic Acid (top) | −98.5 | **not in database** — | — |
| Resmetirom | −95.2 | **not in database** — | — |
| Pioglitazone | −91.4 | Pioglitazone | −1.3282 |
| Metformin | −86.8 | Metformin | −1.3248 |
| Fenofibrate | −82.1 | Fenofibrate | −1.1804 |
| Vorinostat | −77.5 | **Vorinostat (real top candidate)** | **−1.5791** |
| Pictilisib | −73.2 | GDC-0941 (same compound) | −1.4906 |
| Rapamycin | −68.9 | Sirolimus (same compound) | −1.5134 |

The simulated results ranked Obeticholic Acid first by design (it is a known NAFLD drug
target in clinical trials). The real L1000 query ranks **Vorinostat** first — an HDAC
inhibitor with substantially more cell-line coverage in the L1000 database (3,286
profiles) than any compound in this panel, giving it the most evidence for a reversal
signal. The two FXR/THR-beta agonists driving the simulated ranking are not in the
database, so the real result cannot confirm or deny their connectivity.

---

## Connection to the Broader Project

The query signature was built from the same GSE130970 DESeq2 results used in the
reporter metabolite analysis (folder 08), closing the loop between transcriptomics
(folders 01–05), metabolomics (folder 07), metabolic network analysis (folder 08),
and potential therapeutic translation here in folder 09.

---

## File Guide

```
scripts/
    drug_repositioning_cmap.R    Gene list extraction + .grp file generation

results/
    cmap_up_genes.grp            Top 150 NAFLD upregulated genes (CMap query input)
    cmap_down_genes.grp          Top 150 NAFLD downregulated genes (CMap query input)
    drug_repositioning_candidates.csv  8 compounds with real Min NCS scores
                                        (LINCS_L1000_REAL_QUERY_20260730)

plots/
    lincs_drug_repositioning_tau_scores.png  Bar chart: Min NCS per compound
```

---

## Tools

- **R:** dplyr, ggplot2
- **Python:** csv, matplotlib (for post-query CSV update and plot regeneration)
- **DESeq2 input:** GSE130970 results from `03-NAFLD-second-validation-GSE130970/`
- **LINCS L1000 query:** clue.io — query name NAFLD_repositioning, run 30 July 2026
- **LINCS database:** level5_modz, M2 (1,165,466 signatures)
- **Reference:** Subramanian et al. (2017) *Cell* 171(6):1437–1452
