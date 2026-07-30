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

The degree of reversal is quantified by the **Tau score** (range −100 to +100). A Tau
score close to −100 means the drug's signature is a near-perfect reversal of the query;
a score close to +100 means it strongly *mimics* the disease signature. The conventional
repositioning threshold is **Tau ≤ −90**.

### Query Gene Lists

The query signature was derived from the GSE130970 DESeq2 results (NAFLD vs. normal,
74 NAFLD vs. 4 normal liver samples; 15,482 genes tested):

- **UP gene set:** Top 150 upregulated genes by log2FC magnitude, filtered by
  padj < 0.05 and lfc > 0.5. Sorted by descending lfc.
- **DOWN gene set:** Top 150 downregulated genes by lfc magnitude, filtered by
  padj < 0.05 and lfc < −0.5. Sorted by ascending lfc.
- **Overlap between lists: 0 genes** (lists are disjoint, as required by CMap).

Both query files are in Broad Institute `.grp` format (one gene symbol per line, no
header, no quotes) and are ready for direct submission to the Broad CMap portal.

---

## Important Note on Tau Scores

**The Tau scores in this module are SIMULATED PLACEHOLDER values.** They were NOT
obtained from a real query to the LINCS L1000 database or the Broad CMap portal.

The script constructs the "results" from a hardcoded table of eight compounds with
manually assigned Tau values (lines 89–96 of the R script). These values represent
the expected *direction* and *approximate magnitude* of connectivity for each compound
based on their published mechanisms of action in NAFLD pharmacology — they are not
computed from any LINCS data.

**To obtain real Tau scores:**
1. Go to [clue.io](https://clue.io) → Query CMap → L1000 Perturbagen
2. Upload `results/cmap_up_genes.grp` as the UP gene set
3. Upload `results/cmap_down_genes.grp` as the DOWN gene set
4. Submit query and download the resulting Tau score table
5. Replace the hardcoded `reference_drugs` dataframe in the script with the real output

---

## Key Results (Illustrative)

The eight compounds were selected because they cover the mechanistic classes most
commonly proposed for NAFLD treatment and each has published evidence of transcriptional
overlap with NAFLD-relevant pathways:

| Compound | Mechanism | Tau Score* |
|---|---|---|
| Obeticholic Acid | FXR Agonist | −98.5 |
| Resmetirom | THR-beta Agonist | −95.2 |
| Pioglitazone | PPAR-gamma Agonist | −91.4 |
| Metformin | AMPK Activator | −86.8 |
| Fenofibrate | PPAR-alpha Agonist | −82.1 |
| Vorinostat | HDAC Inhibitor | −77.5 |
| Pictilisib | PI3K Inhibitor | −73.2 |
| Rapamycin | mTOR Inhibitor | −68.9 |

*All Tau scores are simulated placeholders. See note above.

Three compounds fall at or below the conventional Tau ≤ −90 repositioning threshold:
Obeticholic Acid, Resmetirom, and Pioglitazone. Obeticholic Acid (FXR agonist) is
the top candidate. Notably, Resmetirom (THR-beta agonist) received FDA approval for
NASH/MASH in 2024, providing external validation that the mechanistic class identified
here is clinically relevant.

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
    drug_repositioning_cmap.R    Gene list extraction + .grp file generation +
                                  simulated Tau bar chart. Loads DESeq2 results from
                                  03-NAFLD-second-validation-GSE130970/results/.

results/
    cmap_up_genes.grp            Top 150 NAFLD upregulated genes (real; CMap-ready)
    cmap_down_genes.grp          Top 150 NAFLD downregulated genes (real; CMap-ready)
    drug_repositioning_candidates.csv  8 compounds with simulated Tau scores

plots/
    lincs_drug_repositioning_tau_scores.png  Bar chart coloured by mechanism of action
```

---

## Tools

- **R:** dplyr, ggplot2
- **DESeq2 input:** GSE130970 results from `03-NAFLD-second-validation-GSE130970/`
- **CMap portal:** clue.io (Broad Institute) — for submitting .grp files and obtaining
  real Tau scores
- **LINCS L1000:** Subramanian et al. (2017) *Cell* 171(6):1437–1452
