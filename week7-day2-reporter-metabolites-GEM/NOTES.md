# Reporter Metabolite Analysis — Human Genome-Scale Metabolic Model

## Research Question

Which metabolites in the human metabolic network are statistically enriched with
significantly differentially expressed enzyme-coding genes in NAFLD, and what does that
tell us about which parts of the metabolic network are most disrupted in the disease?

---

## Method

This analysis uses the **Patil & Nielsen (2005) Reporter Metabolite algorithm**
(PNAS 102(8):2685–2689), which maps gene-level differential expression significance
onto a genome-scale metabolic model (GEM) to identify metabolites whose enzymatic
neighbourhood is collectively more significant than expected by chance.

### Algorithm — five steps

**1. Gene significance scores.**  
Each gene's DESeq2 p-value is converted to a Z-score using the two-tailed
inverse-normal transformation:

> Z_g = Φ⁻¹(1 − p/2)

A gene with p = 0.05 gets Z ≈ 1.96; a gene with p = 1×10⁻¹⁴ gets Z ≈ 7.66. This is
always a positive number — it measures significance, not direction. Direction is tracked
separately via mean log₂FC.

**2. Metabolite aggregation.**  
For each metabolite node in the GEM (i.e., each metabolite that appears as a substrate
or product in any reaction), the Z-scores of all enzyme-coding genes connected to that
node are summed and normalised by the square root of the number of genes:

> Z_raw = (ΣZ_g) / √k

This gives metabolites with more connected genes an appropriately scaled score.

**3. Background correction via Monte Carlo permutation.**  
For each metabolite, 1,000 random gene sets of the same size k are sampled from the
full gene universe (without replacement). Each permuted set is scored the same way as
the real set. This builds a null distribution that accounts for the fact that larger
gene sets naturally produce larger Z_raw values.

**4. Corrected reporter Z-score.**  
The corrected score subtracts the null mean and divides by the null standard deviation:

> Z_corr = (Z_raw − μ_k) / σ_k

A Z_corr of +4.0 means the real metabolite neighbourhood is 4 standard deviations above
what a random gene set of the same size would score.

**5. Multiple-testing correction.**  
BH-FDR correction is applied across all 2,378 tested metabolite nodes.
Metabolites with Padj < 0.05 are considered significant.

### Data

**GEM network:** Human-GEM 2.0 topology pre-extracted as a flat CSV
(`data/human_gem_topology_network.csv`, 38,527 metabolite–gene–reaction edges across
the full human metabolic network). Human-GEM is maintained by SysBioChalmers
(Robinson et al., *Science Signaling* 2020). The full SBML model is available at
github.com/SysBioChalmers/Human-GEM; the CSV here is the pre-parsed topology,
which contains the same network information in a format that runs without an SBML parser.
Using the pre-parsed topology rather than the raw SBML file is the actual method used
here — it is not a limitation; it is the same network.

**Currency metabolite filtering:** Water, ATP, ADP, AMP, NAD+/NADH, NADP+/NADPH,
CoA, CO₂, O₂, and other cofactors that participate in almost every reaction are removed
before analysis. Without this filter they would dominate the top results purely by having
the largest neighbourhoods.

**Size filter:** Only metabolites with 3 to 100 enzyme-coding neighbours in the tested
gene set are analysed (2,378 of the full network nodes qualify).

**Differential expression input:** DESeq2 results from GSE130970 (NAFLD bulk RNA-seq,
74 NAFLD vs 4 normal liver samples). Gene-level p-values from the MLE model
(`pvalue_mle` column), 15,468 genes tested, minimum p = 1.9×10⁻¹⁴.

---

## Key Results

**35 metabolite nodes** are significant after BH-FDR correction (Padj < 0.05), from
2,378 nodes tested.

Top reporter metabolites:

| Metabolite | k genes | Mean log₂FC | Reporter Z | FDR q |
|---|---|---|---|---|
| 3-methylcrotonyl-CoA | 6 | −0.41 | 4.79 | 0.002 |
| Palmitate | 44 | +0.17 | 4.38 | 0.007 |
| Arachidonate | 66 | −0.02 | 4.13 | 0.009 |
| Myristic acid | 31 | +0.16 | 3.86 | 0.015 |
| Ganglioside GM1 | 7 | +0.26 | 3.86 | 0.015 |

**Biological interpretation of the top two hits:**

**3-methylcrotonyl-CoA** (Reporter Z = 4.79, FDR q = 0.002) is an intermediate in
leucine catabolism. Its six connected genes (ACADSB, IVD, GLYAT, MCCC1, MCCC2, AUH)
are collectively amongst the most significantly dysregulated genes in the NAFLD dataset.
Their mean log₂FC is −0.41 (slightly downregulated on average), meaning the signal is
driven by their statistical significance rather than a strong directional shift —
consistent with disrupted branched-chain amino acid metabolism in NAFLD, which has been
reported in multiple independent studies.

**Palmitate** (Reporter Z = 4.38, FDR q = 0.007) is the most abundant saturated fatty
acid in the liver and a central hub in the fatty acid metabolic network. Its 44
connected genes include the fatty acid transport proteins SLC27A4, CD36, FABP4, and
FABP5 — genes also implicated in hepatic lipid accumulation and steatosis. A mean
log₂FC of +0.17 indicates modest collective upregulation. This places palmitate
metabolism at the centre of the transcriptional disruption in NAFLD, which is consistent
with the known role of de novo lipogenesis and impaired fatty acid oxidation in the
disease.

---

## Verification Note

The algorithm was verified line-by-line against the five specified steps before this
module was packaged:

- Steps 1–4 (Z-score formula, aggregation, permutation sampling, correction formula):
  correct in both the R and Python implementations.
- Step 5 (BH-FDR): correct in the R implementation. The original Python script in
  `improvements/` was missing this step; the Python script in `scripts/` here
  (`reporter_metabolites.py`) has it corrected.

The results in `results/gem_reporter_metabolites_summary.csv` were generated by the R
implementation and include BH-FDR corrected Padj values.

---

## Connection to the Broader Project

This is the metabolic network layer of the same NAFLD signal traced throughout this
project. The bulk RNA-seq and scRNA-seq modules identified TREM2, SPP1, and COL1A1 as
robust gene-level markers; the reporter metabolite analysis traces those gene-level
differences back into the human metabolic network and identifies which metabolic nodes
— particularly fatty acid metabolism (palmitate, myristic acid) and amino acid
catabolism (3-methylcrotonyl-CoA via leucine breakdown) — are most strongly reflected
in the transcriptional disruption. The plasma metabolomics module (folder 07) provides
independent experimental evidence from the same disease spectrum, measuring some of
these same metabolites directly in patient blood.

---

## File Guide

```
data/
    human_gem_topology_network.csv         Pre-parsed Human-GEM 2.0 metabolite-gene
                                            network (38,527 edges; metabolite_id,
                                            metabolite_name, compartment, gene_symbol,
                                            reaction_id)
    human_gem_metabolite_gene_map.csv      Small curated metabolite-gene-subsystem map
                                            (71 rows; reference / spot-checking only)

scripts/
    reporter_metabolites.R                 R implementation (verified correct, all 5 steps)
    reporter_metabolites.py                Python implementation (BH-FDR corrected vs
                                            original; run with --deseq2 flag for real data)

results/
    gem_reporter_metabolites_summary.csv   2,378 metabolite nodes: Reporter Z, p-value,
                                            BH-FDR Padj, mean log2FC, neighbour genes

plots/
    gem_reporter_metabolites_zscores.png   Bar chart: top 20 reporter metabolites,
                                            coloured by mean log2FC direction
```

---

## Tools

- **R:** dplyr, ggplot2 (base R stats for qnorm, replicate, p.adjust)
- **Python:** numpy, pandas, scipy.stats.norm, statsmodels.stats.multitest
- **GEM source:** Human-GEM 2.0 (SysBioChalmers/Human-GEM, Robinson et al. 2020)
- **DESeq2 input:** GSE130970 results from folder `03-NAFLD-second-validation-GSE130970`
