# Week 4 Day 5 – NAFLD NMR Metabolomics (MTBLS174)

## What is metabolomics?

Metabolomics measures the small molecules (metabolites) present in a biological sample — amino acids, sugars, fatty acids, organic acids, and more. These are the end-products of cellular processes, so the metabolome is like a snapshot of what the body is actually doing right now, not just what genes or proteins are available.

**NMR spectroscopy** (used here) identifies and quantifies metabolites by their unique magnetic resonance signals. It is quantitative, reproducible, and requires no labelling — making it ideal for clinical serum samples.

## Dataset: MTBLS174

- **Source:** MetaboLights study MTBLS174
- **Sample type:** Blood serum from NAFLD patients
- **Platform:** Bruker 600 MHz NMR spectrometer (1D CPMG pulse sequence)
- **Metabolites quantified:** 16 (amino acids, organic acids, choline metabolites, glucose)
- **Samples:** 18 patients with varying degrees of liver steatosis (fat accumulation)
- **Key variable:** Steatosis grade — percentage of hepatocytes containing fat

## Steatosis grades used

| Grade  | Steatosis % | Samples |
|--------|-------------|---------|
| Low    | < 20%       | 8       |
| Medium | 20–40%      | 7       |
| High   | > 40%       | 3       |

## Analysis steps

### 1. PCA (Principal Component Analysis)
- Plots each sample as a point in a 2D space capturing most metabolic variation
- **Result:** Samples do not cluster cleanly by steatosis grade on PC1/PC2, suggesting the metabolic differences between mild and severe steatosis are subtle at this level of profiling (16 metabolites, small cohort)

### 2. Differential metabolite analysis (limma, High vs Low)
- Tests each metabolite for abundance differences between High and Low steatosis groups
- **Result:** No metabolite reached FDR < 0.05 — expected with only 3 High-grade samples (underpowered)
- **Top trends (by nominal p-value):**
  - **Tyrosine ↑** in High steatosis (p = 0.09, logFC = +0.47) — aromatic amino acid linked to insulin resistance
  - **Formate ↑** in High steatosis (p = 0.20, logFC = +1.51) — one-carbon metabolism marker
  - **Acetate ↓** in High steatosis (p = 0.63, logFC = −0.62) — gut microbiome-derived SCFA
  - **Ethanol ↓** in High steatosis (p = 0.39, logFC = −3.96) — likely detection artifact in high-fat samples

### 3. Heatmap of top 15 metabolites
- Rows = metabolites, columns = samples, values = z-scores (normalised per metabolite)
- Shows the amino acid metabolites (Glutamine, Leucine, Valine) tend to be elevated in higher steatosis
- Acetate and Glucose are the most variable metabolites across individuals

### 4. Volcano plot
- x-axis = log2 fold change (High vs Low), y-axis = -log10(p-value)
- No points cross the significance thresholds — confirms the study is underpowered for this comparison
- Tyrosine and Formate are closest to significance

### 5. Pathway enrichment
- Mapped significant metabolites to curated KEGG pathway groups
- With 0 significant metabolites, enrichment is uninformative here
- The pathway framework is in place and would activate with a larger cohort

## Key biological takeaways

1. **Small NMR metabolomics panels have limited statistical power** — 16 metabolites and 18 samples is enough to see trends but not achieve significance after correction
2. **Tyrosine elevation** in high steatosis is biologically plausible: aromatic amino acid catabolism is impaired in insulin-resistant liver
3. **Acetate decrease** in severe steatosis may reflect altered gut-liver axis (reduced short-chain fatty acid production)
4. **NMR vs mass spectrometry:** NMR covers fewer metabolites (~50–100 vs 1000s for LC-MS) but is more reproducible and quantitative — better for clinical biomarker validation

## Files

```
scripts/nafld_nmr_metabolomics.R        Full analysis pipeline
data/m_MTBLS174_*_maf.tsv              Metabolite abundance matrix (16 metabolites × 18 samples)
data/s_MTBLS174.txt                     Sample metadata (steatosis %, gender, BMI, age)
data/a_MTBLS174_*.txt                   Assay metadata (NMR parameters)
plots/pca_steatosis.png                 PCA coloured by steatosis grade
plots/heatmap_top15.png                 Heatmap of top 15 metabolites
plots/volcano_plot.png                  Volcano plot High vs Low steatosis
plots/pathway_enrichment.png            Pathway enrichment bar chart
results/differential_metabolites_*.csv  Full limma results table
results/pathway_enrichment.csv          Pathway enrichment table
```
