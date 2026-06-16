# Week 2 — Day 3 — DESeq2: Mouse Myocardial Infarction (KCL Dataset)

> Script: `deseq2_mouse_mi.R`
> Data: `KCLModule2_2022/Transcriptomics/abundances/`
> Source: Yang Hong, KCL Module 2 Transcriptomics (2022)
> Repo: https://github.com/sysmedicine/KCLModule2_2022

This analysis follows the same DESeq2 pipeline as Week 2 Day 1 (airway).
This file only covers what is different: the dataset, design, and results.

---

## The Experiment

**Model organism:** Mouse (*Mus musculus*)
**Question:** Which genes change in the heart after a myocardial infarction (heart attack)?

Mice underwent surgical ligation of the left anterior descending (LAD) coronary artery
to induce MI, or sham surgery (same procedure without ligation) as a control.
Hearts were collected at two time points post-surgery.

**8 samples — 2 factors:**

| Sample name | SRR accession | Treatment | Time |
|-------------|---------------|-----------|------|
| sham_1dMI_1 | SRR6068402 | sham | Day 1 |
| sham_1dMI_2 | SRR6068403 | sham | Day 1 |
| MI_1dMI_1   | SRR6068404 | MI   | Day 1 |
| MI_1dMI_2   | SRR6068405 | MI   | Day 1 |
| sham_3dMI_1 | SRR6068406 | sham | Day 3 |
| sham_3dMI_2 | SRR6068407 | sham | Day 3 |
| MI_3dMI_1   | SRR6068408 | MI   | Day 3 |
| MI_3dMI_2   | SRR6068409 | MI   | Day 3 |

---

## Key Differences from the Airway Analysis

### 1. Input: Kallisto `.h5` files via tximport

The airway dataset came pre-loaded in R. Here the raw data is Kallisto output —
pseudoaligned transcript-level abundance estimates stored in binary `.h5` files.

`tximport` aggregates transcript-level estimates to gene level using a `tx2gene`
table (`proteinCodingGenes.Rda`) that maps Ensembl transcript IDs to Ensembl gene IDs.

```r
txi <- tximport(files, type = "kallisto", txOut = FALSE,
                tx2gene = proteinCodingGenes, ignoreTxVersion = TRUE)
```

`round()` is applied to the count matrix before passing it to DESeq2 because
tximport returns non-integer "estimated counts" (fractional, from EM). DESeq2
requires integer counts.

### 2. Two-factor design

```r
design = ~ Treatment + Time
```

The airway design had `~ cell + dex` (cell line as a batch effect, dex as the test).
Here both factors are biologically meaningful:
- **Treatment** — the main contrast of interest (MI vs sham)
- **Time** — a blocking factor that also has its own biological signal (day 1 vs day 3)

Including Time in the model controls for the fact that gene expression changes over
time regardless of MI, so the Treatment effect estimate is cleaner.

### 3. Testing a contrast, not a named coefficient

In the airway analysis, `dex` was the *last* term in the formula, so DESeq2
automatically named the coefficient `dex_trt_vs_untrt`.

Here `Treatment` is the *first* term, but DESeq2 still creates
`Treatment_MI_vs_sham` as a named coefficient (because sham is the reference level).
Both `results(contrast = ...)` and `lfcShrink(coef = "Treatment_MI_vs_sham")` work.

### 4. Stricter significance threshold

The original KCL workshop used `alpha = 0.01` (not 0.05).
With only 2 biological replicates per group, power is low and a tighter FDR threshold
reduces false positives. The volcano plot dashed line is also drawn at padj = 0.01.

### 5. Two PCA plots

Because there are two biological factors (Treatment and Time), we produce separate
PCA plots coloured by each. This lets you see:
- How much variance is explained by MI vs sham
- How much is explained by the day 1 vs day 3 time difference

### 6. Gene symbols on plots

Mouse Ensembl gene IDs (e.g. `ENSMUSG00000051951`) are replaced with readable
gene symbols (e.g. `Xkr4`) on the volcano plot and heatmap row labels, using the
`external_gene_name` column from `proteinCodingGenes.Rda`.

---

## Key Results (MI vs sham, FDR < 1%, |LFC| > 1)

| Gene | Direction | Known role in MI |
|------|-----------|-----------------|
| Spp1 (Osteopontin) | Up | Inflammation, cardiac remodelling |
| Cthrc1 | Up | Collagen secretion, heart tissue repair |
| Hspa1a | Down | Heat shock protein — stress response |
| Hspa1b | Down | Heat shock protein — stress response |
| Ncapg | Up | Cell cycle / DNA condensation |

The small number of significant genes (13 in the original workshop) reflects
the low replicate count (n=2 per group). Real cardiac studies use n≥3.

---

## Output Files

| File | Contents |
|------|---------|
| `plots/dispersion_estimates.png` | Dispersion plot |
| `plots/ma_plot.png` | Raw vs shrunken LFC side by side |
| `plots/volcano_plot.png` | MI vs sham, labelled with gene symbols |
| `plots/pca_treatment.png` | Samples coloured by MI / sham |
| `plots/pca_time.png` | Samples coloured by day 1 / day 3 |
| `plots/sample_distance_heatmap.png` | Euclidean distances between samples |
| `plots/top_genes_heatmap.png` | Top DE genes, z-scored VST counts |
| `deseq2_results_MI_vs_sham.csv` | Full results table |
| `deseq2_significant_genes.csv` | Significant hits only |

---

*Week 2, Day 3 — DESeq2 on mouse MI data using tximport + KCL dataset*
