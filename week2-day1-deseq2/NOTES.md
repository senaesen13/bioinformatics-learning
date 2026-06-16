# Week 2 — Day 1 — DESeq2 Differential Expression Analysis

> Standalone DESeq2 walkthrough using the built-in `airway` dataset.
> Script: `deseq2_analysis.R`

---

## The Experiment

**Dataset:** airway (Bioconductor — Himes et al. 2014, PLOS ONE)
**Organism:** Human airway smooth muscle cells
**Question:** Which genes change when cells are treated with dexamethasone (a steroid)?

8 samples across 4 cell lines, 2 conditions each:

| Sample | Cell line | Treatment |
|--------|-----------|-----------|
| SRR1039508 | N61311 | untreated |
| SRR1039509 | N61311 | dexamethasone |
| SRR1039512 | N052611 | untreated |
| SRR1039513 | N052611 | dexamethasone |
| SRR1039516 | N080611 | untreated |
| SRR1039517 | N080611 | dexamethasone |
| SRR1039520 | N061011 | untreated |
| SRR1039521 | N061011 | dexamethasone |

---

## The Full Pipeline

```
Raw counts (63,677 genes × 8 samples)
  └─ Filter low counts → ~33,000 genes kept
       └─ Size factor normalisation (median-of-ratios)
            └─ Dispersion estimation (empirical Bayes shrinkage)
                 └─ Wald test (GLM per gene)
                      └─ LFC shrinkage (apeglm)
                           └─ Results table + volcano plot
                                └─ VST transformation
                                     └─ PCA + heatmaps
```

---

## Step 1 — Load the Data

The `airway` package ships a **SummarizedExperiment** object — a container that holds three things together:

- **assay** → the count matrix (rows = genes, columns = samples)
- **colData** → sample metadata (which sample is treated vs untreated, which cell line)
- **rowData** → gene metadata (Ensembl gene IDs)

```r
data("airway")
se <- airway
dim(assay(se))   # 63,677 genes × 8 samples
colData(se)      # see the "dex" and "cell" columns
```

The two columns that matter:
- `dex` — the condition: `"trt"` (treated) or `"untrt"` (untreated)
- `cell` — which of the 4 cell lines the sample came from

---

## Step 2 — Create a DESeqDataSet

DESeq2 wraps the data in its own object called a **DESeqDataSet**, and you give it a **design formula** that describes what you want to test.

```r
dds <- DESeqDataSet(se, design = ~ cell + dex)
dds$dex <- relevel(dds$dex, ref = "untrt")
```

### What does `~ cell + dex` mean?

Think of it as telling DESeq2:
> "Account for the fact that different cell lines have different baseline expression levels, then find genes that change because of the dexamethasone treatment."

The variable you want to **test** goes **last** in the formula. `cell` goes first as a **blocking factor** — it soaks up the variation due to cell line differences so it doesn't confuse the treatment signal.

### Why relevel?

By default R picks levels alphabetically, so "trt" would come before "untrt" and fold changes would be calculated backwards (untreated / treated). Setting `ref = "untrt"` ensures fold changes are always **treated ÷ untreated**, which is what we want.

---

## Step 3 — Filter Low-Count Genes

```r
keep <- rowSums(counts(dds)) >= 10
dds  <- dds[keep, ]
```

**Why?** Genes with near-zero counts across all samples are noise — they can't possibly be called significant, but they do count toward the multiple-testing correction, making it harder to find the real hits.

Rule of thumb: keep genes with at least 10 total reads across all samples.

> DESeq2 also does its own automatic filtering later (called "independent filtering") — this manual step just removes the absolute zeros early to save time.

---

## Step 4 — Normalisation (Size Factors)

Raw counts can't be directly compared between samples because **sequencing depth varies** — a sample sequenced twice as deeply will have roughly twice the counts for every gene. That's not biology, it's just a technical difference.

DESeq2 uses the **median-of-ratios** method to correct for this:

```
1. For each gene, calculate its geometric mean across all samples.

2. For each sample, divide each gene's count by that gene's geometric mean.
   → This gives a ratio per gene per sample.

3. The size factor for a sample = the MEDIAN of all those ratios.
   (The median is used because most genes are NOT changing — so the median
   reflects sequencing depth, not biology.)

4. Divide every count in a sample by its size factor → normalised counts.
```

```r
dds <- estimateSizeFactors(dds)
sizeFactors(dds)   # should all be close to 1.0
```

### Why not TPM or RPKM?

If a treatment massively upregulates one giant gene, it inflates the total read count. Dividing by the total would make every other gene look artificially smaller (downregulated). The median-of-ratios method is **robust** to this problem because a single gene can't shift the median.

---

## Step 5 — Dispersion Estimation

Before testing, DESeq2 needs to know how **noisy** each gene is — otherwise it can't tell if a difference between two groups is real or just random fluctuation.

This noise is measured by **dispersion** (α). High dispersion = noisier gene.

RNA-seq counts follow a **Negative Binomial** distribution (not a simple Poisson). The Negative Binomial has an extra parameter (dispersion) that captures variability beyond what Poisson predicts.

**The problem:** with only 8 samples, each gene has very little data to estimate its own dispersion reliably.

**DESeq2's solution — empirical Bayes shrinkage:**

```
1. Estimate a raw per-gene dispersion from the data alone (noisy).

2. Fit a smooth curve of "expected dispersion" vs mean count.
   (Low-count genes tend to be noisier — the curve captures this trend.)

3. Shrink each gene's raw estimate toward the curve.
   → Borrows strength from other genes. Well-measured genes barely move.
     Noisy low-count genes are pulled to a more sensible value.
```

```r
dds <- estimateDispersions(dds)
plotDispEsts(dds)
```

In the dispersion plot:
- **Black dots** = raw per-gene estimates
- **Red line** = fitted trend
- **Blue dots** = final shrunken estimates (used in the model)
- Circled outliers = genes with unusually high variance — DESeq2 keeps their raw estimate rather than shrinking

---

## Step 6 — Wald Test (Testing for Differential Expression)

```r
dds <- DESeq(dds)
```

This one function call runs steps 4, 5, and 6 together. It:
1. Fits a **Generalised Linear Model (GLM)** for each gene
2. Estimates a log2 fold change (β_dex) and its uncertainty (SE)
3. Performs the **Wald test** on each gene

### What is the Wald test?

DESeq2 models each gene as:

```
log(count) = intercept + β_cell + β_dex
```

`β_dex` is the log2 fold change for dexamethasone. We want to know: **is it significantly different from zero?**

The Wald test:
```
1. Estimate β_dex and its standard error (SE) from the data.
2. Compute the Wald statistic: W = β_dex / SE(β_dex)
3. Under H₀ (no effect), W follows a standard normal distribution.
4. Convert W to a two-tailed p-value.
```

A large |W| → small p-value → the gene is likely differentially expressed.

> **Alternative:** DESeq2 also has a **Likelihood Ratio Test (LRT)**, better for comparing more than two groups or complex multi-factor designs. For a simple two-group test, Wald is standard.

---

## Step 7 — Extract Results

```r
res <- results(dds, name = "dex_trt_vs_untrt", alpha = 0.05)
summary(res)
```

### What do the columns mean?

| Column | Meaning |
|--------|---------|
| `baseMean` | Average normalised count across all 8 samples |
| `log2FoldChange` | log2(treated / untreated) — positive = up in treated |
| `lfcSE` | Standard error of the LFC — how uncertain the estimate is |
| `stat` | Wald statistic = LFC / lfcSE |
| `pvalue` | Raw p-value from the Wald test |
| `padj` | **Adjusted p-value** (Benjamini-Hochberg) — accounts for multiple testing |

### Why do we need adjusted p-values?

We test ~33,000 genes simultaneously. Even if NO genes were truly changing, we'd expect 5% × 33,000 = **1,650 false positives** at p < 0.05 by chance alone.

Benjamini-Hochberg (BH) correction controls the **False Discovery Rate (FDR)**:
> "Of all the genes I call significant, at most 5% are expected to be false positives."

This is less conservative than Bonferroni (which controls the probability of even one false positive) but much better suited to genomics where we expect many true positives.

---

## Step 8 — LFC Shrinkage

```r
res_shrunk <- lfcShrink(dds, coef = "dex_trt_vs_untrt", type = "apeglm")
```

**The problem:** genes with very low counts have wildly unstable fold change estimates. A gene with 1 read vs 3 reads gets LFC ≈ 1.58 — but that's based on almost no data.

**apeglm shrinkage** (Zhu et al. 2019):
- Genes with **high counts and consistent signal** → LFC barely changes
- Genes with **low counts or high noise** → LFC shrunk toward zero

The p-values and padj are **unchanged** — only the LFC estimates are corrected.

Use shrunken LFCs for:
- Volcano plots
- Ranked gene lists (e.g. for GSEA)
- Any visualisation involving fold change

Use raw LFCs only if you specifically need the MLE estimate.

---

## Step 9 — MA Plot (Sanity Check)

```r
plotMA(res)         # raw LFC
plotMA(res_shrunk)  # shrunken LFC
```

An MA plot shows:
- **X axis (A):** mean expression (average count)
- **Y axis (M):** log2 fold change

What to look for:
- Most genes should cluster around **LFC = 0** (grey)
- Significant genes are coloured
- After shrinkage, the cloud should be much tighter at the **left** (low-count) end — the inflated fold changes are pulled in

---

## Step 10 — Volcano Plot

```r
ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), colour = sig)) + ...
ggsave("volcano_plot.png")
```

A volcano plot is the standard way to display differential expression results:

```
High on Y axis     → very significant (tiny padj)
Far left or right  → large fold change
Top-right corner   → strongly upregulated AND significant  ← most interesting
Top-left corner    → strongly downregulated AND significant ← most interesting
Bottom             → large fold change but not significant, or vice versa
```

Dashed lines mark the thresholds:
- **Vertical:** |log2FC| = 1 (i.e. 2-fold change)
- **Horizontal:** padj = 0.05

Top 20 significant genes are labelled automatically with `ggrepel` (labels push apart to avoid overlap).

---

## Step 11 — VST Transformation

```r
vst_data <- vst(dds, blind = TRUE)
```

Raw counts are **not suitable** for PCA or heatmaps because:
- Variance increases with mean (high-count genes dominate everything)
- Highly expressed genes visually drown out lowly expressed ones

**VST (Variance Stabilising Transformation)** applies a smooth function that makes variance roughly **constant** across the full range of mean counts.

`blind = TRUE` means the transformation ignores the experimental design — we want to see the raw data structure, not something that's already been shaped by the model.

> **rlog** is an alternative that applies more aggressive shrinkage and is more accurate for very small datasets (<10 samples), but VST is faster and works well here.

---

## Step 12 — PCA

```r
pca_plot <- plotPCA(vst_data, intgroup = c("dex", "cell"), returnData = TRUE)
```

PCA (Principal Component Analysis) takes thousands of gene measurements per sample and reduces them to just 2-3 numbers that capture the biggest sources of variation.

- **PC1** = the direction of *most* variance between samples
- **PC2** = the direction of *second most* variance, uncorrelated with PC1

Each dot is one sample. Samples that behave similarly cluster together.

### How to interpret the PCA

| What you see | What it means |
|---|---|
| Treated and untreated well-separated on PC1 | Treatment is the dominant signal — great |
| Cell lines separated on PC1, treatment on PC2 | Cell line variation is larger than treatment — OK, our design accounts for it |
| Samples from same condition cluster tightly | Low technical noise |
| One sample far from its group | Potential outlier — investigate |

---

## Step 13 — Sample Distance Heatmap

```r
sample_dists <- dist(t(assay(vst_data)))
pheatmap(as.matrix(sample_dists), ...)
```

Computes the **Euclidean distance** between every pair of samples (using all VST-transformed gene values). Visualised as a symmetric heatmap.

- **Dark blue** = samples that are very similar
- **White/light** = samples that are very different
- Treated samples should cluster together; untreated samples should cluster together

A useful sanity check before going further — an unexpected cluster may reveal a mislabelled or contaminated sample.

---

## Step 14 — Top DE Genes Heatmap

```r
mat_scaled <- t(scale(t(mat_top)))   # z-score each gene
pheatmap(mat_scaled, ...)
```

Takes the top 40 differentially expressed genes and shows their expression pattern across all 8 samples.

**Z-scoring** each gene (subtract its mean, divide by its standard deviation) puts all genes on the same scale — otherwise a gene expressed at 10,000 counts would visually dominate one expressed at 100 counts.

- **Red** = higher than average expression (for that gene)
- **Blue** = lower than average expression (for that gene)
- pheatmap clusters rows (genes) and columns (samples) automatically by hierarchical clustering — genes with similar patterns group together

---

## Step 15 — Export Results

```r
write.csv(res_export, "deseq2_results.csv")
write.csv(sig_genes,  "deseq2_significant_genes.csv")
```

Two CSV files are saved:
- `deseq2_results.csv` — full results for all tested genes, sorted by padj
- `deseq2_significant_genes.csv` — significant genes only (padj < 0.05 and |LFC| > 1)

---

## Key Concepts Reference

| Concept | Simple definition |
|---------|-------------------|
| **Count matrix** | Table of raw RNA-seq reads: rows = genes, columns = samples |
| **Size factor** | Per-sample scaling value that corrects for sequencing depth |
| **Median-of-ratios** | DESeq2's normalisation method — robust to highly expressed outlier genes |
| **Dispersion** | How much a gene's counts vary beyond what chance (Poisson) predicts |
| **Empirical Bayes** | Borrow information from all genes to stabilise estimates for individual genes |
| **Wald test** | Test whether LFC / SE is large enough to rule out zero effect |
| **log2FoldChange** | log2(treated / untreated) — positive means up in treated |
| **LFC shrinkage** | Pull noisy fold changes toward zero so plots and rankings are trustworthy |
| **padj** | Adjusted p-value — controls false discovery rate across all genes tested |
| **FDR** | False Discovery Rate: expected % of false positives among your significant calls |
| **VST** | Variance Stabilising Transformation — equalises variance before PCA/heatmaps |
| **PCA** | Reduces thousands of gene dimensions to 2–3 axes capturing maximum variance |
| **PC1** | The single axis that explains the most variation between samples |

---

## Output Files

| File | What it is |
|------|-----------|
| `volcano_plot.png` | Fold change vs significance for all genes |
| `pca_plot.png` | Sample clustering by treatment and cell line |
| `deseq2_results.csv` | Full results table (all tested genes) |
| `deseq2_significant_genes.csv` | Filtered to significant hits only |

---

## Connection to Week 1

In Week 1 Day 4 we ran DESeq2 on the **mouse heart attack** data (MI vs sham). The same core steps applied — the main differences here are:

| | Week 1 Day 4 | Week 2 Day 1 |
|-|---|---|
| Data input | tximport from Kallisto `.h5` files | Built-in `airway` SummarizedExperiment |
| Organism | Mouse (*Mus musculus*) | Human (*Homo sapiens*) |
| Design | `~ Treatment + Time` | `~ cell + dex` |
| Shrinkage | Not applied | apeglm |
| Visualisation | View() only | Volcano + PCA + heatmaps |

---

*Week 2, Day 1 — DESeq2 deep dive on human airway data*
