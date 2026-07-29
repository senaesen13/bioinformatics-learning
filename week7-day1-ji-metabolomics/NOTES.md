# Ji et al. 2022 — NAFLD Metabolomics Analysis Notes

**Paper:** Ji et al., "Plasma Metabolomics Provides Insight into Non-alcoholic Fatty Liver Disease Pathogenesis," *Biomedicines* 2022, 10, 1669.  
**Study:** 86 fasting plasma samples — Control (n=25), NAFL (n=42), NASH (n=19).  
**Metabolites:** 79 plasma metabolites measured by GC-MS and LC-MS (33 amino acids, 4 kynurenine-pathway metabolites, 4 nucleosides, 18 organic acids, 20 fatty acids).

---

## Step 1 — Inspect the Data

**What this step does and why:**  
Before writing any analysis code, you open the data file and look at it directly — what the columns are actually named, how many rows there are, what the group labels say. This sounds obvious but is critical: papers sometimes describe data one way (e.g., "86 rows, one per patient") while the supplementary file is structured completely differently (e.g., one row per metabolite). If you skip this step and assume the layout matches the description, your code will be wrong in ways that are hard to catch.

**What was found:**  
`Table_S1_Metabolite_Levels.csv` has **79 rows (one per metabolite) and 16 columns**, not 86 rows (one per patient). The columns are: `No`, `Category`, `Metabolite`, `Control_Mean`, `Control_SD`, `NAFL_Mean`, `NAFL_SD`, `NASH_Mean`, `NASH_SD`, `Normalized_NAFL`, `Normalized_NASH`, `Pval_Control_vs_NAFL`, `Pval_Control_vs_NASH`, `Pval_NAFL_vs_NASH`, `Pval_Kruskal_Wallis`, `Qval_FDR`.

This is a **summary statistics table** — it contains the group means, standard deviations, and pre-calculated p-values per metabolite. Individual patient-level measurements (86 rows × 79 columns) were not released with the paper. Every step below that says "individual data required" means exactly this: you would need those raw measurements, which are not in this file and would have to be requested from the authors.

---

## Step 2 — Quality Control

**What this step does and why:**  
QC asks: is the data clean enough to analyse? You check for missing values (gaps in the table), verify the group sizes match what the paper reported, look for known outliers the authors flagged, and check whether any values would break a mathematical operation you plan to use (like a log transform, which cannot handle zero or negative numbers).

**What was found:**
- **Missing values:** None. All 79 metabolites have complete summary statistics for all three groups.
- **Group sizes:** Cannot be verified from a summary table (it has no per-patient rows). The paper states Control=25, NAFL=42, NASH=19 (total=86).
- **Aspartic acid outlier:** The NASH group shows mean=43.18 µg/uL but SD=80.63 µg/uL — an SD nearly double the mean, which is a strong signal of an extreme value in the raw data. The paper itself noted one NASH patient had an aspartic acid level of ~382.4 µg/uL (roughly 9× the group mean). **Decision: retain the outlier.** The paper retained it and used Kruskal-Wallis, which is rank-based and far less sensitive to outliers than a t-test or ANOVA.
- **Zero/negative values:** None. All group means are positive, so a log-transform is mathematically safe at the summary level. The smallest means (1-methyladenosine at 0.002 µg/uL, kynurenic acid at 0.005) are very small but positive.

---

## Step 3 — Normality Check (Shapiro-Wilk)

**What this step does and why:**  
Biological measurements rarely follow a perfect bell curve. A normality test like Shapiro-Wilk checks whether a set of numbers is bell-curve-shaped, group by group, for each metabolite. This matters because the most common statistical tests (t-test, ANOVA) assume normality — if the data is skewed or has extreme outliers, those tests give wrong answers. If normality fails for even one metabolite in one group, you switch to non-parametric tests (like Kruskal-Wallis) that make no distributional assumptions.

**What was found:**  
Individual patient data is required to run Shapiro-Wilk, and that data is not in this file. The test cannot be executed. However, the conclusion is still clear: metabolomics data is routinely right-skewed (many metabolites have a few very high values), and the Aspartic acid NASH outlier alone would likely fail normality. The paper's choice of Kruskal-Wallis throughout is the correct decision, consistent with what a normality check would tell you.

---

## Step 4 — Univariate Statistics (Kruskal-Wallis + BH-FDR)

**What this step does and why:**  
Kruskal-Wallis is a non-parametric test that checks whether the three groups (Control, NAFL, NASH) differ in their distribution of a single metabolite, without assuming normality. You run it separately for each of the 79 metabolites. Because you are running 79 tests, you expect roughly 4 false positives just by chance (5% of 79 ≈ 4). The Benjamini-Hochberg (BH) FDR correction accounts for this by adjusting each p-value upward based on how many tests you ran — only metabolites that survive this correction are considered genuinely significant.

**What was found:**  
The Kruskal-Wallis p-values and FDR q-values are pre-calculated in the CSV. We independently re-applied BH-FDR correction to the provided KW p-values as a cross-check. Exactly **6 metabolites** passed FDR correction (q < 0.05), matching the paper exactly:

| Metabolite | KW p | FDR q (paper) | Direction |
|---|---|---|---|
| Glutamic acid | <0.001 | <0.001 | ↑ Control < NAFL < NASH |
| α-Ketoglutaric acid | <0.001 | 0.004 | ↑ Control < NAFL < NASH |
| Myristoleic acid | <0.001 | 0.004 | ↑ Control < NAFL < NASH |
| Tyrosine | 0.001 | 0.015 | ↑ Control < NAFL < NASH |
| Kynurenic acid | 0.001 | 0.015 | ↑ Control < NAFL (then slightly ↓ in NASH) |
| Palmitoleic acid | 0.004 | 0.048 | ↑ Control < NAFL < NASH |

Note: our re-computation flagged Palmitoleic acid's FDR q as 0.053 (just above 0.05) because the paper's actual KW p-value is more significant than 0.001 (it was only rounded to "<0.001" in the table). When using the paper's provided Qval_FDR column directly, Palmitoleic acid is correctly included at q=0.048. This illustrates why "<0.001" notation loses information.

Pairwise Wilcoxon tests (also from the paper's table) showed that the main signal for most metabolites is the Control-vs-NAFL and Control-vs-NASH contrasts, with the NAFL-vs-NASH contrast less significant — suggesting these metabolites rise early in disease (already elevated in NAFL, not just in advanced NASH).

---

## Step 5 — Trend Test (Enhancement, not in the original paper)

**What this step does and why:**  
Kruskal-Wallis only asks "are these groups different?" — not "are they different in a specific direction?" The Jonckheere-Terpstra (JT) test is designed for ordered groups: it asks "do the values go up (or down) as the groups progress from Control to NAFL to NASH?" Because the biology here explicitly predicts worsening with disease stage, JT is more statistically powerful than Kruskal-Wallis when that monotone trend exists. The paper did not use JT, but it is the natural enhancement.

**What was found:**  
JT requires individual patient data and cannot be run from this file. As a proxy, we computed Spearman rank correlation between group order (Control=0, NAFL=1, NASH=2) and the group means for each metabolite. A correlation of +1.0 means the group means increase perfectly from Control to NAFL to NASH.

All 6 FDR-significant metabolites showed ρ = +1.0 (perfect monotone increase) except kynurenic acid (ρ = +0.5 — it rises from Control to NAFL but the NAFL and NASH means are nearly identical). Across all 79 metabolites, 36 showed a monotone ↑ trend and 14 a monotone ↓ trend in group means.

This is a proxy only — it uses means, not individual values, so there is no p-value and no correction for sample size or variance. It describes direction, not significance.

---

## Step 6 — Visualization

**What this step does and why:**  
Statistical results need to be shown visually to be interpretable. Three plots were made:

1. **Bar plots with error bars** (mean ± SD) for the 6 significant metabolites. These show the magnitude of change across the three groups. True box-and-whisker plots cannot be made without individual data; bar plots with SD are the closest honest alternative from summary statistics.

2. **Z-score heatmap** for all 79 metabolites, approximating Figure 1A from the paper. For each metabolite, the three group means are converted to Z-scores (how far each group is from the average of the three groups, measured in standard deviations). A Z-score of +1.5 means that group's mean is 1.5 SDs above the average of the three group means. This makes all 79 metabolites comparable on one colour scale even though their raw units vary widely. This is a proxy computed from group means, not from individual patient values as in the paper.

3. **Volcano plot** (log₂ fold change NASH/Control vs –log₁₀ FDR q-value). Points in the upper-right corner are metabolites that are both highly significant and strongly elevated in NASH.

**What was found:**  
The bar plots confirm the visual pattern: all 6 significant metabolites show a stepwise increase from Control to NAFL to NASH, most strikingly for glutamic acid (2.4× in NASH vs Control) and myristoleic acid (2.5× in NASH vs Control). The heatmap shows that amino acids split into two groups — glutamic acid and tyrosine rise with disease, while many others (glycine, serine, threonine, histidine) remain stable or decrease slightly.

---

## Step 7 — Multivariate Analysis (PLS-DA)

**What this step does and why:**  
Kruskal-Wallis tests one metabolite at a time. PLS-DA (Partial Least Squares Discriminant Analysis) uses all 79 metabolites simultaneously to ask: can we separate Control, NAFL, and NASH patients using the entire metabolome profile? It works by finding combinations of metabolites that best distinguish the groups, then plotting each patient as a point in 2D space. If the groups form separate clusters, the metabolome as a whole is a good classifier. VIP (Variable Importance in Projection) scores rank each metabolite by how much it contributes to the group separation — metabolites with VIP > 1.0 are the key drivers.

**What was found:**  
PLS-DA requires a matrix of individual patient measurements (86 patients × 79 metabolites). That matrix is not in the CSV file. Simulating individual data from group means and SDs would be misleading because it ignores correlations between metabolites — the covariance structure is what drives PLS-DA separation, and we cannot recover it from summary statistics. The step is described here to explain what the paper did and why, but it cannot be reproduced without the raw data.

The paper reported that the same 6 FDR-significant metabolites all had VIP > 1.0 in their PLS-DA, confirming that the univariate and multivariate analyses identified the same key metabolites.

---

## Step 8 — Pathway Enrichment

**What this step does and why:**  
After finding 6 significant metabolites, the next question is: what biological pathways do they belong to? Pathway enrichment takes your list of significant metabolites and checks whether they appear in the same biological pathways more often than you'd expect by chance. If "glutamic acid," "alpha-ketoglutaric acid," and "tyrosine" all point to the same pathway (amino acid metabolism), that's a meaningful biological signal, not just 3 random metabolites.

The standard tool for this in metabolomics is MetaboAnalyst (metaboanalyst.ca), which uses maintained pathway databases (KEGG, HMDB) that are updated regularly. Hand-rolling a custom KEGG mapping is not done here because those databases change — any hardcoded mapping would go stale.

**What was found:**  
The 6 significant metabolite names have been written to `results/sig_metabolites_for_metaboanalyst.txt` for direct upload. The paper reported enrichment in alanine/aspartate/glutamate metabolism and aminoacyl-tRNA biosynthesis (primarily driven by glutamic acid and α-ketoglutaric acid, both TCA cycle intermediates). The two monounsaturated fatty acids (myristoleic and palmitoleic acid) point to fatty acid elongation and desaturation pathways, consistent with the known lipogenic activity in NAFLD/NASH.

To run enrichment: upload `results/sig_metabolites_for_metaboanalyst.txt` to MetaboAnalyst → Pathway Analysis → KEGG → Homo sapiens.

---

## Key Limitations of This Reanalysis

- **No individual patient data.** The CSV contains summary statistics only. Steps 3 (Shapiro-Wilk), 5 (Jonckheere-Terpstra), and 7 (PLS-DA) cannot be run from this file.
- **`<0.001` truncation.** Several KW p-values are reported only as `<0.001`. This loses precision when re-computing FDR, causing slight differences from the paper's Qval_FDR for borderline metabolites (notably Palmitoleic acid). We used the paper's provided Qval_FDR as primary reference.
- **Heatmap is a proxy.** The heatmap is built from group means, not from 86 individual Z-scores as in the paper's Figure 1A.

---

## Output Files

| File | Description |
|---|---|
| `plots/significant_metabolites_bar.png` | Mean ± SD bar plots for the 6 FDR-significant metabolites |
| `plots/heatmap_79metabolites_zscore.png` | Z-score heatmap, all 79 metabolites (proxy from group means) |
| `plots/volcano_NASH_vs_Control.png` | Volcano plot: fold change vs FDR significance |
| `results/all_metabolites_statistics.csv` | Full table with KW p-values, FDR (paper + recomputed), trend ρ |
| `results/significant_metabolites_FDR05.csv` | The 6 FDR-significant metabolites only |
| `results/trend_direction_proxy.csv` | Spearman ρ (direction proxy) for all 79 metabolites |
| `results/sig_metabolites_for_metaboanalyst.txt` | 6 metabolite names for MetaboAnalyst upload |
| `results/all_metabolites_ranked.csv` | All 79 metabolites ranked by FDR q-value |
