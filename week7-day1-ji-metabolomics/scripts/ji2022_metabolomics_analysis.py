#!/usr/bin/env python3
"""
Ji et al. 2022 NAFLD Metabolomics — Full Statistical Analysis
Paper : Biomedicines 2022, 10, 1669 — doi:10.3390/biomedicines10071669
Study : 86 plasma samples — Control (n=25), NAFL (n=42), NASH (n=19)
        79 plasma metabolites measured by GC-MS / LC-MS

DATA SOURCE NOTE
---------------
Table_S1_Metabolite_Levels.csv contains SUMMARY STATISTICS only:
group means, SDs, and pre-calculated p-values. Individual patient-level
measurements were not released with the paper. Steps that strictly
require individual data (Shapiro-Wilk, Jonckheere-Terpstra, PLS-DA)
are clearly labelled with what would be done and why they cannot be
executed here.
"""

import os
import sys
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from scipy import stats
from statsmodels.stats.multitest import multipletests
import warnings
warnings.filterwarnings("ignore")

# ── paths ────────────────────────────────────────────────────────────────────
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BASE_DIR   = os.path.dirname(SCRIPT_DIR)
CSV_PATH   = os.path.join(BASE_DIR, "Table_S1_Metabolite_Levels.csv")
PLOTS_DIR  = os.path.join(BASE_DIR, "plots")
RESULTS_DIR = os.path.join(BASE_DIR, "results")

os.makedirs(PLOTS_DIR, exist_ok=True)
os.makedirs(RESULTS_DIR, exist_ok=True)


# ═══════════════════════════════════════════════════════════════════════════
# STEP 1 — INSPECT THE DATA
# ═══════════════════════════════════════════════════════════════════════════
print("\n" + "="*70)
print("STEP 1 — DATA INSPECTION")
print("="*70)

raw = pd.read_csv(CSV_PATH)

print("\nColumn names:")
for c in raw.columns:
    print(f"  {c}")

print(f"\nDimensions: {raw.shape[0]} rows × {raw.shape[1]} columns")

print("\nFirst 5 rows:")
print(raw.head(5).to_string())

print("\nUnique categories:")
for cat, n in raw["Category"].value_counts().items():
    print(f"  {cat}: {n} metabolites")

print("\n*** CRITICAL FINDING ***")
print("This CSV contains SUMMARY STATISTICS (group means, SDs, pre-computed")
print("p-values), NOT individual patient measurements.")
print("Reported sample sizes (Control n=25, NAFL n=42, NASH n=19) are stated")
print("in the paper; they cannot be verified from this file alone.")
print("Steps requiring individual-level data are clearly flagged below.")


# ═══════════════════════════════════════════════════════════════════════════
# STEP 2 — QC
# ═══════════════════════════════════════════════════════════════════════════
print("\n" + "="*70)
print("STEP 2 — QUALITY CONTROL")
print("="*70)

df = raw.copy()

# ── 2a. Missingness ─────────────────────────────────────────────────────────
numeric_cols = ["Control_Mean","Control_SD","NAFL_Mean","NAFL_SD",
                "NASH_Mean","NASH_SD","Normalized_NAFL","Normalized_NASH"]
print("\nMissing values per summary column:")
for c in numeric_cols:
    n_miss = df[c].isna().sum()
    print(f"  {c}: {n_miss} missing")

# ── 2b. Group sizes from the paper ──────────────────────────────────────────
N_CTRL, N_NAFL, N_NASH = 25, 42, 19
N_TOTAL = N_CTRL + N_NAFL + N_NASH
print(f"\nGroup sizes per paper: Control={N_CTRL}, NAFL={N_NAFL}, NASH={N_NASH} → total={N_TOTAL}")
print("(Cannot verify independently from summary-statistics table)")

# ── 2c. Aspartic acid outlier ────────────────────────────────────────────────
asp = df[df["Metabolite"] == "Aspartic acid"].iloc[0]
print(f"\nAspartic acid summary:")
print(f"  Control  mean={asp['Control_Mean']:.2f} SD={asp['Control_SD']:.2f}")
print(f"  NAFL     mean={asp['NAFL_Mean']:.2f}  SD={asp['NAFL_SD']:.2f}")
print(f"  NASH     mean={asp['NASH_Mean']:.2f} SD={asp['NASH_SD']:.2f}")

# Infer outlier: the paper flagged one NASH sample ~382.4 ug/uL
# With n=19 NASH, mean=43.18, SD=80.63:
# If one sample = X, mean of rest = (43.18*19 - X)/18
# SD inflation of ~80.63 with n=19 is consistent with one extreme outlier
inferred_x = asp["NASH_Mean"] * N_NASH  # total sum
print(f"\n  NASH SD={asp['NASH_SD']:.2f} with mean={asp['NASH_Mean']:.2f} is highly inflated.")
print(f"  This is consistent with the outlier (~382.4 µg/uL) the paper itself")
print(f"  described in one NASH patient.")
print(f"  DECISION: RETAIN the outlier (paper's own decision; Kruskal-Wallis")
print(f"  is rank-based and robust to this outlier).")

# ── 2d. Near-zero or negative means (log-transform risk) ────────────────────
tiny_means = df[numeric_cols[:6]].min(axis=1)
n_zero = (tiny_means <= 0).sum()
n_tiny = (tiny_means < 0.001).sum()
print(f"\nMetabolites with any group mean ≤ 0 : {n_zero}")
print(f"Metabolites with any group mean < 0.001: {n_tiny}")
print("Smallest group means:")
df["min_mean"] = df[["Control_Mean","NAFL_Mean","NASH_Mean"]].min(axis=1)
print(df.nsmallest(5, "min_mean")[["Metabolite","Control_Mean","NAFL_Mean","NASH_Mean"]].to_string(index=False))
print("Note: All group means are positive — log-transform is safe at the")
print("summary-statistics level, though individual values could be near zero.")


# ── Parse p-value columns (handle '<0.001' strings) ─────────────────────────
def parse_pval(s):
    s = str(s).strip()
    if s.startswith("<"):
        return 0.001   # conservative upper bound
    try:
        return float(s)
    except ValueError:
        return np.nan

for col in ["Pval_Control_vs_NAFL","Pval_Control_vs_NASH","Pval_NAFL_vs_NASH",
            "Pval_Kruskal_Wallis","Qval_FDR"]:
    df[col] = df[col].apply(parse_pval)


# ═══════════════════════════════════════════════════════════════════════════
# STEP 3 — NORMALITY CHECK
# ═══════════════════════════════════════════════════════════════════════════
print("\n" + "="*70)
print("STEP 3 — NORMALITY CHECK (SHAPIRO-WILK)")
print("="*70)
print("""
INDIVIDUAL DATA REQUIRED — NOT POSSIBLE FROM THIS FILE.

Shapiro-Wilk tests each group of raw measurements for each metabolite.
Because this CSV contains only means and SDs, the raw values are unavailable
and the test cannot be run.

What would it show?
  Metabolomics data are routinely right-skewed (many low-abundance analytes,
  occasional extreme values like the Aspartic acid outlier).  A log-transform
  usually normalises the data, but Shapiro-Wilk would confirm this formally.
  Even one non-normal metabolite in any group justifies using Kruskal-Wallis
  (non-parametric) instead of ANOVA.  The paper chose non-parametric tests,
  which is consistent with the expectation of non-normality.

Conclusion: non-parametric tests are justified by convention and paper choice.
""")


# ═══════════════════════════════════════════════════════════════════════════
# STEP 4 — UNIVARIATE STATISTICS
# ═══════════════════════════════════════════════════════════════════════════
print("\n" + "="*70)
print("STEP 4 — UNIVARIATE STATISTICS (KRUSKAL-WALLIS + BH-FDR)")
print("="*70)

# Re-derive FDR from the KW p-values in the table.
# NOTE: "<0.001" values were converted to 0.001 (upper bound), so our
# re-computed FDR may differ slightly from the paper's for borderline metabolites.
# We therefore report BOTH: our re-computed FDR and the paper's Qval_FDR.
kw_pvals = df["Pval_Kruskal_Wallis"].values
reject, fdr_recomputed, _, _ = multipletests(kw_pvals, method="fdr_bh")
df["FDR_recomputed"] = fdr_recomputed

# Use paper's Qval_FDR as primary reference (avoids the <0.001 rounding issue)
df["FDR_paper"] = df["Qval_FDR"]   # already parsed above

print(f"\nKruskal-Wallis p-values sourced from paper's Table S1.")
print(f"BH-FDR re-computed independently using statsmodels multipletests().")
print(f"Note: '<0.001' KW p-values were approximated as 0.001 for re-computation.")
print(f"This flattens ranks, so we cross-check with the paper's Qval_FDR column.")
print(f"\nMetabolites significant after FDR correction (q < 0.05):")

# Use paper's Qval_FDR as the primary significance filter
sig = df[df["FDR_paper"] < 0.05].sort_values("FDR_paper")
for _, row in sig.iterrows():
    trend = ""
    if row["Control_Mean"] <= row["NAFL_Mean"] <= row["NASH_Mean"]:
        trend = "↑ Ctrl < NAFL < NASH"
    elif row["Control_Mean"] >= row["NAFL_Mean"] >= row["NASH_Mean"]:
        trend = "↓ Ctrl > NAFL > NASH"
    else:
        trend = "non-monotone"
    print(f"  {row['Metabolite']:<30s}  KW p={row['Pval_Kruskal_Wallis']:.4f}  "
          f"FDR q(paper)={row['FDR_paper']:.4f}  FDR q(recomp)={row['FDR_recomputed']:.4f}  {trend}")

print(f"\nTotal FDR-significant metabolites (paper's Qval_FDR < 0.05): {(df['FDR_paper'] < 0.05).sum()}")

# Paper's stated list: glutamic acid, tyrosine, kynurenic acid,
# alpha-ketoglutaric acid, myristoleic acid, palmitoleic acid
PAPER_SIG = {"Glutamic acid","Tyrosine","Kynurenic acid",
             "a-Ketoglutaric acid","Myristoleic acid","Palmitoleic acid"}
our_sig = set(sig["Metabolite"].tolist())
print(f"\nSanity check vs paper's 6 reported significant metabolites:")
for m in sorted(PAPER_SIG):
    match = any(m.lower() in x.lower() for x in our_sig)
    print(f"  {'[MATCH]' if match else '[MISS ]'}  {m}")

# Pairwise Wilcoxon (from paper's pre-calculated values, already BH-corrected)
print("\nPairwise Wilcoxon p-values for FDR-significant metabolites")
print("(sourced from paper's Table S1 — already reported):")
cols = ["Metabolite","Pval_Control_vs_NAFL","Pval_Control_vs_NASH","Pval_NAFL_vs_NASH"]
print(df[df["FDR_recomputed"] < 0.05][cols].to_string(index=False))

# Save results
df.to_csv(os.path.join(RESULTS_DIR, "all_metabolites_statistics.csv"), index=False)
sig.to_csv(os.path.join(RESULTS_DIR, "significant_metabolites_FDR05.csv"), index=False)
print(f"\nFull statistics saved → results/all_metabolites_statistics.csv")
print(f"Significant metabolites saved → results/significant_metabolites_FDR05.csv")


# ═══════════════════════════════════════════════════════════════════════════
# STEP 5 — TREND TEST (JONCKHEERE-TERPSTRA)
# Enhancement — not in the original paper; clearly labelled as such
# ═══════════════════════════════════════════════════════════════════════════
print("\n" + "="*70)
print("STEP 5 — TREND TEST (ENHANCEMENT — NOT IN ORIGINAL PAPER)")
print("="*70)

print("""
INDIVIDUAL DATA REQUIRED for a formal Jonckheere-Terpstra test.

What JT does: unlike Kruskal-Wallis, which only asks "are the groups
different?", JT asks "do the groups follow a specific ordered direction?".
Here the biology predicts Control < NAFL < NASH (worsening liver disease).
JT is more powerful than KW when this monotone direction is expected.

APPROXIMATION FROM SUMMARY DATA:
For each metabolite, we compute Spearman correlation between group order
(Control=0, NAFL=1, NASH=2) and group means.  This captures direction
but not significance (no individual variance or sample-size information).
It is labelled 'Trend direction (proxy)' and flagged as approximate.
""")

group_order = np.array([0, 1, 2])

def trend_direction(row):
    means = np.array([row["Control_Mean"], row["NAFL_Mean"], row["NASH_Mean"]])
    rho, _ = stats.spearmanr(group_order, means)
    return rho

df["Trend_rho"] = df.apply(trend_direction, axis=1)

# Refresh sig now that Trend_rho is in df
sig = df[df["FDR_paper"] < 0.05].sort_values("FDR_paper").copy()

monotone_up   = ((df["Control_Mean"] <= df["NAFL_Mean"]) &
                 (df["NAFL_Mean"]    <= df["NASH_Mean"])).sum()
monotone_down = ((df["Control_Mean"] >= df["NAFL_Mean"]) &
                 (df["NAFL_Mean"]    >= df["NASH_Mean"])).sum()
print(f"Metabolites with monotone ↑ trend (Ctrl ≤ NAFL ≤ NASH): {monotone_up}")
print(f"Metabolites with monotone ↓ trend (Ctrl ≥ NAFL ≥ NASH): {monotone_down}")

print("\nFDR-significant metabolites — trend direction:")
for _, row in sig.iterrows():
    print(f"  {row['Metabolite']:<30s}  Spearman ρ = {row['Trend_rho']:+.3f}")

trend_results = df[["Metabolite","Category","Trend_rho","FDR_recomputed","FDR_paper"]].copy()
trend_results.to_csv(os.path.join(RESULTS_DIR, "trend_direction_proxy.csv"), index=False)
print("\nTrend direction (proxy) saved → results/trend_direction_proxy.csv")


# ═══════════════════════════════════════════════════════════════════════════
# STEP 6 — VISUALIZATION
# ═══════════════════════════════════════════════════════════════════════════
print("\n" + "="*70)
print("STEP 6 — VISUALIZATION")
print("="*70)

# ── Color palette ────────────────────────────────────────────────────────────
GROUP_COLORS = {"Control": "#4E9AC7", "NAFL": "#F5A623", "NASH": "#D0021B"}
GROUPS = ["Control", "NAFL", "NASH"]

# ── 6a. Bar plots for 6 FDR-significant metabolites ─────────────────────────
# (True boxplots require individual patient data — unavailable in this CSV)
sig_mets = sig["Metabolite"].tolist()

fig, axes = plt.subplots(2, 3, figsize=(14, 9))
axes = axes.flatten()

for ax, met in zip(axes, sig_mets):
    row = df[df["Metabolite"] == met].iloc[0]
    means = [row["Control_Mean"], row["NAFL_Mean"], row["NASH_Mean"]]
    sds   = [row["Control_SD"],   row["NAFL_SD"],   row["NASH_SD"]]

    x = np.arange(3)
    bars = ax.bar(x, means, yerr=sds, capsize=6,
                  color=[GROUP_COLORS[g] for g in GROUPS],
                  edgecolor="black", linewidth=0.8, width=0.6)
    ax.set_xticks(x)
    ax.set_xticklabels(GROUPS, fontsize=10)
    ax.set_title(met, fontsize=11, fontweight="bold")
    ax.set_ylabel("Mean ± SD (µg/uL)", fontsize=9)
    q = row["FDR_recomputed"]
    ax.text(0.98, 0.97, f"FDR q = {q:.4f}", transform=ax.transAxes,
            ha="right", va="top", fontsize=9, color="#555555")

fig.suptitle("FDR-Significant Metabolites — Ji et al. 2022\n"
             "(mean ± SD; individual-level data not available for boxplots)",
             fontsize=13, fontweight="bold", y=1.01)
plt.tight_layout()
out = os.path.join(PLOTS_DIR, "significant_metabolites_bar.png")
plt.savefig(out, dpi=150, bbox_inches="tight")
plt.close()
print(f"\nBar plots saved → plots/significant_metabolites_bar.png")

# ── 6b. Z-score heatmap (all 79 metabolites) ────────────────────────────────
# Z-score computed across the 3 group means per metabolite.
# This is a proxy heatmap (from summary data), clearly labelled.

hm_data = df[["Metabolite","Category","Control_Mean","NAFL_Mean","NASH_Mean"]].copy()
hm_data = hm_data.set_index("Metabolite")

# Compute pseudo-Z-score across the 3 group means per row
mean_matrix = hm_data[["Control_Mean","NAFL_Mean","NASH_Mean"]].values.astype(float)
row_means = mean_matrix.mean(axis=1, keepdims=True)
row_stds  = mean_matrix.std(axis=1, keepdims=True)
row_stds[row_stds == 0] = 1   # avoid divide-by-zero
z_matrix = (mean_matrix - row_means) / row_stds

z_df = pd.DataFrame(z_matrix,
                    index=hm_data.index,
                    columns=["Control", "NAFL", "NASH"])
z_df["Category"] = hm_data["Category"].values

# Sort by category then by NASH-Control direction
z_df["sort_key"] = z_df["NASH"] - z_df["Control"]
z_df = z_df.sort_values(["Category", "sort_key"], ascending=[True, False])

# Category row-colour strip
cat_palette = {
    "Amino acids": "#6BAED6",
    "Kynurenine pathway metabolites and nucleosides": "#74C476",
    "Organic acids": "#FD8D3C",
    "Fatty acids": "#9E9AC8"
}

fig, ax = plt.subplots(figsize=(8, 18))
hm_mat = z_df[["Control","NAFL","NASH"]].values

im = ax.imshow(hm_mat, aspect="auto", cmap="RdBu_r", vmin=-1.5, vmax=1.5)

ax.set_xticks([0, 1, 2])
ax.set_xticklabels(["Control", "NAFL", "NASH"], fontsize=11, fontweight="bold")
ax.set_yticks(range(len(z_df)))
ax.set_yticklabels(z_df.index, fontsize=7)

# Mark FDR-significant metabolites
for i, met in enumerate(z_df.index):
    if met in our_sig:
        ax.annotate("*", xy=(2.55, i), fontsize=10, color="red",
                    ha="left", va="center", fontweight="bold")

plt.colorbar(im, ax=ax, label="Pseudo Z-score (group means)", shrink=0.4)
ax.set_title("Z-Score Heatmap — All 79 Metabolites\n"
             "Ji et al. 2022  ·  * = FDR < 0.05\n"
             "(computed from group means — proxy for Figure 1A)",
             fontsize=11, fontweight="bold")

# Category legend
patches = [mpatches.Patch(color=v, label=k) for k, v in cat_palette.items()]
ax.legend(handles=patches, loc="lower right", fontsize=7,
          title="Category", bbox_to_anchor=(1.35, 0.0))

plt.tight_layout()
out = os.path.join(PLOTS_DIR, "heatmap_79metabolites_zscore.png")
plt.savefig(out, dpi=150, bbox_inches="tight")
plt.close()
print(f"Z-score heatmap saved → plots/heatmap_79metabolites_zscore.png")

# ── 6c. Volcano-style plot: -log10(FDR) vs fold change (NASH / Control) ─────
df["log2FC_NASH_Ctrl"] = np.log2(df["NASH_Mean"] / df["Control_Mean"])
df["neg_log10_FDR"]    = -np.log10(df["FDR_recomputed"])

fig, ax = plt.subplots(figsize=(9, 6))

cat_colors_map = {
    "Amino acids": "#6BAED6",
    "Kynurenine pathway metabolites and nucleosides": "#74C476",
    "Organic acids": "#FD8D3C",
    "Fatty acids": "#9E9AC8"
}

for cat, grp in df.groupby("Category"):
    sig_mask = grp["FDR_recomputed"] < 0.05
    color = cat_colors_map.get(cat, "grey")
    ax.scatter(grp.loc[~sig_mask, "log2FC_NASH_Ctrl"],
               grp.loc[~sig_mask, "neg_log10_FDR"],
               color=color, alpha=0.6, s=40, label=cat, edgecolors="none")
    ax.scatter(grp.loc[sig_mask, "log2FC_NASH_Ctrl"],
               grp.loc[sig_mask, "neg_log10_FDR"],
               color=color, alpha=1.0, s=80, edgecolors="black", linewidths=0.8)

# Label significant metabolites
for _, row in df[df["FDR_recomputed"] < 0.05].iterrows():
    ax.annotate(row["Metabolite"], (row["log2FC_NASH_Ctrl"], row["neg_log10_FDR"]),
                textcoords="offset points", xytext=(5, 3), fontsize=8)

ax.axhline(-np.log10(0.05), color="red", linestyle="--", linewidth=1, label="FDR = 0.05")
ax.axvline(0, color="grey", linestyle=":", linewidth=0.8)
ax.set_xlabel("log₂ Fold Change (NASH / Control)", fontsize=11)
ax.set_ylabel("-log₁₀(FDR q-value)", fontsize=11)
ax.set_title("NASH vs Control — Metabolite Significance\nJi et al. 2022", fontsize=12)
ax.legend(fontsize=8, loc="upper left", ncol=1)
plt.tight_layout()
out = os.path.join(PLOTS_DIR, "volcano_NASH_vs_Control.png")
plt.savefig(out, dpi=150, bbox_inches="tight")
plt.close()
print(f"Volcano plot saved → plots/volcano_NASH_vs_Control.png")


# ═══════════════════════════════════════════════════════════════════════════
# STEP 7 — MULTIVARIATE ANALYSIS (PLS-DA)
# ═══════════════════════════════════════════════════════════════════════════
print("\n" + "="*70)
print("STEP 7 — MULTIVARIATE ANALYSIS (PLS-DA)")
print("="*70)

print("""
INDIVIDUAL DATA REQUIRED — NOT POSSIBLE FROM THIS FILE.

PLS-DA (Partial Least Squares Discriminant Analysis) builds a model that
separates groups using all 79 metabolites simultaneously.  It needs a
matrix of shape (86 subjects × 79 metabolites) — one row per patient.
The paper log-transformed and autoscaled (mean-centred, divided by SD)
each metabolite column, then fit a 2-component PLS-DA.

What it would show:
  - A scores plot (2D scatter, one point per patient) revealing whether
    Control / NAFL / NASH form distinct clusters.
  - VIP (Variable Importance in Projection) scores per metabolite.
    VIP > 1.0 flags the metabolites that drive the separation.  The paper
    found the same 6 FDR-significant metabolites all had VIP > 1.0.

Why not simulate?
  Simulating 86 subjects from group means and SDs would ignore the
  covariance structure (correlations between metabolites), producing
  a misleading PLS-DA that cannot replicate the paper's results.

Action: individual data can be requested from the corresponding author
(contact info in the paper).  Once obtained, run scikit-learn PLSRegression
or the mixOmics R package as the paper used.
""")


# ═══════════════════════════════════════════════════════════════════════════
# STEP 8 — PATHWAY ENRICHMENT
# ═══════════════════════════════════════════════════════════════════════════
print("\n" + "="*70)
print("STEP 8 — PATHWAY ENRICHMENT")
print("="*70)

# Export significant metabolite names for MetaboAnalyst
sig_names = sig["Metabolite"].tolist()
print("\nFDR-significant metabolites for MetaboAnalyst input:")
for m in sig_names:
    print(f"  {m}")

with open(os.path.join(RESULTS_DIR, "sig_metabolites_for_metaboanalyst.txt"), "w") as f:
    f.write("\n".join(sig_names) + "\n")

# Export all metabolites ranked by FDR for over-representation analysis
df_export = df[["Metabolite","Category","FDR_recomputed","log2FC_NASH_Ctrl"]].copy()
df_export = df_export.sort_values("FDR_recomputed")
df_export.to_csv(os.path.join(RESULTS_DIR, "all_metabolites_ranked.csv"), index=False)

print(f"""
Pathway enrichment cannot be run inline without a maintained KEGG/HMDB
mapping (which changes over time — hand-rolled mappings go stale).

How to proceed:
  1. Go to: https://www.metaboanalyst.ca/MetaboAnalyst/upload/EnrichUploadView.xhtml
  2. Upload: results/sig_metabolites_for_metaboanalyst.txt  (6 names)
  3. Select 'Pathway Analysis' → KEGG → Homo sapiens
  4. The paper found enrichment in alanine/aspartate/glutamate metabolism
     and aminoacyl-tRNA biosynthesis (driven mainly by glutamic acid /
     alpha-ketoglutaric acid / the two monounsaturated fatty acids).

For a broader background set, upload:
  results/all_metabolites_ranked.csv  (79 metabolites, ranked by FDR)

Exported files:
  results/sig_metabolites_for_metaboanalyst.txt
  results/all_metabolites_ranked.csv
""")


# ═══════════════════════════════════════════════════════════════════════════
# SUMMARY
# ═══════════════════════════════════════════════════════════════════════════
print("\n" + "="*70)
print("ANALYSIS COMPLETE — FILE SUMMARY")
print("="*70)
print(f"  plots/significant_metabolites_bar.png")
print(f"  plots/heatmap_79metabolites_zscore.png")
print(f"  plots/volcano_NASH_vs_Control.png")
print(f"  results/all_metabolites_statistics.csv")
print(f"  results/significant_metabolites_FDR05.csv")
print(f"  results/trend_direction_proxy.csv")
print(f"  results/sig_metabolites_for_metaboanalyst.txt")
print(f"  results/all_metabolites_ranked.csv")
