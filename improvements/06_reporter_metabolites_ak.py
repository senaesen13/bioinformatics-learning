#!/usr/bin/env python3
"""
==============================================================================
Script: 06_reporter_metabolites_ak.py
Author: Dr. Ali Kaynar (King's College London)
Description: Patil & Nielsen Reporter Metabolite Algorithm with GEM Model Input (Python)
==============================================================================
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
from scipy.stats import norm

def get_root_dir():
    script_path = os.path.abspath(__file__)
    return os.path.dirname(os.path.dirname(script_path))

root_dir = get_root_dir()
out_dir = os.path.join(root_dir, "improvements", "results_06_reporter_metabolites_python")
os.makedirs(out_dir, exist_ok=True)

parser = argparse.ArgumentParser(description="GEM-based Reporter Metabolite Analysis")
parser.add_argument("--gem", type=str, default=None, help="Path to GEM model file (.xml, .csv, .yml)")
args, unknown = parser.parse_known_args()

print("============================================================")
print("AK Pipeline 06 (Python): GEM-Based Reporter Metabolite Analysis")
print("Root Directory:", root_dir)
print("Output Directory:", out_dir)
print("============================================================\n")

# 1. Resolve GEM Input Path
gem_path = args.gem
if not gem_path:
    local_xml = os.path.join(root_dir, "improvements", "data", "Human-GEM.xml")
    csv_default = os.path.join(root_dir, "improvements", "data", "human_gem_topology_network.csv")
    if os.path.exists(local_xml):
        gem_path = local_xml
    elif os.path.exists(csv_default):
        gem_path = csv_default
    else:
        gem_path = "/Users/k2254978/Desktop/Work/ToolBoxProgs/Human_GEM_2.0/model/Human-GEM.xml"

print(f"[INFO] GEM Model Input File: {gem_path}")

if gem_path.endswith(".xml"):
    from parse_gem import parse_human_gem_xml
    csv_path = os.path.join(root_dir, "improvements", "data", "human_gem_topology_network.csv")
    gem_df = parse_human_gem_xml(gem_path, csv_path)
elif gem_path.endswith(".csv"):
    gem_df = pd.read_csv(gem_path)
else:
    raise ValueError(f"Unsupported GEM model format: {gem_path}")

# 2. Currency Metabolites Filter
currency_metabolites = {"h2o", "atp", "adp", "amp", "nad+", "nadh", "nadp+", "nadph",
                        "h+", "pi", "ppi", "coa", "co2", "o2", "hco3-", "na+", "k+", "cl-", "water", "oxygen", "phosphate"}

gem_df_clean = gem_df.dropna(subset=["metabolite_name", "gene_symbol"]).copy()
gem_df_clean = gem_df_clean[~gem_df_clean["metabolite_name"].str.lower().isin(currency_metabolites)]

print(f"Loaded {len(gem_df_clean)} GEM network associations across {gem_df_clean['metabolite_name'].nunique()} metabolites.")

# 3. Load Differential Expression Data
res_file = os.path.join(root_dir, "week4-day3-nafld-gse130970-rnaseq", "results", "gse130970_results.csv")
if not os.path.exists(res_file):
    res_file = os.path.join(root_dir, "week4-day1-nafld-bulk-rnaseq", "results", "deseq2_results.csv")

def find_col(df_cols, candidates):
    for c in candidates:
        if c in df_cols:
            return c
    return None

if os.path.exists(res_file):
    print(f"[INFO] Loading differential gene expression results from: {res_file}")
    df_raw = pd.read_csv(res_file)
    cols = list(df_raw.columns)
    
    sym_col = find_col(cols, ["gene_symbol", "symbol", "Gene", "gene"])
    pval_col = find_col(cols, ["pvalue_mle", "pvalue", "p.value", "pval"])
    lfc_col = find_col(cols, ["lfc_mle", "log2FoldChange", "lfc", "lfc_apeglm"])
    
    df_clean = df_raw.dropna(subset=[sym_col, pval_col]).copy()
    df_clean["gene_symbol"] = df_clean[sym_col].astype(str)
    df_clean["pvalue"] = np.maximum(df_clean[pval_col].values, 1e-300)
    df_clean["lfc"] = df_clean[lfc_col].values if lfc_col else 0.0
    
    df_clean = df_clean.sort_values(by="pvalue").groupby("gene_symbol").first().reset_index()
else:
    print("[INFO] Generating synthetic gene dataset...")
    all_genes = list(set(gem_df_clean["gene_symbol"].tolist() + [f"GENE_{i:04d}" for i in range(1, 501)]))
    np.random.seed(42)
    df_clean = pd.DataFrame({
        "gene_symbol": all_genes,
        "pvalue": np.random.uniform(1e-5, 0.5, size=len(all_genes)),
        "lfc": np.random.normal(0, 1.5, size=len(all_genes))
    })

# 4. Calculate Gene Z-scores
df_clean["gene_z"] = norm.ppf(1.0 - df_clean["pvalue"].values / 2.0)
gene_z_dict = dict(zip(df_clean["gene_symbol"], df_clean["gene_z"]))
gene_lfc_dict = dict(zip(df_clean["gene_symbol"], df_clean["lfc"]))
gene_universe = list(gene_z_dict.keys())

# 5. Patil & Nielsen Reporter Algorithm across GEM Nodes
met_groups = gem_df_clean.groupby("metabolite_name")["gene_symbol"].unique().to_dict()
qual_mets = {m: genes for m, genes in met_groups.items() if 3 <= len([g for g in genes if g in gene_z_dict]) <= 100}

print(f"Running Reporter Metabolite Analysis across {len(qual_mets)} qualified GEM nodes...")

n_perm = 1000
results = []

for idx, (met, genes) in enumerate(qual_mets.items()):
    g_set = [g for g in genes if g in gene_z_dict]
    k = len(g_set)
    
    z_scores = [gene_z_dict[g] for g in g_set]
    lfcs = [gene_lfc_dict[g] for g in g_set]
    
    z_raw = np.sum(z_scores) / np.sqrt(k)
    avg_lfc = np.mean(lfcs)
    
    np.random.seed(42 + idx)
    perm_z = [np.sum(np.random.choice(list(gene_z_dict.values()), size=k, replace=False)) / np.sqrt(k) for _ in range(n_perm)]
    
    mu_k = np.mean(perm_z)
    sigma_k = np.std(perm_z)
    z_corr = (z_raw - mu_k) / sigma_k if sigma_k > 0 else z_raw
    p_rep = 1.0 - norm.cdf(z_corr)
    
    results.append({
        "Metabolite": met,
        "Neighbor_Genes_Count": k,
        "Mean_Log2FC": round(avg_lfc, 3),
        "Z_Raw": round(z_raw, 3),
        "Reporter_Z_Score": round(z_corr, 3),
        "Pvalue": float(f"{p_rep:.4e}"),
        "Neighbor_Genes": "; ".join(g_set[:8])
    })

res_df = pd.DataFrame(results).sort_values(by="Reporter_Z_Score", ascending=False)
print(res_df.head(15).to_string(index=False))

# 6. Save CSV & PNG Plot
out_csv = os.path.join(out_dir, "gem_reporter_metabolites_python.csv")
res_df.to_csv(out_csv, index=False)
print(f"\nSaved CSV: {out_csv}")

try:
    import matplotlib.pyplot as plt
    
    top_df = res_df.head(20)
    plt.figure(figsize=(8, 6.5))
    colors = ["#D81B60" if lfc > 0 else "#1E88E5" for lfc in top_df["Mean_Log2FC"]]
    plt.barh(top_df["Metabolite"], top_df["Reporter_Z_Score"], color=colors, edgecolor="black")
    plt.xlabel("Corrected Reporter Z-score", fontsize=11, fontweight="bold")
    plt.ylabel("Metabolite Node", fontsize=11, fontweight="bold")
    plt.title("Top 20 GEM Reporter Metabolites (Human-GEM Input - Python)", fontsize=12, fontweight="bold")
    plt.gca().invert_yaxis()
    plt.tight_layout()
    
    out_png = os.path.join(out_dir, "gem_reporter_metabolites_python.png")
    plt.savefig(out_png, dpi=300)
    plt.close()
    print(f"Saved PNG Plot: {out_png}")
except Exception as e:
    print(f"[NOTE] Matplotlib plot skipped: {e}")

print("[SUCCESS] Python Pipeline 06 completed successfully using GEM model input!")
