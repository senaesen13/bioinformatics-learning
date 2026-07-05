# Week 4 Day 2 — NAFLD RNA-seq Validation Cohort

## Dataset: GSE135251

| Field | Value |
|---|---|
| GEO Accession | GSE135251 |
| Publication | Govaere et al. (2020), *Science Translational Medicine* |
| Platform | Illumina HiSeq (GPL18573) |
| Tissue | Liver biopsy (human) |
| Total samples | 216 |
| Normal (healthy controls) | 10 |
| NAFLD (all stages NAFL + NASH F0–F4) | 206 |
| Alignment | GRCh38; raw HTSeq counts, one file per sample |
| Gene identifiers | Ensembl IDs (no version suffixes) |

**Sample breakdown (as in paper):**
control=10 · NAFL=51 · NASH_F0-F1=34 · NASH_F2=53 · NASH_F3=54 · NASH_F4=14

This dataset is the independent validation cohort for the discovery analysis in
`week4-day1-nafld-bulk-rnaseq/` (GSE162694, Suppli et al. 2021).

---

## Grouping Strategy

All NAFL and NASH samples (regardless of fibrosis stage F0–F4) are collapsed into
a single **"NAFLD"** group and compared against healthy controls (**"Normal"**).
This matches the binary Normal vs NAFLD design used in the discovery cohort.

The disease status column in GEO metadata is `disease:ch1`:
- `"disease: Control"` → Normal (n = 10)
- `"disease: NAFLD"` → NAFLD (n = 206)

---

## Methodology

Identical pipeline to the discovery cohort:

1. **Protein-coding filter** — applied before count filtering. Same gene list derived
   from biomaRt in the discovery cohort (`protein_coding_gene_list.rds`).
   Keeps only `gene_biotype == protein_coding` Ensembl IDs.

2. **Mean-count filter** — `rowMeans(counts) >= 10`.

3. **DESeq2** — design `~ condition` (reference = Normal). Both MLE (un-shrunk) and
   apeglm-shrunk LFC are computed and stored in the results table.

4. **Significance threshold** — `padj < 0.01` AND `|log2FC_MLE| > 2`.
   MLE LFC is used for threshold filtering; apeglm LFC is reported in the output table
   (consistent with the discovery cohort approach).

### Filtering Counts

| Stage | Genes |
|---|---|
| Starting count matrix (HTSeq GRCh38) | 64,253 |
| After protein-coding filter | 16,390 |
| After mean-count filter (rowMeans ≥ 10) | **13,166** |

---

## Significant Gene Results

**Threshold:** padj < 0.01 AND |MLE log2FC| > 2

| Direction | Count |
|---|---|
| Up in NAFLD | 63 |
| Down in NAFLD | 40 |
| **Total** | **103** |

The validation cohort yields more significant genes (103 vs 24 in discovery) despite
having only 10 controls, because the NAFLD group here is enriched for advanced fibrosis
(F3+F4 = 68/206 samples), which increases average effect size relative to controls.

---

## TREM2 / SPP1 / GPNMB

| Gene | MLE log2FC | padj | In sig list? |
|---|---|---|---|
| **TREM2** | +2.52 | 2.2e-07 | **YES** |
| SPP1 | +1.41 | 2.8e-03 | NO — LFC below 2× |
| GPNMB | +0.91 | 6.0e-03 | NO — LFC below 2× |

**TREM2** replicates as significant with concordant direction (upregulated in NAFLD
in both cohorts, similar magnitude: +2.48 discovery, +2.52 validation).

**SPP1 and GPNMB** are statistically significant (padj < 0.01 in discovery; padj < 0.05
in validation) and consistently upregulated, but their fold changes stay below 2× in
both cohorts. This is expected for bulk RNA-seq across mixed fibrosis stages — both
genes are expressed in macrophage/stellate cell subpopulations that are diluted within
whole-liver biopsies. They are genuine NAFLD markers but the bulk RNA-seq fold change
does not reach the stringent |LFC| > 2 threshold.

---

## Cross-Cohort Validation Summary

**7 genes significant in both cohorts** (Ensembl ID match):

| Gene | disc log2FC | disc padj | val log2FC | val padj | Concordant |
|---|---|---|---|---|---|
| TREM2 | +2.48 | 5.2e-13 | +2.52 | 2.2e-07 | YES |
| CHI3L1 | +2.23 | 1.3e-09 | +2.19 | 1.1e-04 | YES |
| STMN2 | +2.68 | 4.0e-08 | +2.57 | 3.6e-03 | YES |
| PDZK1IP1 | +2.20 | 1.4e-07 | +2.23 | 2.2e-04 | YES |
| PRAMEF10 | +3.52 | 2.5e-10 | +3.23 | 3.2e-05 | YES |
| NR4A1 | +2.04 | 8.5e-08 | −3.72 | 6.4e-35 | NO |
| RGS1 | +2.04 | 1.2e-06 | −4.52 | 1.6e-49 | NO |

NR4A1 and RGS1 show opposite directions between cohorts (up in discovery, strongly
down in validation). Both regulate immune cell function; the discordance likely reflects
cohort-level differences in inflammatory cell composition. The large and significant
downregulation in the validation cohort is noteworthy and may reflect hepatocyte
dedifferentiation effects in advanced fibrosis that dominate in GSE135251.

| Metric | Value |
|---|---|
| Discovery sig genes | 24 |
| Validation sig genes | 103 |
| Overlap | 7 genes |
| Overlap % of discovery | 29.2% |
| Fisher's exact OR | 55.52 |
| Fisher's exact p-value | 4.7e-10 |
| Direction concordance | 71.4% (5/7) |
| Spearman ρ (sig overlap, n=7) | 0.929 (p = 0.007) |
| Pearson r (sig overlap, n=7) | 0.624 (p = 0.134 — underpowered at n=7) |
| Pearson r (genome-wide, 13,100 genes) | 0.323 |

The overlap is highly enriched above chance (OR = 55.5, p = 4.7e-10). The low raw
overlap count (7 genes) is a consequence of the stringent |LFC| > 2 threshold
producing small gene lists. At the genome-wide level (13,100 shared tested genes),
the two cohorts show moderate positive correlation (Pearson r = 0.32), consistent
with independent patient populations studying the same disease but with different
fibrosis stage distributions.

---

## GSEA Results (Validation Cohort)

**Method:** Genes ranked by apeglm-shrunken log2FC (NAFLD / Normal), descending.
clusterProfiler `gseKEGG` (organism = "hsa") + `GSEA` with MSigDB Hallmark (H collection).
Both use BH-adjusted p-value cutoff < 0.05, minGSSize = 15, maxGSSize = 500.

**KEGG:** 87 significant pathways

| Direction | Top pathways (padj) |
|---|---|
| Activated in NAFLD | ECM-receptor interaction (4.4e-04), Integrin signaling (1.6e-03), Cornified envelope formation (3.7e-03), Focal adhesion (3.9e-02), Cholesterol metabolism (3.0e-02) |
| Suppressed in NAFLD | Ribosome (3.5e-08), Coronavirus disease (2.3e-06), IL-17 signaling (2.3e-04), Amphetamine addiction (5.3e-04), MAPK signaling (5.3e-04) |

Note: Several immune/inflammatory pathways (IL-17, TNF, osteoclast differentiation) are
**suppressed** in this cohort, opposite to the discovery cohort — see cross-cohort comparison below.

**Hallmark:** 18 significant gene sets

| Direction | Gene set | NES | padj |
|---|---|---|---|
| Activated | CHOLESTEROL_HOMEOSTASIS | +2.12 | 1.6e-05 |
| Activated | APICAL_JUNCTION | +1.94 | 9.2e-06 |
| Activated | MYOGENESIS | +1.68 | 3.3e-03 |
| Activated | MTORC1_SIGNALING | +1.58 | 6.9e-03 |
| Activated | MITOTIC_SPINDLE | +1.57 | 8.9e-03 |
| Suppressed | TNFA_SIGNALING_VIA_NFKB | −3.04 | 5.0e-09 |
| Suppressed | HYPOXIA | −2.13 | 8.8e-07 |
| Suppressed | KRAS_SIGNALING_UP | −2.00 | 9.1e-06 |
| Suppressed | UV_RESPONSE_UP | −1.75 | 1.9e-03 |
| Suppressed | TGF_BETA_SIGNALING | −1.76 | 9.3e-03 |

---

## Cross-Cohort Pathway Comparison (GSEA)

Discovery = GSE162694 (Day 1) · Validation = GSE135251 (Day 2)
Overlap = pathways/gene sets significant (padj < 0.05) in BOTH cohorts.

### KEGG: 45 pathways significant in both · 15 concordant · 30 discordant

#### Concordant activated in NAFLD (both cohorts NES > 0)

| Pathway | NES disc | NES val | padj disc | padj val |
|---|---|---|---|---|
| ECM-receptor interaction | +1.70 | +2.06 | 1.5e-04 | 4.4e-04 |
| Integrin signaling | +1.81 | +1.81 | 3.6e-08 | 1.6e-03 |
| Focal adhesion | +1.69 | +1.49 | 7.1e-07 | 3.9e-02 |
| Cornified envelope formation | +1.74 | +1.84 | 6.5e-06 | 3.7e-03 |
| Cholesterol metabolism | +1.50 | +1.73 | 2.6e-02 | 3.0e-02 |
| Fructose and mannose metabolism | +1.57 | +1.90 | 2.6e-02 | 8.6e-03 |
| Cytoskeleton in muscle cells | +1.68 | +1.51 | 1.2e-06 | 3.1e-02 |
| Human papillomavirus infection | +1.52 | +1.41 | 9.4e-06 | 4.7e-02 |

**No concordant suppressed KEGG pathways** were found in both cohorts.

#### Major discordant KEGG pathways (opposite direction)

| Pathway | NES disc | NES val | Likely explanation |
|---|---|---|---|
| IL-17 signaling pathway | +2.20 | −2.23 | Inflammatory signal inverted |
| TNF signaling pathway | +1.92 | −1.81 | Inflammatory signal inverted |
| NF-kappa B signaling pathway | +1.87 | −1.47 | Inflammatory signal inverted |
| Steroid hormone biosynthesis | −2.96 | +1.62 | Opposing metabolic shifts |
| Osteoclast differentiation | +1.89 | −1.69 | Immune infiltrate composition |
| Th17 cell differentiation | +1.91 | −1.37 | T-cell infiltrate difference |
| Non-alcoholic fatty liver disease | +1.43 | −1.63 | Pathway contains mixed signals |

### Hallmark: 10 gene sets significant in both · 2 concordant · 8 discordant

#### Concordant activated in NAFLD (both cohorts)

| Gene set | NES disc | NES val | padj disc | padj val |
|---|---|---|---|---|
| MYOGENESIS | +1.85 | +1.68 | 5.9e-10 | 3.3e-03 |
| APICAL_JUNCTION | +1.64 | +1.94 | 8.7e-06 | 9.2e-06 |

**No concordant suppressed Hallmark gene sets** were found in both cohorts.

#### Major discordant Hallmark gene sets

| Gene set | NES disc | NES val |
|---|---|---|
| TNFA_SIGNALING_VIA_NFKB | +2.14 | −3.04 |
| EPITHELIAL_MESENCHYMAL_TRANSITION | +1.98 | −1.41 |
| INFLAMMATORY_RESPONSE | +1.95 | −1.48 |
| HYPOXIA | +1.84 | −2.13 |
| KRAS_SIGNALING_UP | +1.69 | −2.00 |
| TGF_BETA_SIGNALING | +1.55 | −1.76 |

### Interpretation

**Robust cross-cohort signal (concordant in both):** ECM remodelling pathways
(ECM-receptor interaction, Integrin signaling, Focal adhesion), Cholesterol metabolism,
Myogenesis, and Apical junction are consistently activated in NAFLD across both datasets.
These reflect fibrosis-associated extracellular matrix remodelling and hepatic lipid
handling — the biological core of NAFLD.

**Large-scale discordance of inflammatory pathways:** IL-17, TNF, NF-kB, TNFA/NF-kB
(Hallmark), and related immune gene sets are activated in the discovery cohort but
suppressed in the validation cohort. This mirrors the gene-level discordance already
observed for NR4A1 and RGS1 (see Cross-Cohort Validation Summary above). Probable
causes:

1. **Cell-composition confound.** GSE135251 is heavily enriched for advanced-fibrosis
   samples (NASH F3/F4). In advanced NASH, hepatocytes are replaced by fibrotic ECM;
   the bulk RNA signal is dominated by hepatocyte/ECM transcriptomes, diluting the
   immune cell signal. The 10 control samples in GSE135251 are transcriptionally "pure"
   liver, so the ratio flips.

2. **Unbalanced control group.** Ten controls vs 206 NAFLD makes DESeq2 size-factor
   estimation sensitive to outliers in the control group, which can affect gene ranking
   and hence GSEA direction.

3. **These immune pathways require single-cell or deconvolution approaches** to be
   reliably interpreted across bulk RNA-seq datasets with different fibrosis-stage
   compositions.

**Conclusion:** The concordant ECM/adhesion/cholesterol pathways are the most
dataset-independent NAFLD signatures at the pathway level and represent the most
interpretable cross-cohort findings. Inflammatory pathway enrichment results are
cohort-composition-dependent and should be interpreted with caution in bulk RNA-seq.

---

## Output Files

| File | Description |
|---|---|
| `scripts/validation_gse135251.R` | Full pipeline: download → filter → DESeq2 → cross-cohort |
| `scripts/gsea_analysis.R` | GSEA: KEGG + Hallmark, apeglm LFC ranking, plots + results |
| `plots/validation_overlap.png` | Bar chart: discovery-only / both / validation-only gene counts |
| `plots/lfc_correlation.png` | LFC scatter (13,100 shared genes), purple = both sig, yellow = markers |
| `plots/gsea_kegg_dotplot.png` | KEGG GSEA dotplot (padj<0.05, split by direction) |
| `plots/gsea_kegg_ridgeplot.png` | KEGG GSEA ridgeplot (core-enrichment LFC distributions) |
| `plots/gsea_hallmark_dotplot.png` | Hallmark GSEA dotplot |
| `plots/gsea_hallmark_ridgeplot.png` | Hallmark GSEA ridgeplot |
| `results/gse135251_results.csv` | Full DESeq2 results (MLE + apeglm LFC, gene symbol) |
| `results/gsea_kegg_results.csv` | KEGG GSEA results (87 significant pathways) |
| `results/gsea_hallmark_results.csv` | Hallmark GSEA results (18 significant gene sets) |
| `results/validation_overlap_summary.csv` | Gene overlap stats, Fisher's test, concordance, correlation |
| `results/protein_coding_gene_list.rds` | Protein-coding Ensembl IDs (biomaRt, for pipeline reproducibility) |
| `data/GSE135251/` | Raw HTSeq count files — gitignored, re-downloadable via script |

---

## Tools and Citations

- **GEOquery** (Davis & Meltzer, 2007) — GEO data download
- **DESeq2** (Love, Huber & Anders, 2014) — differential expression
- **apeglm** (Zhu, Ibrahim & Love, 2018) — LFC shrinkage
- **clusterProfiler** (Wu et al., 2021) — GSEA (KEGG + custom gene sets)
- **msigdbr** — MSigDB Hallmark gene sets (H collection)
- **enrichplot** — dotplot and ridgeplot visualisation
- **org.Hs.eg.db** — Ensembl → gene symbol / Entrez ID mapping
- **ggplot2** + **ggrepel** — visualisation
