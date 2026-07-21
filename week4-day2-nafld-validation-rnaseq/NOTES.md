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

> **Threshold update (2026-07-14):** Significance thresholds were changed from
> `padj < 0.01 & |MLE LFC| > 2` to `padj < 0.05 & |MLE LFC| > 1` to explore a
> broader gene set (matching the update in the discovery cohort). All significant gene
> counts, tables, and statistics below reflect the new thresholds. GSEA was not re-run
> (it ranks all genes by LFC, not just the significant subset).

> **GSEA ranking update (2026-07-20):** The GSEA ranking metric was changed from
> `lfc_apeglm` to `sign(lfc_apeglm) * -log10(pvalue_mle)`. GSEA was re-run for KEGG
> and Hallmark. KEGG: 66 pathways (was 87). Hallmark: 16 gene sets (was 18). The GSEA
> Results section below reflects the new ranking. Note: the Cross-Cohort Pathway
> Comparison section was computed with the previous ranking and has not been re-run.

Identical pipeline to the discovery cohort:

1. **Protein-coding filter** — applied before count filtering. Same gene list derived
   from biomaRt in the discovery cohort (`protein_coding_gene_list.rds`).
   Keeps only `gene_biotype == protein_coding` Ensembl IDs.

2. **Mean-count filter** — `rowMeans(counts) >= 10`.

3. **DESeq2** — design `~ condition` (reference = Normal). Both MLE (un-shrunk) and
   apeglm-shrunk LFC are computed and stored in the results table.

4. **Significance threshold** — `padj < 0.05` AND `|log2FC_MLE| > 1`.
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

**Threshold:** padj < 0.05 AND |MLE log2FC| > 1

| Direction | Count |
|---|---|
| Up in NAFLD | 645 |
| Down in NAFLD | 413 |
| **Total** | **1058** |

The validation cohort yields more significant genes (1058 vs 485 in discovery) and
more balanced up/down representation. The larger down-regulated set in GSE135251
reflects its heavy enrichment for advanced fibrosis (NASH F3/F4 = 68/206 samples):
in advanced disease, hepatocyte-specific transcription programmes are suppressed as
parenchyma is replaced by fibrotic tissue, producing strong down-regulation signals
that are averaged out in the discovery cohort's more even fibrosis-stage distribution.

---

## TREM2 / SPP1 / GPNMB

| Gene | MLE log2FC | padj | In sig list? |
|---|---|---|---|
| **TREM2** | +2.52 | 2.2e-07 | **YES** |
| **SPP1** | +1.41 | 2.8e-03 | **YES** |
| GPNMB | +0.91 | 6.0e-03 | NO — LFC just below 1× |

**TREM2 and SPP1** are now significant in the validation cohort. Both replicate the
discovery cohort direction (upregulated in NAFLD): TREM2 at similar magnitude (+2.48
discovery, +2.52 validation); SPP1 with modest attenuation (+1.87 discovery, +1.41
validation).

**GPNMB** is statistically significant (padj = 6.0e-03) and upregulated, but its MLE
LFC (+0.91) falls just below the |LFC| > 1 cutoff. This is consistent with bulk
RNA-seq dilution effects — GPNMB is a macrophage/stellate cell marker that is
expressed in a subpopulation diluted by the predominant hepatocyte signal in
whole-liver biopsies.

---

## Cross-Cohort Validation Summary

**139 genes significant in both cohorts** (Ensembl ID match):

> **See also:** This Day 1 vs Day 2 comparison is consolidated alongside the other two
> pairwise comparisons (Day 3 vs Day 1 and Day 3 vs Day 2) in
> `week4-day4-nafld-cross-cohort-comparison/results/pairwise_comparison.md` (Comparison C),
> for a single-folder view of all three cohort combinations.

### Key concordantly up-regulated overlap genes (selected)

| Gene | disc log2FC | val log2FC | Notes |
|---|---|---|---|
| TREM2 | +2.48 | +2.52 | Kupffer/macrophage NAFLD marker |
| PRAMEF10 | +3.52 | +3.23 | Consistently largest FC in both |
| STMN2 | +2.68 | +2.57 | Stathmin-2; neural/hepatic stress |
| MMP9 | +2.28 | +1.40 | ECM remodelling |
| CHI3L1 | +2.23 | +2.19 | Fibrosis-associated lectin |
| PDZK1IP1 | +2.20 | +2.23 | MAP17; tubular/cholestatic signal |
| JUNB | +2.18 | −1.55 | Discordant (see below) |
| LIME1 | +1.99 | +2.34 | Immune adaptor |
| FOS | +1.92 | −6.28 | Discordant IEG (see below) |
| SPP1 | +1.87 | +1.41 | Macrophage/stellate marker |
| COL1A1 | +1.87 | +1.28 | Fibrosis collagen |
| FASN | +1.36 | +2.56 | Fatty acid synthesis |

### Key discordant overlap genes (up in discovery, down in validation)

| Gene | disc log2FC | val log2FC | Likely explanation |
|---|---|---|---|
| FOSB | +1.86 | −7.53 | Immediate-early gene; inverted in advanced fibrosis |
| FOS | +1.92 | −6.28 | Immediate-early gene; same reason |
| PPP1R3G | +1.10 | −3.80 | Glycogen metabolism; hepatocyte loss in F3/F4 |
| GADD45G | +1.26 | −3.16 | Stress/apoptosis response |
| DUSP1 | +1.30 | −3.23 | MAPK phosphatase; inflammatory regulation |
| NR4A1 | +2.04 | −3.72 | Nuclear receptor; immune/metabolic axis |
| RGS1 | +2.04 | −4.52 | G-protein regulator; immune cell composition |
| FOSB/FOS/JUN pattern | — | — | All immediate-early genes reversed in GSE135251 |

NR4A1, RGS1, FOS, FOSB, JUN, JUNB and several DUSPs are up in discovery but
strongly down in validation. The IEG pattern (FOS, FOSB, JUN, JUNB, DUSP1) reversed
in GSE135251 likely reflects hepatocyte transcriptional collapse in advanced fibrosis
(the validation cohort is NASH F3/F4-enriched) overriding the inflammatory IEG
induction signal.

### Summary statistics

| Metric | Value |
|---|---|
| Discovery sig genes | 485 |
| Validation sig genes | 1058 |
| Overlap | 139 genes |
| Overlap % of discovery | 28.7% |
| Fisher's exact OR | 5.11 |
| Fisher's exact p-value | 9.5e-43 |
| Direction concordance | 73.4% (102/139) |
| Spearman ρ (sig overlap, n=139) | 0.039 (p = 0.65) |
| Pearson r (sig overlap, n=139) | −0.048 (p = 0.58) |
| Pearson r (genome-wide, 13,100 genes) | 0.323 |

The overlap is highly enriched above chance (OR = 5.1, p = 9.5e-43) — with 139
genes the enrichment remains strong even though the OR is lower than at the |LFC|>2
threshold (which gave OR = 55.5 for just 7 genes). The near-zero Pearson and Spearman
within the 139 overlap genes reflects the discordant subset: all 485 discovery genes
are upregulated (lfc_mle > 1), but a subset (~37 of 139 in overlap) are strongly
downregulated in validation, and those genes (FOS, FOSB, RGS1) have very large
negative validation LFCs (−4 to −8) that dominate the correlation calculation. The
direction concordance of 73.4% is the more informative metric: it shows that the
majority of genes call the same direction in both cohorts. At the genome-wide level
(13,100 shared tested genes), Pearson r = 0.323 — unchanged from the strict threshold,
as expected (the full LFC distribution is not affected by the significance filter).

---

## GSEA Results (Validation Cohort)

**Method:** Genes ranked by `sign(lfc_apeglm) * -log10(pvalue_mle)` (NAFLD / Normal), descending.
clusterProfiler `gseKEGG` (organism = "hsa") + `GSEA` with MSigDB Hallmark (H collection).
BH-adjusted p-value cutoff < 0.05, minGSSize = 15, maxGSSize = 500.

**KEGG:** 66 significant pathways (was 87 with lfc_apeglm ranking)

| Direction | Top pathways (padj) |
|---|---|
| Activated in NAFLD | Lysosome biogenesis (6.6e-04), Carbon metabolism (7.6e-03), Motor proteins (7.6e-03), Biosynthesis of amino acids (1.7e-02), Integrin signaling (2.6e-02) |
| Suppressed in NAFLD | Ribosome (1.7e-08), Coronavirus disease (1.7e-08), IL-17 signaling (2.6e-02), Osteoclast differentiation (2.6e-02), FoxO signaling (3.8e-02) |

Note: Several immune/inflammatory pathways (IL-17, osteoclast differentiation) remain
**suppressed** in this cohort, opposite to the discovery cohort — see cross-cohort comparison below.

**Hallmark:** 16 significant gene sets (was 18 with lfc_apeglm ranking)

| Direction | Gene set | NES | padj |
|---|---|---|---|
| Activated | APICAL_JUNCTION | +1.66 | 2.7e-03 |
| Activated | MTORC1_SIGNALING | +1.60 | 3.5e-03 |
| Activated | BILE_ACID_METABOLISM | +1.58 | 1.2e-02 |
| Activated | MITOTIC_SPINDLE | +1.45 | 4.4e-02 |
| Suppressed | TNFA_SIGNALING_VIA_NFKB | −2.56 | 5.0e-09 |
| Suppressed | KRAS_SIGNALING_UP | −1.89 | 1.5e-05 |
| Suppressed | HYPOXIA | −1.61 | 3.8e-03 |
| Suppressed | TGF_BETA_SIGNALING | −1.61 | 4.4e-02 |
| Suppressed | UV_RESPONSE_UP | −1.58 | 1.2e-02 |

Notable change from previous ranking: **BILE_ACID_METABOLISM** switched from suppressed
to activated (+1.58). CHOLESTEROL_HOMEOSTASIS (previously #1 activated at NES +2.12)
drops below the padj<0.05 cutoff. TNFA_SIGNALING_VIA_NFKB remains the most strongly
suppressed set in this cohort (NES −2.56), consistent with the previous ranking.

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
| `scripts/gsea_analysis.R` | GSEA: KEGG + Hallmark, sign(lfc_apeglm)×-log10(pvalue_mle) ranking, plots + results |
| `plots/validation_overlap.png` | Bar chart: discovery-only / both / validation-only gene counts |
| `plots/lfc_correlation.png` | LFC scatter (13,100 shared genes), purple = both sig, yellow = markers |
| `plots/gsea_kegg_dotplot.png` | KEGG GSEA dotplot (padj<0.05, split by direction) |
| `plots/gsea_kegg_ridgeplot.png` | KEGG GSEA ridgeplot (core-enrichment LFC distributions) |
| `plots/gsea_hallmark_dotplot.png` | Hallmark GSEA dotplot |
| `plots/gsea_hallmark_ridgeplot.png` | Hallmark GSEA ridgeplot |
| `results/gse135251_results.csv` | Full DESeq2 results (MLE + apeglm LFC, gene symbol) |
| `results/gsea_kegg_results.csv` | KEGG GSEA results (66 significant pathways; re-run with signed significance ranking) |
| `results/gsea_hallmark_results.csv` | Hallmark GSEA results (16 significant gene sets; re-run with signed significance ranking) |
| `results/validation_overlap_summary.csv` | Gene overlap stats, Fisher's test, concordance, correlation (padj<0.05, \|LFC\|>1) |
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
