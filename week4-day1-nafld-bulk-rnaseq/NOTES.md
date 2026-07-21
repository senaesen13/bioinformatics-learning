# Week 4 Day 1 — NAFLD Bulk RNA-seq: Discovery Cohort

## Dataset: GSE162694

| Field | Value |
|---|---|
| GEO Accession | GSE162694 |
| Publication | Suppli et al. (2021), *Hepatology* |
| Platform | Illumina NovaSeq 6000 (GPL21290) |
| Tissue | Liver biopsy (human) |
| Total samples | 143 |
| Normal (healthy histology) | 31 |
| NAFLD (all fibrosis stages F0–F4) | 112 |
| Gene identifiers | Ensembl IDs (no version suffixes) |

**Sample breakdown by fibrosis stage:**
F0 = 35 · F1 = 30 · F2 = 27 · F3 = 8 · F4 = 12 · Normal = 31

> Validation cohort analysis (GSE135251) is in `week4-day2-nafld-validation-rnaseq/` — see that folder for methodology and results.
>
> Cross-cohort pathway comparison (GSEA concordance) is in `week4-day2-nafld-validation-rnaseq/NOTES.md` under **"Cross-Cohort Pathway Comparison (GSEA)"**. Short summary: ECM remodelling (ECM-receptor interaction, Integrin signaling, Focal adhesion) and Cholesterol metabolism are concordantly activated in NAFLD across both cohorts; inflammatory pathways (TNF, IL-17, NF-kB) are discordant — likely a cell-composition confound driven by fibrosis-stage differences between the two datasets.

---

## Methodology

> **Threshold update (2026-07-14):** Significance thresholds were changed from
> `padj < 0.01 & |MLE LFC| > 2` to `padj < 0.05 & |MLE LFC| > 1` to explore a
> broader gene set. All significant gene counts, tables, and statistics below reflect
> the new thresholds. GSEA was not re-run (it ranks all genes by LFC, not just the
> significant subset).

> **GSEA ranking update (2026-07-20):** The GSEA ranking metric was changed from
> `lfc_apeglm` to `sign(lfc_apeglm) * -log10(pvalue_mle)`. This signed significance
> score weights genes by both effect size direction and statistical confidence, and is
> a more robust standard ranking approach. GSEA was re-run for KEGG and Hallmark.
> KEGG: 164 pathways (was 167). Hallmark: 32 gene sets (was 29). See updated GSEA
> results below.

### Step 1 — Gene Biotype Annotation (applied before count filtering)

Protein-coding genes were identified using **biomaRt** `useEnsembl()` (Ensembl current via asia mirror, release ~111). For each Ensembl ID in the count matrix, `gene_biotype` was retrieved via `getBM()`. Only genes annotated as `protein_coding` were retained.

This step was applied **before** any count-level filtering — removing lncRNAs, pseudogenes, miRNAs, snoRNAs, snRNAs, etc. up front.

**Fallback order if biomaRt unavailable:** EnsDb.Hsapiens.v86 → AnnotationHub EnsDb.

### Step 2 — Mean-Count Filter

After protein-coding filtering, genes with `rowMeans(counts) < 10` were removed. This replaces the old `>10 counts in ≥3 samples` filter and better captures consistently expressed genes in heterogeneous sample sets.

### Step 3 — DESeq2 Differential Expression

**Design:** `~ condition` (reference level = Normal)  
**LFC shrinkage:** apeglm (Zhu, Ibrahim & Love, 2018)  
**Results table:** contains both MLE (un-shrunk) LFC and apeglm-shrunk LFC  

### Step 4 — Significance Threshold

**Filtering criterion for significant gene list:** `padj < 0.05` AND `|log2FoldChange_MLE| > 1`

The MLE (un-shrunk) LFC is used for the threshold, because apeglm shrinkage aggressively contracts fold-change estimates toward zero, causing biologically large effects to fall below a hard LFC cutoff. The apeglm LFC is reported in the output table (it provides better point estimates), but is not used for inclusion/exclusion filtering.

---

## Filtering Counts

| Stage | Genes |
|---|---|
| Starting count matrix | 31,683 |
| After protein-coding filter | 16,390 |
| After mean-count filter (rowMeans ≥ 10) | **15,637** |

**Biotypes removed (top categories):** lncRNA 7,155 · processed_pseudogene 2,587 · transcribed_unprocessed_pseudogene 642 · misc_RNA 633 · transcribed_processed_pseudogene 407 · snoRNA 218 · snRNA 209 · miRNA 94 · plus 2,748 Ensembl IDs unrecognised by biomaRt (also dropped).

---

## Significant Genes — padj < 0.05 AND |MLE LFC| > 1

| Direction | Count |
|---|---|
| Up in NAFLD | 470 |
| Down in NAFLD | 15 |
| **Total** | **485** |

**Context:** At padj < 0.05 there are 5,631 significant genes, but the median |MLE LFC| among them is modest; the |LFC| > 1 threshold captures genes with at least a 2-fold change. The 485 significant genes reflect the biology of comparing healthy liver to NAFLD across *all* fibrosis stages combined — the F0/F1 samples (n=65 of 112 NAFLD) are transcriptionally close to normal and dilute the signal from advanced disease. The predominance of up-regulated genes (470 up vs 15 down) is expected: NAFLD involves activation of fibrosis, inflammation, and lipid-handling programmes more than suppression of constitutively expressed hepatocyte genes at these fold-change magnitudes.

---

## TREM2 / SPP1 / GPNMB

| Gene | MLE log2FC | padj | In sig list? |
|---|---|---|---|
| **TREM2** | +2.48 | 5.2e-13 | **YES** |
| **SPP1** | +1.87 | 1.1e-07 | **YES** |
| **GPNMB** | +1.19 | 2.3e-09 | **YES** |

All three canonical NAFLD macrophage markers are now significant under the relaxed thresholds. Their MLE log2FC values (1.19–2.48) are genuine, not driven by outliers — apeglm-shrunken estimates are consistent (+2.41, +1.77, +1.14 respectively). The fold changes reflect dilution by whole-liver bulk RNA-seq across mixed fibrosis stages; scRNA-seq in macrophage subpopulations would show substantially larger effects. See day2 for validation cohort concordance.

---

## GSEA Results

**Method:** Genes ranked by `sign(lfc_apeglm) * -log10(pvalue_mle)` (NAFLD / Normal), descending.
clusterProfiler `gseKEGG` (organism = "hsa") + `GSEA` with MSigDB Hallmark (H collection).
BH-adjusted p-value cutoff < 0.05, minGSSize = 15, maxGSSize = 500.

**KEGG:** 164 significant pathways

| Direction | Top pathways (padj) |
|---|---|
| Activated in NAFLD | Rheumatoid arthritis (1.1e-08), Integrin signaling (1.1e-08), Cytoskeleton in muscle cells (1.1e-08), Phagocytosis (1.8e-08), Antigen processing and presentation (1.3e-07) |
| Suppressed in NAFLD | Valine/leucine/isoleucine degradation (2.0e-06), Drug metabolism — CYP450 (4.5e-06), Chemical carcinogenesis — DNA adducts (5.5e-06), Biosynthesis of cofactors (1.0e-05), Steroid hormone biosynthesis (1.5e-05) |

Compared to the previous lfc_apeglm ranking: immune/infection pathways (Rheumatoid arthritis, Phagocytosis, Antigen presentation) rise to the top activated set alongside ECM/Integrin signals. The suppressed set is now dominated by amino acid catabolism and xenobiotic metabolism — consistent with impaired hepatocyte metabolic function in NAFLD.

**Hallmark:** 32 significant gene sets

| Direction | Gene set | NES | padj |
|---|---|---|---|
| Activated | ALLOGRAFT_REJECTION | +2.35 | 8.3e-10 |
| Activated | TNFA_SIGNALING_VIA_NFKB | +2.34 | 8.3e-10 |
| Activated | EPITHELIAL_MESENCHYMAL_TRANSITION | +2.28 | 8.3e-10 |
| Activated | APOPTOSIS | +2.10 | 8.3e-10 |
| Activated | MYOGENESIS | +2.04 | 8.3e-10 |
| Activated | INFLAMMATORY_RESPONSE | +1.98 | 8.3e-10 |
| Suppressed | BILE_ACID_METABOLISM | −1.87 | 5.8e-05 |
| Suppressed | PROTEIN_SECRETION | −1.70 | 4.4e-04 |
| Suppressed | XENOBIOTIC_METABOLISM | −1.39 | 3.2e-03 |

The direction pattern is broadly consistent with the previous ranking. ALLOGRAFT_REJECTION and APOPTOSIS rise to the top activated set (reflecting immune infiltration and hepatocyte stress). BILE_ACID_METABOLISM remains the top suppressed set, confirming impaired hepatocyte bile acid synthesis as a robust NAFLD signal.

---

## Output Files

| File | Description |
|---|---|
| `scripts/deseq2_analysis.R` | Full discovery pipeline (biotype filter → DESeq2 → plots) |
| `scripts/gsea_analysis.R` | GSEA: KEGG + Hallmark, sign(lfc_apeglm)×-log10(pvalue_mle) ranking, plots + results |
| `plots/pca.png` | PCA of VST-normalised counts (protein-coding genes) |
| `plots/volcano.png` | Volcano plot: MLE LFC x-axis, MLE padj colour, thresholds padj<0.05 & \|LFC\|>1, TREM2/SPP1/GPNMB labelled |
| `plots/gsea_kegg_dotplot.png` | KEGG GSEA dotplot (164 pathways, padj<0.05, split by direction) |
| `plots/gsea_kegg_ridgeplot.png` | KEGG GSEA ridgeplot (core-enrichment LFC distributions) |
| `plots/gsea_hallmark_dotplot.png` | Hallmark GSEA dotplot (32 gene sets) |
| `plots/gsea_hallmark_ridgeplot.png` | Hallmark GSEA ridgeplot |
| `results/deseq2_results.csv` | Full DESeq2 results (MLE + apeglm LFC, gene symbol) |
| `results/significant_genes.csv` | 485 significant genes (padj<0.05, \|MLE LFC\|>1) |
| `results/gsea_kegg_results.csv` | KEGG GSEA results (164 significant pathways) |
| `results/gsea_hallmark_results.csv` | Hallmark GSEA results (32 significant gene sets) |
| `results/protein_coding_gene_list.rds` | 16,390 protein-coding Ensembl IDs (shared with day2 pipeline) |

---

## Tools and Citations

- **GEOquery** (Davis & Meltzer, 2007) — GEO data download
- **biomaRt** (Durinck et al., 2009) — gene biotype annotation
- **DESeq2** (Love, Huber & Anders, 2014) — differential expression
- **apeglm** (Zhu, Ibrahim & Love, 2018) — LFC shrinkage
- **clusterProfiler** (Wu et al., 2021) — GSEA (KEGG + custom gene sets)
- **msigdbr** — MSigDB Hallmark gene sets (H collection)
- **enrichplot** — dotplot and ridgeplot visualisation
- **org.Hs.eg.db** — Ensembl → gene symbol / Entrez ID mapping
- **ggplot2** + **ggrepel** — visualisation
