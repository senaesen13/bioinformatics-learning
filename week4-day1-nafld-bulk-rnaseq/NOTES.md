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

**Filtering criterion for significant gene list:** `padj < 0.01` AND `|log2FoldChange_MLE| > 2`

The MLE (un-shrunk) LFC is used for the threshold, because apeglm shrinkage aggressively contracts fold-change estimates toward zero, causing biologically large effects to fall below a hard |LFC| > 2 cutoff. The apeglm LFC is reported in the output table and used for the volcano plot x-axis (it provides better point estimates), but is not used for inclusion/exclusion filtering.

---

## Filtering Counts

| Stage | Genes |
|---|---|
| Starting count matrix | 31,683 |
| After protein-coding filter | 16,390 |
| After mean-count filter (rowMeans ≥ 10) | **15,637** |

**Biotypes removed (top categories):** lncRNA 7,155 · processed_pseudogene 2,587 · transcribed_unprocessed_pseudogene 642 · misc_RNA 633 · transcribed_processed_pseudogene 407 · snoRNA 218 · snRNA 209 · miRNA 94 · plus 2,748 Ensembl IDs unrecognised by biomaRt (also dropped).

---

## Significant Genes — padj < 0.01 AND |MLE LFC| > 2

| Direction | Count |
|---|---|
| Up in NAFLD | 24 |
| Down in NAFLD | 0 |
| **Total** | **24** |

**vs. expected ~1500:** The discrepancy reflects the nature of this dataset. At padj < 0.01 there are 5,173 significant genes, but the median |LFC| among them is only 0.43; the 99th percentile is 1.72. The |LFC| > 2 threshold captures only the most extreme 0.5% of effect sizes. This is consistent with comparing healthy liver to NAFLD across *all* fibrosis stages combined — the F0/F1 samples (n=65 of 112 NAFLD) are transcriptionally close to normal and dilute the signal from advanced disease. The expected ~1500 figure was not achievable without either lowering the LFC threshold (|LFC| > 1 gives 470 genes) or restricting to severe fibrosis only.

---

## TREM2 / SPP1 / GPNMB

| Gene | MLE log2FC | padj | In sig list? |
|---|---|---|---|
| **TREM2** | +2.48 | 5.2e-13 | **YES** |
| SPP1 | +1.87 | 1.1e-07 | NO — LFC below 2× |
| GPNMB | +1.19 | 2.3e-09 | NO — LFC below 2× |

**TREM2** is significant at the strict threshold.

**SPP1 and GPNMB** are statistically significant (padj < 0.01) and consistently upregulated, but their fold changes do not reach 2×. This is biologically expected: both genes are expressed in macrophage and stellate cell subpopulations that are diluted within whole-liver biopsies, especially at early fibrosis stages. See day2 for validation cohort results.

---

## Output Files

| File | Description |
|---|---|
| `scripts/deseq2_analysis.R` | Full discovery pipeline (biotype filter → DESeq2 → plots) |
| `scripts/gsea_analysis.R` | GSEA: KEGG + Hallmark, apeglm LFC ranking, plots + results |
| `plots/pca.png` | PCA of VST-normalised counts (protein-coding genes) |
| `plots/volcano.png` | Volcano plot: apeglm LFC x-axis, MLE padj colour, TREM2/SPP1/GPNMB labelled |
| `plots/gsea_kegg_dotplot.png` | KEGG GSEA dotplot (167 pathways, padj<0.05, split by direction) |
| `plots/gsea_kegg_ridgeplot.png` | KEGG GSEA ridgeplot (core-enrichment LFC distributions) |
| `plots/gsea_hallmark_dotplot.png` | Hallmark GSEA dotplot (29 gene sets) |
| `plots/gsea_hallmark_ridgeplot.png` | Hallmark GSEA ridgeplot |
| `results/deseq2_results.csv` | Full DESeq2 results (MLE + apeglm LFC, gene symbol) |
| `results/significant_genes.csv` | 24 significant genes (padj<0.01, \|MLE LFC\|>2) |
| `results/gsea_kegg_results.csv` | KEGG GSEA results (167 significant pathways) |
| `results/gsea_hallmark_results.csv` | Hallmark GSEA results (29 significant gene sets) |
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
