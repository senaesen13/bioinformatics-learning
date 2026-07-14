# Week 4 Day 3 — NAFLD Bulk RNA-seq: GSE130970 (Third Cohort)

## Dataset: GSE130970

| Field | Value |
|---|---|
| GEO Accession | GSE130970 |
| Platform | Illumina HiSeq 4000 (GPL16791) |
| Tissue | Human liver biopsy |
| Total samples | 78 |
| Normal (NAS = 0) | 4 |
| NAFLD (NAS 1–6) | 74 |
| Gene identifiers | Entrez gene IDs (salmon/tximport pipeline) |
| Quantification | Salmon, processed with tximport |

**Sample breakdown by NAFLD activity score (NAS):**
NAS 0 = 4 · NAS 1 = 5 · NAS 2 = 9 · NAS 3 = 18 · NAS 4 = 16 · NAS 5 = 18 · NAS 6 = 8

**Sample breakdown by fibrosis stage:**
F0 = 25 · F1 = 28 · F2 = 9 · F3 = 14 · F4 = 2

> This is the third independent NAFLD cohort, used alongside GSE162694 (Day 1, discovery)
> and GSE135251 (Day 2, validation). Cross-cohort pairwise comparisons are in
> `week4-day4-nafld-cross-cohort-comparison/`.

---

## Grouping Strategy

GSE130970 does not supply an explicit "healthy vs NAFLD" label. Disease severity is
captured by the **NAFLD activity score (NAS)**, which is the sum of steatosis grade
(0–3), lobular inflammation (0–2), and cytological ballooning (0–2).

- **Normal (control):** NAS = 0 (no active NAFLD histology; n = 4)
- **NAFLD:** NAS > 0 (n = 74)

This mapping is the most biologically defensible for a binary comparison: NAS = 0
means no histological evidence of active NAFLD. An alternative grouping (steatosis ≥ 2
vs steatosis = 0) gives a similar 8/70 split, but NAS is more clinically standard.

**Important caveat:** With only 4 control samples, individual gene-level power is
severely limited. DESeq2 replaces outliers only when `minReplicatesForReplace = 7`
(the default), so small-n groups are handled conservatively. The 180 significant genes
reflect the genes with the clearest biological effect despite this limitation. GSEA
(which uses the full ranked list, not the significant subset) is the more reliable
output from this dataset.

---

## Methodology

Identical pipeline to Days 1 and 2:

1. **Gene Biotype Filter:** org.Hs.eg.db `GENETYPE == "protein-coding"` applied to
   Entrez gene IDs (the native ID type in this dataset — no Ensembl conversion needed).

2. **Mean-Count Filter:** `rowMeans(counts) >= 10`.

3. **DESeq2:** design `~ condition` (reference = Normal). Both MLE (un-shrunk) and
   apeglm-shrunk LFC computed and stored.

4. **Significance Threshold:** `padj < 0.05` AND `|log2FC_MLE| > 1`.

5. **Sample alignment:** The count matrix column names (`440349.1.X_1` etc.) are
   internal sample identifiers that match the GEO series matrix `title` column
   exactly (aligned by numeric prefix). No samples were dropped.

### Filtering Counts

| Stage | Genes |
|---|---|
| Starting count matrix (Salmon/tximport, Entrez IDs) | 19,585 |
| After protein-coding filter (org.Hs.eg.db GENETYPE) | 18,591 |
| After mean-count filter (rowMeans ≥ 10) | **15,482** |

---

## Significant Gene Results

**Threshold:** padj < 0.05 AND |MLE log2FC| > 1

| Direction | Count |
|---|---|
| Up in NAFLD | 133 |
| Down in NAFLD | 47 |
| **Total** | **180** |

The limited significant gene count (180 vs 485 in Day 1 and 1058 in Day 2) directly
reflects the statistical power constraint of n=4 controls. The genome-wide LFC
correlation with the other cohorts remains informative (see Day 4), confirming that
the biological signal is present even when individual gene-level significance is lost.

---

## TREM2 / SPP1 / GPNMB

| Gene | MLE log2FC | apeglm log2FC | padj | In sig list? |
|---|---|---|---|---|
| TREM2 | +1.51 | +0.04 | 1.4e-01 | NO — padj > 0.05 |
| SPP1 | +0.18 | +0.02 | 8.97e-01 | NO — LFC and padj both below threshold |
| GPNMB | +0.30 | +0.04 | 6.05e-01 | NO — padj > 0.05 |

None of the three canonical NAFLD macrophage markers reach significance in this dataset.

Two compounding reasons:
1. **Statistical power:** 4 control samples provide very limited reference-group
   variance estimation. Even with n=74 NAFLD, the Normal group is the bottleneck.
2. **apeglm heavy shrinkage:** The dramatic gap between MLE LFC (+1.51 for TREM2)
   and apeglm (+0.04) indicates that apeglm's prior — fitted to all genes in the
   dataset — is aggressively shrinking genes with uncertain effect size estimates.
   This is the prior working as intended when data is sparse.

TREM2 shows a positive MLE LFC (+1.51) consistent with the other cohorts, but the
low apeglm value and borderline padj (0.14) prevent it from meeting the threshold.
In the genome-wide LFC correlation with Day 1 and Day 2, TREM2 shows the expected
positive direction in all cohorts.

---

## GSEA Results

**Method:** Genes ranked by apeglm-shrunken log2FC (NAFLD / Normal), descending.
Gene IDs are Entrez (native to this dataset — no Ensembl→Entrez conversion needed).
clusterProfiler `gseKEGG` + `GSEA` with MSigDB Hallmark. BH-adjusted p < 0.05,
minGSSize = 15, maxGSSize = 500.

### KEGG: 94 significant pathways

| Direction | Top pathways (padj) |
|---|---|
| Activated in NAFLD | Cytoskeleton in muscle cells (1.7e-08), Motor proteins (1.7e-08), Cardiac muscle contraction (6.0e-08), Dilated cardiomyopathy (9.4e-06), Hypertrophic cardiomyopathy (4.1e-05) |
| Suppressed in NAFLD | Tryptophan metabolism (4.2e-04), One carbon pool by folate (2.0e-03), Chemical carcinogenesis–DNA adducts (2.2e-03), Linoleic acid metabolism (2.3e-03), Xenobiotic metabolism by CYP450 (8.2e-03) |

The top activated pathways (cardiac/muscle cytoskeleton) may appear surprising in
a liver dataset. This signal likely reflects the transcriptional overlap between
hepatic stellate cell activation (a core NAFLD fibrosis mechanism) and smooth
muscle/cytoskeletal gene programmes — stellate cells upregulate myosin heavy chains,
troponins, and actomyosin regulatory genes as they become myofibroblasts. This
is a known confound in liver GSEA when fibrosis-stage F1/F2 samples dominate.

The suppressed metabolic pathways (tryptophan catabolism, folate one-carbon metabolism,
fatty acid metabolism, xenobiotic CYP450) are genuine hepatic function signatures:
in NAFLD, hepatocyte-specific metabolic programmes are progressively impaired.

### Hallmark: 17 significant gene sets

| Direction | Gene set | NES | padj |
|---|---|---|---|
| Activated | MYOGENESIS | +2.55 | 5.0e-09 |
| Activated | APICAL_JUNCTION | +1.94 | 6.6e-04 |
| Activated | ALLOGRAFT_REJECTION | +1.85 | 4.0e-03 |
| Activated | TNFA_SIGNALING_VIA_NFKB | +1.73 | 2.0e-02 |
| Activated | INTERFERON_GAMMA_RESPONSE | +1.71 | 2.0e-02 |
| Activated | IL6_JAK_STAT3_SIGNALING | +1.79 | 2.9e-02 |
| Activated | APOPTOSIS | +1.73 | 3.2e-02 |
| Activated | EPITHELIAL_MESENCHYMAL_TRANSITION | +1.63 | 3.2e-02 |
| Activated | P53_PATHWAY | +1.61 | 3.6e-02 |
| Activated | INTERFERON_ALPHA_RESPONSE | +1.71 | 3.8e-02 |
| Suppressed | BILE_ACID_METABOLISM | −1.66 | 2.0e-02 |

**MYOGENESIS** and **APICAL_JUNCTION** are concordant with both Day 1 and Day 2
(the only two Hallmark gene sets concordant across all three cohorts).
**TNFA_SIGNALING_VIA_NFKB** is activated here and in Day 1, but suppressed in Day 2
(the same discordance pattern seen between days 1 and 2).
**BILE_ACID_METABOLISM** suppressed is a characteristic hepatocyte
signature — bile acid synthesis genes are impaired in NAFLD and their suppression
is a reliable marker of hepatocyte dysfunction.

---

## Output Files

| File | Description |
|---|---|
| `scripts/deseq2_analysis.R` | Full pipeline: metadata → count loading → DESeq2 → plots |
| `scripts/gsea_analysis.R` | GSEA: KEGG + Hallmark (Entrez IDs, no conversion needed) |
| `plots/pca.png` | PCA of VST-normalised counts (protein-coding genes) |
| `plots/volcano.png` | Volcano plot: MLE LFC, thresholds padj<0.05 & \|LFC\|>1, TREM2/SPP1/GPNMB |
| `plots/gsea_kegg_dotplot.png` | KEGG GSEA dotplot (94 pathways) |
| `plots/gsea_kegg_ridgeplot.png` | KEGG GSEA ridgeplot |
| `plots/gsea_hallmark_dotplot.png` | Hallmark GSEA dotplot (17 gene sets) |
| `plots/gsea_hallmark_ridgeplot.png` | Hallmark GSEA ridgeplot |
| `results/gse130970_results.csv` | Full DESeq2 results (MLE + apeglm LFC, gene symbol) |
| `results/significant_genes.csv` | 180 significant genes (padj<0.05, \|MLE LFC\|>1) |
| `results/gsea_kegg_results.csv` | KEGG GSEA results (94 significant pathways) |
| `results/gsea_hallmark_results.csv` | Hallmark GSEA results (17 significant gene sets) |

---

## Tools and Citations

- **GEOquery** (Davis & Meltzer, 2007) — GEO data download
- **DESeq2** (Love, Huber & Anders, 2014) — differential expression
- **apeglm** (Zhu, Ibrahim & Love, 2018) — LFC shrinkage
- **org.Hs.eg.db** — Entrez gene biotype annotation and symbol mapping
- **clusterProfiler** (Wu et al., 2021) — GSEA (KEGG + custom gene sets)
- **msigdbr** — MSigDB Hallmark gene sets (H collection)
- **enrichplot** — dotplot and ridgeplot visualisation
- **ggplot2** + **ggrepel** — visualisation
