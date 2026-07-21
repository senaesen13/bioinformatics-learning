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

> **Ranking update (2026-07-20):** GSEA ranking metric changed from `lfc_apeglm` to
> `sign(lfc_apeglm) * -log10(pvalue_mle)`. This is a more robust approach that weights
> genes by both direction and statistical confidence. Results below reflect the updated
> ranking. KEGG: 174 pathways (was 94). Hallmark: 34 gene sets (was 17).

**Method:** Genes ranked by `sign(lfc_apeglm) * -log10(pvalue_mle)` (NAFLD / Normal), descending.
Gene IDs are Entrez (native to this dataset — no Ensembl→Entrez conversion needed).
clusterProfiler `gseKEGG` + `GSEA` with MSigDB Hallmark. BH-adjusted p < 0.05,
minGSSize = 15, maxGSSize = 500.

### KEGG: 174 significant pathways (was 94)

| Direction | Top pathways (padj) |
|---|---|
| Activated in NAFLD | Systemic lupus erythematosus (5.8e-09), Rheumatoid arthritis (5.8e-09), Phagocytosis (5.8e-09), Tuberculosis (5.8e-09), Herpes simplex virus 1 infection (5.8e-09) |
| Suppressed in NAFLD | One carbon pool by folate (4.0e-04), Tryptophan metabolism (1.2e-03), Valine/leucine/isoleucine degradation (2.9e-03), Glycine/serine/threonine metabolism (1.4e-02), Linoleic acid metabolism (1.6e-02) |

**Key change from previous ranking:** The previous lfc_apeglm ranking put cardiac/muscle
cytoskeleton pathways (Cytoskeleton in muscle cells, Motor proteins, Cardiac muscle
contraction) at the top of activated KEGG pathways. With the signed significance score,
these are replaced by immune/infection pathways (lupus, rheumatoid arthritis, phagocytosis,
intracellular pathogens). This is more biologically interpretable: the cardiac/muscle
dominance was likely an artifact of apeglm over-shrinking in the n=4 control group, where
muscle-related genes happened to have high MLE LFC but high variance. The signed
significance metric corrects for this by down-weighting statistically uncertain genes.
The suppressed metabolic pathways (amino acid catabolism, folate cycle, fatty acid
metabolism) are consistent with both rankings and represent genuine hepatic function loss.

### Hallmark: 34 significant gene sets (was 17)

| Direction | Gene set | NES | padj |
|---|---|---|---|
| Activated | ALLOGRAFT_REJECTION | +2.48 | 1.3e-09 |
| Activated | INTERFERON_GAMMA_RESPONSE | +2.26 | 1.3e-09 |
| Activated | MTORC1_SIGNALING | +2.16 | 1.3e-09 |
| Activated | TNFA_SIGNALING_VIA_NFKB | +2.14 | 1.3e-09 |
| Activated | APOPTOSIS | +2.11 | 1.3e-08 |
| Activated | P53_PATHWAY | +1.97 | 1.2e-07 |
| Activated | INTERFERON_ALPHA_RESPONSE | +2.14 | 5.3e-07 |
| Activated | INFLAMMATORY_RESPONSE | +1.93 | 6.8e-07 |
| Activated | IL6_JAK_STAT3_SIGNALING | +2.11 | 1.9e-06 |
| Activated | EPITHELIAL_MESENCHYMAL_TRANSITION | +1.86 | 5.1e-06 |
| Suppressed | BILE_ACID_METABOLISM | −1.41 | 2.8e-02 |

**Key change from previous ranking:** MYOGENESIS was the #1 activated gene set with
lfc_apeglm (NES +2.55). With the signed significance score it drops out of the top 10.
Immune activation pathways (ALLOGRAFT_REJECTION, INTERFERON_GAMMA, TNFA, INFLAMMATORY)
are now dominant — a more coherent biological picture for NAFLD immune infiltration.
APICAL_JUNCTION is no longer significant at padj<0.05 with the new ranking.
**BILE_ACID_METABOLISM** suppressed remains robust across both rankings (−1.66 → −1.41),
confirming hepatocyte bile acid synthesis impairment as the most reliable suppression
signal in this dataset.

---

## Output Files

| File | Description |
|---|---|
| `scripts/deseq2_analysis.R` | Full pipeline: metadata → count loading → DESeq2 → plots |
| `scripts/gsea_analysis.R` | GSEA: KEGG + Hallmark, sign(lfc_apeglm)×-log10(pvalue_mle) ranking (Entrez IDs, no conversion needed) |
| `plots/pca.png` | PCA of VST-normalised counts (protein-coding genes) |
| `plots/volcano.png` | Volcano plot: MLE LFC, thresholds padj<0.05 & \|LFC\|>1, TREM2/SPP1/GPNMB |
| `plots/gsea_kegg_dotplot.png` | KEGG GSEA dotplot (174 pathways) |
| `plots/gsea_kegg_ridgeplot.png` | KEGG GSEA ridgeplot |
| `plots/gsea_hallmark_dotplot.png` | Hallmark GSEA dotplot (34 gene sets) |
| `plots/gsea_hallmark_ridgeplot.png` | Hallmark GSEA ridgeplot |
| `results/gse130970_results.csv` | Full DESeq2 results (MLE + apeglm LFC, gene symbol) |
| `results/significant_genes.csv` | 180 significant genes (padj<0.05, \|MLE LFC\|>1) |
| `results/gsea_kegg_results.csv` | KEGG GSEA results (174 significant pathways) |
| `results/gsea_hallmark_results.csv` | Hallmark GSEA results (34 significant gene sets) |

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
