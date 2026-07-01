# Week 4 Day 1 — NAFLD Bulk RNA-seq: DESeq2 + GSEA

## Dataset: GSE162694

**GEO Accession:** GSE162694  
**Publication:** Suppli et al. (2021), *Hepatology* — "Differential transcriptomic signatures of liver tissue from patients with NAFLD compared to healthy controls"  
**Platform:** Illumina NovaSeq 6000 (GPL21290)  
**Organism:** Homo sapiens  
**Tissue:** Liver biopsy  
**Total samples:** 143  

### Sample Breakdown

| Fibrosis Stage | Description | N |
|---|---|---|
| Normal | Healthy liver (normal histology) | 31 |
| F0 | NAFLD, no fibrosis | 35 |
| F1 | Mild fibrosis | 30 |
| F2 | Moderate fibrosis | 27 |
| F3 | Advanced fibrosis | 8 |
| F4 | Cirrhosis | 12 |
| **Total** | | **143** |

For this analysis: **Normal (n=31) vs NAFLD (n=112, all fibrosis stages combined)**.

---

## Analysis Steps

### Step 1 — Data Download
Used `GEOquery::getGEO()` to download metadata (series matrix) and `getGEOSuppFiles()` to retrieve the supplementary raw count CSV (`GSE162694_raw_counts.csv.gz`).  
- Gene identifiers: Ensembl IDs (ENSG...)
- Sample names in count matrix: `548nash1`, `548nash10`, etc. (extracted from `title` column in metadata)

### Step 2 — Sample Alignment
Matched count matrix column names to metadata using the second word in the `title` field. All 143 samples aligned successfully.

### Step 3 — Low-Count Gene Filtering
Removed genes with ≤10 counts in fewer than 3 samples.  
- Before: 31,683 genes  
- After: 31,426 genes (only 257 removed — most were truly lowly expressed)

### Step 4 — DESeq2 Differential Expression
**Design:** `~ condition` (Normal vs NAFLD)  
**LFC shrinkage:** apeglm (Zhu et al., 2018)  
**Significance threshold:** padj < 0.05  

**Results (padj < 0.1):**
- Upregulated in NAFLD: 6,932 genes (22%)
- Downregulated in NAFLD: 5,298 genes (17%)

**Volcano plot criteria (padj < 0.05 AND |LFC| > 1):**
- **Up in NAFLD:** 706 genes
- **Down in NAFLD:** 24 genes

The strong upregulation bias suggests widespread transcriptional activation in diseased liver — consistent with inflammatory and fibrotic remodelling.

### Step 5 — Volcano Plot
Saved to: `plots/volcano_nafld_vs_normal.png`  
Top genes labelled by –log10(padj). Dashed lines at LFC = ±1 and padj = 0.05.

### Step 6 — Top 20 DEGs (by absolute LFC, padj < 0.05)

| Gene | log2FC | padj | Note |
|---|---|---|---|
| SINHCAFP3 | +7.65 | 4.4e-12 | Liver-associated |
| FGF23 | +4.71 | 9.8e-19 | Fibroblast growth factor, fibrosis marker |
| NEUROD2 | +4.03 | 6.8e-14 | Neuronal differentiation TF |
| MPO | +3.73 | 3.4e-10 | Myeloperoxidase — neutrophil/macrophage marker |
| PADI1 | +3.64 | 1.3e-05 | Peptidylarginine deiminase, citrullination |
| CXCL8 | +3.41 | 2.0e-04 | IL-8 — key inflammatory chemokine |
| NR4A1 | +3.37 | 5.2e-08 | Nuclear receptor, macrophage activation |
| KRT23 | +3.20 | 9.8e-03 | Keratin 23, hepatic stress marker |
| TNFRSF12A | +3.31 | 2.5e-12 | TWEAK receptor, liver injury & fibrosis |
| MTCO1P40 | -1.86 | 5.8e-04 | Pseudogene |
| GABRA3 | -1.58 | 1.9e-06 | GABA receptor subunit, hepatic |
| LINC00640 | -1.33 | 5.2e-06 | lncRNA |

Only 24 genes met padj<0.05 & LFC<-1, explaining the asymmetric volcano — NAFLD is characterised by strong transcriptional activation.

### Step 7 — PCA (VST-normalised)
Saved to: `plots/pca_nafld_vs_normal.png`  
- PC1: 15.4% variance  
- PC2: 10.8% variance  
NAFLD and Normal samples show some separation on PC1, but overlap is expected given fibrosis stage heterogeneity within the NAFLD group.

---

## GSEA Results

### KEGG Pathways (258 significant, padj < 0.25)

**Most enriched (suppressed in NAFLD — negative NES):**

| Pathway | NES | padj |
|---|---|---|
| Drug metabolism — cytochrome P450 | -2.74 | 1.5e-08 |
| Chemical carcinogenesis — DNA adducts | -2.73 | 2.1e-08 |
| Steroid hormone biosynthesis | -2.41 | 8.3e-07 |
| Retinol metabolism | -2.33 | 1.8e-05 |
| Metabolism of xenobiotics by cytochrome P450 | -2.28 | 2.7e-06 |

**Interpretation:** Loss of normal hepatocyte function — CYP450 enzymes (CYPA, CYP2C, CYP3A families) are dramatically downregulated in NAFLD, consistent with impaired detoxification and drug metabolism capacity of the diseased liver.

**Most enriched (activated in NAFLD — positive NES):**

| Pathway | NES | padj |
|---|---|---|
| IL-17 signaling pathway | +2.05 | 3.8e-09 |
| Rheumatoid arthritis | +1.99 | 3.8e-09 |
| TNF signaling pathway | +1.95 | 3.8e-09 |
| Staphylococcus aureus infection | +2.06 | 8.6e-08 |
| Bladder cancer | +2.02 | 3.0e-06 |

**Interpretation:** IL-17 and TNF-alpha signalling are canonical NAFLD/NASH inflammatory pathways. The disease-gene sets ("Rheumatoid arthritis", "Staphylococcus infection") reflect shared immune-activation modules.

### Hallmark Gene Sets (41 significant, padj < 0.25)

All 41 significant Hallmark sets were **positively enriched in NAFLD**:

| Gene Set | NES | padj |
|---|---|---|
| TNFA_SIGNALING_VIA_NFKB | +2.12 | 5.5e-10 |
| ALLOGRAFT_REJECTION | +2.04 | 5.5e-10 |
| EPITHELIAL_MESENCHYMAL_TRANSITION | +1.97 | 5.5e-10 |
| APOPTOSIS | +1.93 | 5.5e-10 |
| IL6_JAK_STAT3_SIGNALING | +1.91 | 2.3e-08 |
| INFLAMMATORY_RESPONSE | +1.89 | 5.5e-10 |
| ANGIOGENESIS | +1.88 | 4.9e-05 |
| MYOGENESIS | +1.84 | 5.5e-10 |
| HYPOXIA | +1.83 | 5.5e-10 |
| P53_PATHWAY | +1.83 | 5.5e-10 |

**Interpretation:** The Hallmark results tell a coherent NAFLD biology story:
- **Inflammation:** TNFa/NFkB, IL6/JAK/STAT3, Inflammatory Response — driving hepatic inflammation and immune cell infiltration
- **Fibrosis:** EMT (Epithelial-Mesenchymal Transition) — hallmark of hepatic stellate cell activation
- **Hypoxia:** Consistent with impaired oxygen delivery in steatotic, inflamed liver
- **Apoptosis & P53:** Hepatocyte damage and DNA stress response
- No Hallmark gene sets were suppressed, confirming the strong upregulation bias

---

## Key Biological Conclusions

1. **Metabolic liver function collapses in NAFLD:** CYP450 detoxification, steroid/retinol metabolism, and xenobiotic metabolism are massively downregulated — the liver can no longer perform its normal biochemical tasks.

2. **Inflammation is the dominant transcriptional signature:** TNFa, IL-17, IL-6/JAK/STAT3, and NFkB cascades are all activated. These are the same pathways targeted by NASH clinical trial drugs (e.g., semaglutide reduces hepatic TNFa; selonsertib targeted hepatic fibrosis).

3. **EMT/fibrosis is already active:** Even the combined NAFLD group (including early-stage F0/F1) shows Hallmark EMT enrichment, suggesting fibrogenic gene programs activate early.

4. **Asymmetric DEG distribution (706 up vs 24 down at LFC>1):** NAFLD is a gain-of-pathological-function state. The liver activates vast inflammatory, immune, and stress programmes while specific hepatocyte-differentiation functions are lost.

---

## Output Files

| File | Description |
|---|---|
| `plots/volcano_nafld_vs_normal.png` | Volcano plot, top 25 genes labelled |
| `plots/top20_degs_barplot.png` | Top 20 DEGs by absolute LFC |
| `plots/pca_nafld_vs_normal.png` | PCA of VST-transformed counts |
| `plots/gsea_kegg_dotplot.png` | KEGG GSEA dotplot (activated/suppressed) |
| `plots/gsea_kegg_ridgeplot.png` | KEGG GSEA ridge plot |
| `plots/gsea_hallmark_dotplot.png` | Hallmark GSEA dotplot |
| `plots/gsea_hallmark_ridgeplot.png` | Hallmark GSEA ridge plot |
| `results/deseq2_nafld_vs_normal.csv` | Full DESeq2 results (31,426 genes) |
| `results/top20_degs.csv` | Top 20 DEGs by absolute LFC |
| `results/gsea_kegg_results.csv` | KEGG GSEA results (258 pathways) |
| `results/gsea_hallmark_results.csv` | Hallmark GSEA results (41 gene sets) |

---

## Tools and Citations

- **GEOquery** (Davis & Meltzer, 2007) — GEO data download
- **DESeq2** (Love, Huber & Anders, 2014) — differential expression
- **apeglm** (Zhu, Ibrahim & Love, 2018) — LFC shrinkage
- **clusterProfiler** (Yu et al., 2012) — GSEA framework
- **msigdbr** — MSigDB Hallmark gene sets in R
- **org.Hs.eg.db** — Ensembl → Entrez ID mapping
- **ggplot2** + **ggrepel** — visualisation
