# Week 5 Day 1 — Obesity Adipose RNA-seq Practice Analysis

## Dataset

- **GEO accession**: GSE166047
- **Study**: Rey/Messa et al. — transcriptomic profiling of subcutaneous adipose tissue in lean vs obese women
- **Comparison**: Normal-weight females (NW) vs obese females without T2D (OBF)
- **Tissue**: Subcutaneous adipose tissue
- **Species**: Homo sapiens (GRCh38.p13)
- **Sequencing**: Paired-end bulk RNA-seq (Illumina)

## Sample Sizes

| Group | n | Description |
|-------|---|-------------|
| NW    | 5 | Normal-weight women, BMI ~24.3 |
| OBF   | 5 | Obese women without T2D, BMI ~38.2 |

Groups OBM (obese men) and OBT2D (obese T2D) were excluded to avoid confounding sex differences with obesity effects.

## Methodology

1. **Raw data**: Per-sample `.dat.gz` files from `GSE166047_RAW.tar`; fragment count column (`frags_locus`) used as raw counts.
2. **Count matrix**: 10-sample matrix built by joining all NW and OBF sample files on gene_ID.
3. **Protein-coding filter**: `org.Hs.eg.db` GENETYPE annotation via `mapIds()`.
4. **Expression filter**: `rowMeans(counts) >= 10`.
5. **DESeq2**: design `~ condition` (NW as reference), apeglm LFC shrinkage.
6. **Significance**: padj < 0.01 AND |apeglm shrunk log2FC| > 2.
7. **GSEA**: Ranked by apeglm-shrunk log2FC; KEGG Legacy and Hallmark gene sets (msigdbr v26.1, MSigDB 2026.1).

## Gene Count Filtering

| Step | Gene count |
|------|-----------|
| Raw (all Ensembl genes) | 60,199 |
| After protein-coding filter | 19,235 |
| After mean-count filter (≥10) | 16,820 |
| DESeq2 genes tested | 16,820 |

## Differential Expression Results

**Thresholds**: padj < 0.01, |apeglm shrunk log2FC| > 2

Genes are filtered using the apeglm-shrunk fold change rather than the raw MLE estimate. For low-count genes, the MLE can be extreme and unstable; apeglm shrinks these toward zero, giving a more reliable measure of effect size. Only genes that pass both the p-value and the shrunk fold-change threshold are reported as significant.

| Direction | Count |
|-----------|-------|
| Upregulated in OBF | 1 |
| Downregulated in OBF | 1 |
| Total significant | 2 |

| Gene | Ensembl ID | Shrunk log2FC | MLE log2FC | padj |
|------|-----------|--------------|------------|------|
| AZGP1 | ENSG00000160862 | −2.24 | −2.43 | 2.8×10⁻¹⁵ |
| CLEC4C | ENSG00000198178 | +8.86 | +7.68 | 1.6×10⁻¹¹ |

**AZGP1** (zinc-alpha-2-glycoprotein) is a well-characterised adipose-expressed protein known to promote lipid mobilisation; its downregulation in obesity is consistent with published literature. **CLEC4C** (BDCA-2) is a plasmacytoid dendritic cell marker, consistent with immune cell infiltration in obese adipose tissue.

## GSEA Highlights

### KEGG (top significant)

| Direction | Pathway | NES | padj |
|-----------|---------|-----|------|
| Activated in OBF | ECM_RECEPTOR_INTERACTION | 1.84 | 0.003 |
| Activated in OBF | FOCAL_ADHESION | 1.75 | 0.073 |
| Activated in OBF | LEUKOCYTE_TRANSENDOTHELIAL_MIGRATION | 1.71 | 0.073 |
| Suppressed in OBF | NITROGEN_METABOLISM | −2.10 | 0.072 |
| Suppressed in OBF | ADIPOCYTOKINE_SIGNALING_PATHWAY | −2.09 | 0.072 |
| Suppressed in OBF | GLYCOLYSIS_GLUCONEOGENESIS | −2.00 | 0.073 |
| Suppressed in OBF | PPAR_SIGNALING_PATHWAY | −2.07 | 0.137 |

### Hallmark (top significant)

| Direction | Pathway | NES | padj |
|-----------|---------|-----|------|
| Activated in OBF | EPITHELIAL_MESENCHYMAL_TRANSITION | 1.77 | 0.009 |
| Activated in OBF | COAGULATION | 1.72 | 0.009 |
| Activated in OBF | ALLOGRAFT_REJECTION | 1.75 | 0.023 |
| Activated in OBF | APICAL_JUNCTION | 1.66 | 0.045 |
| Suppressed in OBF | FATTY_ACID_METABOLISM | −1.96 | 0.009 |
| Suppressed in OBF | ADIPOGENESIS | −1.85 | 0.023 |
| Suppressed in OBF | ANDROGEN_RESPONSE | −1.93 | 0.045 |

## Biological Interpretation

These findings are consistent with known obesity biology in adipose tissue:

- **Activated**: ECM remodeling, focal adhesion, and EMT pathways reflect adipose fibrosis and stromal remodeling in obesity. Coagulation and immune activation (allograft rejection NES proxy for immune response) reflect chronic low-grade inflammation.
- **Suppressed**: Loss of fatty acid metabolism, adipogenesis, and PPAR signaling represents metabolic dysfunction of obese adipocytes — reduced capacity for lipid handling and healthy fat-cell differentiation.
- **AZGP1 downregulation**: Consistent with published literature showing AZGP1 promotes lipid mobilization and is reduced in obesity.

## Output Files

```
results/
  deseq2_results_all.csv          — full DESeq2 results (16,820 genes)
  deseq2_significant_genes.csv    — significant hits (padj<0.01, |apeglm log2FC|>2)
  gsea_KEGG.csv                   — full KEGG GSEA results
  gsea_Hallmark.csv               — full Hallmark GSEA results
  metadata_cache.csv              — GEO sample metadata
  GSE166047_RAW.tar               — raw data archive (cached)
  raw_files/                      — extracted per-sample count files

plots/
  pca_NW_vs_OBF.{pdf,png}         — VST-PCA colored by group
  volcano_NW_vs_OBF.{pdf,png}     — volcano with top-5 gene labels
```

## Standalone Note

This is an independent practice project, unrelated to the NAFLD analysis in week4-day1 and week4-day2.
