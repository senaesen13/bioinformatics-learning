# Week 4 Day 1 — Bulk RNA-seq / Microarray Analysis: NAFLD Liver (GSE89632)

---

## 1. Dataset

**Title:** Gene expression profiling of human liver in NAFLD

**GEO Accession:** [GSE89632](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE89632)

| Property | Value |
|---|---|
| Organism | *Homo sapiens* |
| Tissue | Human liver biopsy |
| Platform | GPL14951 — Illumina HumanHT-12 V4.0 expression BeadChip |
| Data type | Microarray (log2 fluorescence intensity) |
| Total samples | 63 |
| Healthy controls (HC) | 24 |
| Simple steatosis (SS) | 20 |
| NASH | 19 |

Samples are liver biopsies from patients undergoing bariatric surgery or elective abdominal surgery. Histological grading was used to assign NAFLD stage (HC → SS → NASH). All samples are adult human, tissue from a single-channel Illumina BeadChip array.

---

## 2. Original Paper Citation

> Ahrens M, Ammerpohl O, von Schönfels W, Kolarova J, Bens S, Itzel T, Teufel A, Herrmann A, Brosch M, Hinrichsen H, Erhart W, Egberts J, Sipos B, Schreiber S, Häsler R, Stickel F, Becker T, Krawczak M, Röcken C, Siebert R, Schafmayer C, Hampe J. **DNA methylation analysis in nonalcoholic fatty liver disease suggests distinct disease-specific and remodeling signatures after bariatric surgery.** *Cell Metabolism.* 2013 Oct 1;18(4):296–302. doi:10.1016/j.cmet.2013.07.004

The paper used both DNA methylation arrays and gene expression arrays (this dataset) to identify epigenetically regulated genes in NAFLD progression. The expression data (GSE89632) captures the transcriptional landscape across the NAFLD severity spectrum.

---

## 3. Analysis Tool: limma

**Package:** `limma` (Linear Models for Microarray Data)

> Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, Smyth GK. **limma powers differential expression analyses for RNA-sequencing and microarray studies.** *Nucleic Acids Research.* 2015;43(7):e47. doi:10.1093/nar/gkv007

### Why limma and not DESeq2

This is a common question — both produce log2 fold-changes and adjusted p-values, but they are designed for fundamentally different data types:

| | **limma** | **DESeq2** |
|---|---|---|
| Input data | Log2-normalised continuous intensities (microarray) | Raw integer counts (RNA-seq) |
| Statistical model | Linear model + empirical Bayes variance shrinkage | Negative binomial GLM |
| Assumption | Approximately Gaussian on log2 scale | Overdispersed discrete counts |
| Variance estimation | Borrows strength across genes via eBayes | Dispersion estimated per-gene, shrunk toward trend |

**GSE89632 is microarray data** (Illumina BeadChip fluorescence intensities). These are continuous log2 values, not read counts. Feeding them into DESeq2 would violate its core statistical assumptions (negative binomial distribution requires non-negative integers). limma was built specifically for this data type.

That said, `limma-voom` is now also routinely used on RNA-seq count data as an alternative to DESeq2, because voom transforms counts to log-CPM with precision weights, restoring the Gaussian assumption. For this dataset, we use plain `limma` on the already-log2-normalised matrix.

### Model specification

```r
design  <- model.matrix(~ 0 + condition, data = pData(eset))  # group-means
fit     <- lmFit(ex_gene, design)
contrasts <- makeContrasts(
  NASH_vs_HC = NASH - HC,
  SS_vs_HC   = SS   - HC,
  levels     = design
)
fit2    <- contrasts.fit(fit, contrasts)
fit2    <- eBayes(fit2)   # moderated t-statistics via empirical Bayes
```

**Significance thresholds:** FDR (Benjamini-Hochberg) < 0.05 AND |log2FC| > 0.5

---

## 4. Results Summary

### Differential Expression

| Comparison | Up-regulated | Down-regulated | Total DEGs |
|---|---|---|---|
| NASH vs HC | 1,462 | 1,514 | **2,976** |
| SS vs HC | 1,415 | 1,674 | **3,089** |

### Top DEGs — NASH vs Healthy Controls

| Gene | logFC | Direction | Biological role |
|---|---|---|---|
| FOSB | −4.69 | Down | AP-1 transcription factor, stress response |
| MIR21 | −2.11 | Down | microRNA-21; pro-fibrotic regulator |
| AXUD1 | −1.82 | Down | AXIN1-induced, Wnt target |
| JUNB | −2.86 | Down | AP-1 component, inflammation |
| SOCS3 | −3.21 | Down | JAK-STAT inhibitor; loss → cytokine hyperactivation |
| MYC | −2.21 | Down | Proliferation/oncogene |
| TYMS | +1.56 | Up | Thymidylate synthase; cell cycle |
| FMO1 | +2.34 | Up | Flavin monooxygenase; drug/lipid metabolism |

### Key GSEA Findings

GSEA was run with the **MSigDB Hallmark** collection (50 gene sets) and **KEGG pathways** using the limma t-statistic as the ranking metric. Both comparisons (NASH vs HC and SS vs HC) showed 39/50 Hallmark sets significant at FDR < 0.25, consistent with broad transcriptional rewiring in diseased liver.

**Top enriched pathways in NASH vs HC (positive NES = activated in NASH):**

| Pathway | NES | Biological interpretation |
|---|---|---|
| TNFα signalling via NFκB | positive | Core inflammatory activation |
| IL6/JAK/STAT3 signalling | positive | Cytokine-driven hepatocyte stress |
| IL2/STAT5 signalling | positive | Immune cell co-activation |
| Hypoxia | positive | Lipotoxicity and tissue oxygen stress |
| Inflammatory response | positive | Broad innate immune activation |

**Interpretation:** The NASH transcriptome is dominated by inflammatory rewiring — NFκB and JAK-STAT pathways are robustly activated, consistent with the transition from simple fat accumulation (SS) to hepatocellular injury (NASH). The hypoxia signature is notable: hepatic lipid overload impairs mitochondrial oxygen consumption, creating a hypoxic microenvironment that amplifies inflammation.

The downregulation of SOCS3 (a JAK-STAT brake) is mechanistically important — loss of this feedback inhibitor would sustain IL-6 and TNFα signalling even after the initial stimulus resolves.

---

## 5. Output Files

### Plots (`plots/`)

| File | Content |
|---|---|
| `00_sample_distributions.png` | Per-sample boxplot to verify normalisation |
| `01_PCA.png` | PCA coloured by HC / SS / NASH |
| `02_volcano_NASH_vs_HC.png` | Volcano plot, top DEGs labelled |
| `03_volcano_SS_vs_HC.png` | Volcano plot for SS vs HC |
| `04_heatmap_top40_NASH.png` | Row-scaled heatmap, top 40 NASH DEGs |
| `05_GSEA_hallmark_NASH_vs_HC.png` | Hallmark dot plot — NASH |
| `06_GSEA_KEGG_NASH_vs_HC.png` | KEGG dot plot — NASH |
| `07_GSEA_hallmark_SS_vs_HC.png` | Hallmark dot plot — SS |
| `08_GSEA_KEGG_SS_vs_HC.png` | KEGG dot plot — SS |
| `09_gseaplot_NASH_hallmark_*.png` | Running-sum plots, top 5 Hallmark sets (NASH) |
| `10_gseaplot_NASH_KEGG_*.png` | Running-sum plots, top 5 KEGG pathways (NASH) |

### Results (`results/`)

| File | Content |
|---|---|
| `GSE89632_eset.rds` | Full ExpressionSet object (R) |
| `GSE89632_sample_metadata.csv` | 63 × 78 sample annotation |
| `GSE89632_feature_metadata.csv` | 20,819 probe annotations |
| `DEG_NASH_vs_HC.csv` | Full ranked DE table, NASH vs HC |
| `DEG_SS_vs_HC.csv` | Full ranked DE table, SS vs HC |
| `top20_DEG_NASH_vs_HC.csv` | Top 20 significant DEGs, NASH |
| `top20_DEG_SS_vs_HC.csv` | Top 20 significant DEGs, SS |
| `GSEA_hallmark_NASH_vs_HC.csv` | MSigDB Hallmark GSEA, NASH |
| `GSEA_hallmark_SS_vs_HC.csv` | MSigDB Hallmark GSEA, SS |
| `GSEA_KEGG_NASH_vs_HC.csv` | KEGG GSEA, NASH |
| `GSEA_KEGG_SS_vs_HC.csv` | KEGG GSEA, SS |
| `top20_DEG_NASH_vs_HC.csv` | Top 20 NASH DEGs |
| `top20_DEG_SS_vs_HC.csv` | Top 20 SS DEGs |

---

## 6. References

1. **Dataset:** Ahrens M, et al. DNA methylation analysis in nonalcoholic fatty liver disease suggests distinct disease-specific and remodeling signatures after bariatric surgery. *Cell Metabolism.* 2013;18(4):296–302.

2. **limma:** Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, Smyth GK. limma powers differential expression analyses for RNA-sequencing and microarray studies. *Nucleic Acids Research.* 2015;43(7):e47.

3. **clusterProfiler:** Wu T, Hu E, Xu S, Chen M, Guo P, Dai Z, Feng T, Zhou L, Tang W, Zhan L, Fu X, Liu S, Bo X, Yu G. clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. *The Innovation.* 2021;2(3):100141.

4. **GSEA:** Subramanian A, Tamayo P, Mootha VK, Mukherjee S, Ebert BL, Gillette MA, Paulovich A, Pomeroy SL, Golub TR, Lander ES, Mesirov JP. Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. *PNAS.* 2005;102(43):15545–15550.

5. **MSigDB:** Liberzon A, et al. The Molecular Signatures Database Hallmark Gene Set Collection. *Cell Systems.* 2015;1(6):417–425.

6. **GEO:** Barrett T, et al. NCBI GEO: archive for functional genomics data sets—update. *Nucleic Acids Research.* 2013;41:D991–D995.
