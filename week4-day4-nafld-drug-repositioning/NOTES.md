# Week 4 Day 4 — NAFLD Drug Repositioning

## What is Drug Repositioning?

Drug repositioning (also called drug repurposing) finds new therapeutic uses for existing approved drugs. Instead of developing a drug from scratch (10–15 years, $1–2B), we ask: "which existing drug has a transcriptional effect that *opposes* the disease signature?" If a drug reverses the gene expression changes caused by a disease, it is a candidate to treat that disease.

**The Connectivity Map (CMap) approach:**  
1. Characterise the disease transcriptome: which genes are up/down in NAFLD vs healthy liver?  
2. For each drug, retrieve its gene expression signature from a reference database.  
3. Find drugs whose signature is *negatively correlated* with the disease signature — i.e., genes the drug **lowers** are the same genes that go **up** in disease.

---

## Data

**Disease signature:** DESeq2 results from Week 4 Day 1 (GSE130970, 15 NAFLD vs 15 Normal human liver biopsies).  
- 21,545 protein-coding genes ranked by `log2FC × −log10(padj)`  
- 597 significant DEGs (padj<0.05, |lFC|>1): **581 upregulated / 16 downregulated**  
- The NAFLD signature is heavily dominated by upregulated genes (inflammatory, fibrogenic, lipogenic pathways)

**Drug database:** MSigDB C2 Chemical and Genetic Perturbations (CGP), accessed via `msigdbr` v26.1:  
- 3,555 gene sets from curated perturbation experiments  
- Each gene set = genes up or down in response to a drug/chemical/genetic perturbation  
- Filtered to gene sets with 15–500 members: **2,749 tested**

---

## Methods

### Ranking Metric

```
rank_metric = log2FoldChange × −log10(padj_clamped)
```

Positive = upregulated in NAFLD; Negative = downregulated. Genes at the top drove the disease; genes at the bottom are what NAFLD suppresses.

### GSEA with fgsea

`fgsea` (Korotkevich et al.) runs fast preranked GSEA using an adaptive multi-level sampling algorithm. 10,000 permutations, FDR correction by BH.

### Drug Repositioning Logic (Type A — primary)

Because the NAFLD signature has 581 up vs 16 down genes, the most informative candidates come from **_DN gene sets with positive NES**:

```
Drug_X_DN (genes drug downregulates) ──GSEA──► positive NES
= these "downregulated-by-drug" genes are enriched AMONG NAFLD-upregulated genes
= Drug X suppresses what NAFLD overexpresses
= Drug X OPPOSES NAFLD ✓
```

This is equivalent to the CMap "negative correlation" between drug and disease.

---

## Results

### Scale of analysis
| Metric | Value |
|--------|-------|
| Gene sets tested | 2,749 |
| Significant (FDR<0.05) | 675 |
| Type-A drug candidates (_DN, FDR<0.05) | 176 |

### Top Drug/Perturbation Candidates

| Gene Set | NES | FDR | Interpretation |
|----------|-----|-----|----------------|
| DORN_ADENOVIRUS_INFECTION_*_DN | 1.51 | 0.006 | Adenovirus suppresses genes elevated in NAFLD (viral-induced immune/metabolic reprogramming) |
| HERNANDEZ_MITOTIC_ARREST_BY_DOCETAXEL_2_DN | 1.48 | 0.023 | **Docetaxel** (taxane): downregulates pro-fibrotic and cell-cycle genes shared with NAFLD signature |
| HAN_JNK_SIGNALING_DN | 1.46 | 0.023 | **JNK inhibition**: JNK (c-Jun N-terminal kinase) drives NAFLD→NASH progression; its inhibition reverses the signature |
| GESERICK_TERT_TARGETS_DN | 1.48 | 0.020 | **TERT targets**: telomerase-related proliferative genes up in NAFLD |
| URS_ADIPOCYTE_DIFFERENTIATION_DN | 1.45 | 0.029 | Genes suppressed during adipocyte differentiation overlap NAFLD signature |
| MA_RAT_AGING_DN | 1.44 | 6.6e-6 | Aging-related genes elevated in NAFLD (NAFLD accelerates hepatic ageing) |
| XU_HGF_SIGNALING_NOT_VIA_AKT1_48HR_DN | 1.44 | 0.034 | HGF (hepatocyte growth factor) signaling suppresses pro-inflammatory/NAFLD genes |

### Key Biological Interpretation

**JNK signaling (HAN_JNK_SIGNALING_DN):**  
JNK1/2 activation is a central driver of NAFLD progression. Free fatty acids and lipotoxic intermediates (ceramides, diacylglycerols) activate JNK in hepatocytes and Kupffer cells, driving insulin resistance, steatohepatitis, and stellate cell activation. JNK inhibitors (e.g. CC-930/tanzisertib) have been tested in NASH clinical trials, validating this hit.

**Docetaxel (HERNANDEZ_MITOTIC_ARREST_BY_DOCETAXEL_2_DN):**  
Docetaxel and other taxanes suppress cell-cycle entry genes (CDK1, CCNB1, TOP2A) that are aberrantly elevated in NAFLD-related hepatocyte stress. While taxanes are too toxic for NAFLD, CDK inhibitors derived from this biology are active areas of research.

**HGF / MET signaling (XU_HGF_SIGNALING_NOT_VIA_AKT1_48HR_DN):**  
HGF activates MET/PI3K→STAT3 pathways that suppress pro-inflammatory hepatic gene expression. HGF supplementation or MET activators have hepatoprotective effects in NAFLD animal models.

**Adenovirus infection genes:**  
Adenovirus 36 infection has been linked to obesity and fatty liver in human and animal studies. The overlap between adenovirus-suppressed genes and NAFLD-upregulated genes reflects viral manipulation of the same metabolic circuits disrupted in NAFLD.

### Genes Mimicking NAFLD Signature (positive NES, _UP sets)
The gene sets most positively enriched in the NAFLD ranked list reflect processes that *drive* NAFLD:
- Inflammatory signalling (TGF-β, TNF, IL-6 targets)
- Fibrosis / ECM remodelling (TGF-β1 response sets)
- Lipid accumulation / de novo lipogenesis
- Hepatocellular carcinoma progression signatures

These validate our NAFLD signature is capturing the known pathobiology.

---

## Files

### Script
- `scripts/nafld_drug_repositioning.R` — full pipeline

### Plots
| File | Description |
|------|-------------|
| `01_nes_volcano.png` | NES vs −log10(p): all 2749 gene sets, candidates highlighted |
| `02_top_drug_candidates_bar.png` | Top 25 Type-A candidates by NES |
| `03_nafld_mimicking_sets.png` | Top 20 gene sets mimicking NAFLD signature |
| `04_gsea_enrichment_top_candidates.png` | GSEA running-score plots, top 6 candidates |
| `05_bubble_plot.png` | NES × FDR × size bubble chart |
| `06_leading_edge_heatmap.png` | Leading-edge genes shared across top candidates |
| `07_deg_volcano_leading_edge.png` | DEG volcano with leading-edge genes highlighted |

### Results
| File | Description |
|------|-------------|
| `fgsea_all_CGP_results.csv` | Full GSEA results (2,749 gene sets) |
| `top_drug_candidates.csv` | Filtered candidate drug gene sets |
| `drug_summary_paired.csv` | UP/DN gene sets paired per drug |
| `significant_degs.csv` | Input DEGs (padj<0.05, |lFC|>1) |

---

## References

- Subramanian A et al. (2005). Gene set enrichment analysis. *PNAS*, 102(43), 15545–15550.  
- Korotkevich G et al. (2021). Fast gene set enrichment analysis. *bioRxiv*.  
- Lamb J et al. (2006). The Connectivity Map: using gene-expression signatures to connect small molecules, genes, and disease. *Science*, 313(5795), 1929–1935.  
- Morizane Y et al. (2011). c-Jun N-terminal kinase in liver disease. *Clin J Gastroenterol*, 4, 65–72.
