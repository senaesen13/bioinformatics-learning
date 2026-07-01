# Week 4 Day 4 — Drug Repositioning for NAFLD

## What is drug repositioning and why did we do it?

Drug repositioning means finding new uses for drugs that already exist. Instead of inventing a brand-new NAFLD drug (which takes 15 years and billions of dollars), we asked: are there approved drugs whose effect on gene expression looks like the *opposite* of NAFLD? If a drug reverses the changes NAFLD causes in the liver, it might treat it.

## What data did we use?

- **Input:** DEGs from Day 1 bulk RNA-seq (GSE162694, Normal vs NAFLD liver biopsies) — 597 significant genes, 581 upregulated in NAFLD
- **Drug database:** MSigDB C2 CGP — 3,555 gene sets, each describing which genes go up or down in response to a drug or chemical
- **Tool:** fgsea (fast GSEA) in R, run against all drug gene sets

## How does it work?

We ranked all ~21,000 genes from most upregulated to most downregulated in NAFLD. Then, for each drug in the database, we asked: are the genes that drug *turns off* clustered near the top of our NAFLD list? If yes, that drug is turning off the same genes NAFLD turns on — so the drug opposes the disease. In GSEA language, this shows up as a positive NES (normalised enrichment score) for the drug's "\_DN" gene set. The higher the NES, the stronger the opposition to the NAFLD signature.

## Key results

| Drug / Perturbation | NES | What it does |
|---|---|---|
| JNK signalling inhibition | 1.46 | Blocking JNK (a stress kinase) suppresses liver inflammation — JNK inhibitors are already being tested in NASH trials |
| Docetaxel (taxane) | 1.48 | A chemotherapy drug that shuts down cell division genes, many of which are abnormally active in NAFLD liver cells |
| HGF signalling | 1.44 | Hepatocyte growth factor protects liver cells and calms the pro-inflammatory genes elevated in NAFLD |
| Adipocyte differentiation genes | 1.45 | Genes switched off during fat cell formation overlap with genes overactive in NAFLD, pointing to shared metabolic circuitry |
| Aging-related genes (rat) | 1.44 | Genes that decline with liver aging are the same genes suppressed in NAFLD, confirming NAFLD accelerates liver ageing |

## What this means

The NAFLD liver has hundreds of genes switched on that shouldn't be — mainly inflammation, fibrosis, and cell stress genes. These results suggest that drugs blocking JNK signalling or HGF-pathway activators could potentially reverse that pattern, and they give us a shortlist of approved drugs worth testing in NAFLD cell or animal models.

## References

- Suppli MP et al. (2021). Hepatic transcriptome signatures in patients with varying degrees of nonalcoholic fatty liver disease compared with healthy normal-weight individuals. *Am J Physiol Gastrointest Liver Physiol*, 322(4), G439–G452. GEO: GSE162694.
- Liberzon A et al. (2015). The Molecular Signatures Database (MSigDB) hallmark gene set collection. *Cell Syst*, 1(6), 417–425.
