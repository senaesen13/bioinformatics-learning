# Week 1 — Day 4 — Bulk RNA-seq Pipeline (KCL Workshop)

> Hands-on pipeline following the KCL Module 2 Transcriptomics workshop.
> Workshop code adapted from: KCL Module 2 Transcriptomics
> Original source: https://github.com/sysmedicine/KCLModule2_2022
> Original author: Yang Hong

---

## The Experiment

Real mouse heart attack (Myocardial Infarction) data.
Goal: Find which genes change after a heart attack.

8 samples:
- sham_1dMI_1 / sham_1dMI_2  → healthy mouse, day 1
- MI_1dMI_1   / MI_1dMI_2    → heart attack mouse, day 1
- sham_3dMI_1 / sham_3dMI_2  → healthy mouse, day 3
- MI_3dMI_1   / MI_3dMI_2    → heart attack mouse, day 3

---

## The Full Pipeline

```
Kallisto output files (abundance.h5)
  └─ tximport → count matrix (21,906 genes x 8 samples)
       └─ filter low counts → 19,191 genes remaining
            └─ DESeq2 → differential expression results
                 └─ 13 significant genes identified
                      └─ GSEA → biological pathways enriched
```

---

## Step 1 — Set Working Directory

```r
setwd("/Users/senaesen/Desktop/bioinfo-learning/KCLModule2_2022/Transcriptomics")
list.files()
```

---

## Step 2 — Load Libraries

```r
library(tximport)
library(DESeq2)
library(ggplot2)
library(dplyr)
library(EnrichmentBrowser)
library(fgsea)
library(AnnotationDbi)
library(org.Mm.eg.db)
```

---

## Step 3 — Load Protein Coding Genes

```r
load("proteinCodingGenes.Rda")
head(proteinCodingGenes)
```

102,245 mouse protein-coding genes with Ensembl IDs and gene names.

---

## Step 4 — Load Kallisto Files

```r
sra_accession = rev(dir("abundances/"))
files = sort(file.path("abundances", sra_accession, "abundance.h5"))
names(files) = c("sham_1dMI_1", "sham_1dMI_2",
                 "MI_1dMI_1",   "MI_1dMI_2",
                 "sham_3dMI_1", "sham_3dMI_2",
                 "MI_3dMI_1",   "MI_3dMI_2")

txi.kallisto <- tximport(files,
                         type = "kallisto",
                         txOut = F,
                         tx2gene = proteinCodingGenes,
                         ignoreTxVersion = T)
```

---

## Step 5 — Create Matrices

```r
mat.count  <- txi.kallisto$counts     # raw counts
mat.tpm    <- txi.kallisto$abundance  # TPM normalised
dim(mat.count)  # 21906 x 8

mat.design <- data.frame(
  SRA_accession = sra_accession,
  Sample_name   = colnames(mat.count),
  Treatment     = sapply(strsplit(colnames(mat.count), split = "_"), function(x) x[1]),
  Time          = sapply(strsplit(colnames(mat.count), split = "_"), function(x) x[2])
)
```

---

## Step 6 — Filter Low Quality Genes

```r
keep       <- rowSums(mat.count) > 10
mat.count  <- mat.count[keep, ]
mat.tpm    <- mat.tpm[keep, ]
dim(mat.count)  # 19191 x 8
```

Removed 2,715 genes with fewer than 10 total reads — unreliable noise.

---

## Step 7 — DESeq2 Analysis

```r
# Create DESeq2 object
dset <- DESeqDataSetFromMatrix(
  countData = round(mat.count),
  colData   = mat.design,
  design    = ~ Treatment + Time
)

# Run DESeq2
ds2 <- DESeq(dset, fitType = "parametric", parallel = T)

# Extract results (MI vs sham)
res = DESeq2::results(ds2,
                      cooksCutoff = T,
                      alpha = 0.01,
                      lfcThreshold = 0,
                      contrast = c("Treatment", "MI", "sham"))
res = data.frame(res)
res$padj[is.na(res$padj)] = 1

# Get significant genes (padj < 0.01)
res.interest <- res[res$padj < 0.01 & !is.na(res$padj), ]

# Add real gene names
res.interest$gene_name <- proteinCodingGenes$external_gene_name[
  match(rownames(res.interest), proteinCodingGenes$ensembl_gene_id)]

View(res.interest)
```

### Results columns explained

| Column | Meaning |
|--------|---------|
| baseMean | Average expression across all samples |
| log2FoldChange | How much gene changed (+ up in MI, - down in MI) |
| padj | Adjusted p-value (corrected for multiple testing) |

### 13 significant genes found

| Gene | Direction | Known role |
|------|-----------|-----------|
| Spp1 | Up | Inflammation, cardiac remodelling |
| Cthrc1 | Up | Heart tissue repair |
| Hspa1a | Down | Heat shock protein |
| Hspa1b | Down | Heat shock protein |
| Ncapg | Up | Cell division |

---

## Step 8 — GSEA (Gene Set Enrichment Analysis)

### Why GSEA?

Individual genes → what changed
GSEA → WHY it changed (which biological processes)

GSEA ranks all genes by DESeq2 stat score and checks if genes in known pathways cluster at the top (upregulated) or bottom (downregulated).

```r
# Load KEGG mouse pathways (367 pathways)
kegg.gs <- getGenesets(org = "mmu", db = "kegg",
                       cache = T, return.type = "list")

# Convert pathway IDs to Ensembl IDs
kegg.gs.list = list()
count = 0
for(pathway in kegg.gs){
  count = count + 1
  kegg.gs.list[[names(kegg.gs[count])]] <- AnnotationDbi::select(
    org.Mm.eg.db, keys = pathway,
    keytype = "ENTREZID", columns = "ENSEMBL")$ENSEMBL
}

# Prepare ranked gene list
kegg_input = sort(deframe(cbind(Symbol = rownames(res), res["stat"])))
kegg_input = kegg_input[is.finite(kegg_input)]
names(kegg_input) = proteinCodingGenes[
  match(names(kegg_input), proteinCodingGenes$ensembl_gene_id), 3]
kegg_input = kegg_input[unique(names(kegg_input))]

# Run GSEA
fgseaRes = fgsea(kegg.gs.list, kegg_input)

# View significant pathways
fgseaResTidy <- fgseaRes %>% as_tibble() %>% arrange(desc(NES))
fgseaResTidy %>%
  filter(pval < 0.05) %>%
  arrange(pval) %>%
  dplyr::select(pathway, NES, padj)
```

NES = Normalized Enrichment Score
- Positive = pathway upregulated after heart attack
- Negative = pathway downregulated after heart attack

---

## Key Concepts

| Concept | Meaning |
|---------|---------|
| Count matrix | Genes x samples table with read counts |
| TPM | Normalised expression — transcripts per million |
| Design matrix | Table describing each sample's conditions |
| log2FoldChange | How much a gene changed (log2 scale) |
| padj | Adjusted p-value — corrected for multiple testing |
| NES | Normalized Enrichment Score from GSEA |
| KEGG | Database of known biological pathways |

---

## Code Credit

Workshop code adapted from KCL Module 2 Transcriptomics (2022)
Source: https://github.com/sysmedicine/KCLModule2_2022
Original author: Yang Hong
Usage: personal learning and educational purposes only

---

*Week 1, Day 4 — bulk RNA-seq pipeline on mouse heart attack data*
