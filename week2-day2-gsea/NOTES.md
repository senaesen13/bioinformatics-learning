# Week 2 — Day 2 — Gene Set Enrichment Analysis with clusterProfiler

> Builds directly on the Week 2 Day 1 DESeq2 results (airway dataset).
> Script: `gsea_analysis.R`

---

## Why Go Beyond Individual Genes?

DESeq2 gives you a list of individual differentially expressed genes. That's useful,
but a single gene rarely acts alone. Biological processes involve coordinated changes
across groups of related genes — a pathway.

Enrichment analysis asks a higher-level question:
> "Do genes belonging to a known pathway change **as a group** in my experiment?"

Two methods covered here:

| Method | Input | Logic |
|--------|-------|-------|
| **ORA** (Over-Representation Analysis) | Significant genes only | Are pathway genes over-represented in my hit list vs chance? |
| **GSEA** (Gene Set Enrichment Analysis) | All genes, ranked | Do pathway genes cluster at the top or bottom of a ranked list? |

---

## Databases Used

### Gene Ontology (GO)
A structured vocabulary of biological terms organised into three categories:

| Ontology | Abbreviation | Example term |
|----------|-------------|--------------|
| Biological Process | BP | "inflammatory response" |
| Molecular Function | MF | "cytokine receptor binding" |
| Cellular Component | CC | "extracellular matrix" |

GO terms are hierarchical — broad terms (e.g. "immune system process") contain
many specific child terms (e.g. "T cell activation"). We use `simplify()` to
remove redundant parent/child terms from our results.

### KEGG
Kyoto Encyclopedia of Genes and Genomes — manually curated pathway maps.
Examples: "Cytokine-cytokine receptor interaction", "MAPK signalling pathway".
More conservative than GO but higher biological confidence.

---

## Step 1 — Re-run DESeq2

The script re-runs the minimal DESeq2 pipeline from Week 2 Day 1 to get
fresh results rather than depending on a saved CSV file.

```r
dds <- DESeqDataSet(airway, design = ~ cell + dex)
dds$dex <- relevel(dds$dex, ref = "untrt")
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)
res        <- results(dds, name = "dex_trt_vs_untrt", alpha = 0.05)
res_shrunk <- lfcShrink(dds, coef = "dex_trt_vs_untrt", type = "apeglm")
```

One extra step compared to Day 1: the **Wald statistic** (`stat` column) is
merged back from the raw results into the shrunken results data frame.
`lfcShrink()` drops the `stat` column, but we need it for GSEA ranking.

```r
res_df <- as.data.frame(res_shrunk) %>%
  left_join(
    as.data.frame(res) %>% dplyr::select(ensembl_id, wald_stat = stat),
    by = "ensembl_id"
  )
```

---

## Step 2 — ID Conversion: Ensembl → Entrez

The airway dataset uses **Ensembl IDs** (e.g. `ENSG00000000003`).

| Function | Accepts |
|----------|---------|
| `enrichGO()` | Ensembl IDs directly |
| `enrichKEGG()` | Entrez IDs only |
| `gseKEGG()` | Entrez IDs only |

We convert using `AnnotationDbi::select()` with `org.Hs.eg.db` (the human
genome annotation package).

```r
id_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = res_df$ensembl_id,
  keytype = "ENSEMBL",
  columns = c("ENTREZID", "SYMBOL")
)
```

About 7% of Ensembl IDs fail to map — this is normal. Some Ensembl IDs
(e.g. pseudogenes, novel transcripts) don't have a corresponding Entrez ID.

---

## Step 3 — Define Gene Lists

### For ORA: a hit list + a universe

```r
sig_ensembl <- filter(res_df, padj < 0.05, abs(log2FoldChange) > 1)
universe    <- res_df$ensembl_id   # all genes we tested
```

The **universe** is important: it defines the background for the statistical test.
Using only the genes you tested (not the whole genome) is more accurate because
genes you couldn't detect shouldn't count as "not enriched".

### For GSEA: a ranked list of all genes

```r
gsea_ranks <- res_df %>%
  arrange(desc(wald_stat)) %>%
  { setNames(.$wald_stat, .$ENTREZID) }
```

We rank by the **Wald statistic** (LFC / SE from DESeq2), not just fold change.
The Wald stat captures both the magnitude *and* the confidence of the change:
- A gene with LFC = 3 but high noise → moderate Wald stat
- A gene with LFC = 1 but very consistent → high Wald stat

Positive values → upregulated in treated. Negative → downregulated.

---

## Part A — ORA (Over-Representation Analysis)

### What it tests

Given:
- **N** total genes in the universe
- **K** genes from pathway X in the universe
- **n** genes in my hit list
- **k** genes from pathway X in my hit list

Is **k** larger than expected by chance?
This is a **hypergeometric test** (like a one-sided Fisher's exact test).

```
p = P(X ≥ k) under hypergeometric(N, K, n)
```

### Limitation

ORA treats all significant genes as equal. It ignores fold change magnitude —
a gene with LFC = 0.01 (just barely over the threshold) counts the same as
LFC = 5. The choice of significance threshold also affects results.

---

## Step 4 — GO ORA

```r
ego <- enrichGO(
  gene          = sig_ensembl,
  universe      = universe_ensembl,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENSEMBL",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05
)
```

**Result: 227 enriched GO-BP terms** (before simplification)

Top hits:
| GO term | Description | Adjusted p |
|---------|-------------|-----------|
| GO:0003013 | circulatory system process | 8.7e-07 |
| GO:0008015 | blood circulation | 3.1e-06 |
| GO:0007169 | receptor tyrosine kinase signalling | 1.9e-05 |
| GO:0001944 | vasculature development | 1.9e-05 |

These make biological sense: dexamethasone is a glucocorticoid that affects
vascular tone and inflammatory signalling in smooth muscle cells.

### simplify()

```r
ego_simplified <- simplify(ego, cutoff = 0.7, by = "p.adjust")
```

GO terms are highly redundant — "blood circulation" and "circulatory system process"
cover almost the same genes. `simplify()` removes terms whose gene sets overlap
above a similarity threshold, keeping the most significant representative.

### Plots produced

| File | What it shows |
|------|--------------|
| `go_ora_dotplot.png` | Top terms: size = gene count, colour = adjusted p |
| `go_ora_barplot.png` | Top terms as horizontal bars |
| `go_ora_cnetplot.png` | Network: GO terms linked to the genes driving them |

---

## Step 5 — KEGG ORA

```r
ekegg <- enrichKEGG(
  gene     = sig_entrez,
  universe = universe_entrez,
  organism = "hsa"   # hsa = Homo sapiens
)
```

**Result: 15 enriched KEGG pathways**

Top hits:
| Pathway | Description | Adjusted p |
|---------|-------------|-----------|
| hsa04060 | Cytokine-cytokine receptor interaction | 0.00026 |
| hsa04923 | Regulation of lipolysis in adipocytes | 0.00026 |
| hsa04750 | Inflammatory mediator regulation of TRP channels | 0.002 |

Cytokine signalling is a known dexamethasone target — it suppresses
pro-inflammatory cytokines (e.g. IL-6, TNF). This is expected.

---

## Part B — GSEA (Gene Set Enrichment Analysis)

### What it tests

GSEA walks down a ranked gene list (best → worst) and asks:
> "As I go through this ranked list, do the genes from pathway X accumulate
> faster than random? Or slower?"

It keeps a **running enrichment score (ES)** that increases when it hits a
pathway gene and decreases when it hits a non-pathway gene.

The **maximum deviation** from zero is the ES. This is then normalised by the
pathway size → **NES (Normalised Enrichment Score)**.

```
NES > 0  →  pathway genes cluster at the top  → upregulated in treated
NES < 0  →  pathway genes cluster at the bottom → downregulated in treated
```

Statistical significance is assessed by **permuting the gene list** many times
and seeing how extreme the real NES is compared to random chance.

### Why GSEA beats ORA

- Uses **all genes** — no arbitrary significance cutoff
- Catches **subtle but consistent** pathway changes (many genes each moving a little)
- Preserves **directionality** — tells you if a pathway is up or down
- Less sensitive to threshold choices

### When ORA beats GSEA

- Faster and simpler to interpret
- Better when you have a small, high-confidence hit list
- Easier to explain to non-statisticians

---

## Step 6 — GO GSEA

```r
gsea_go <- gseGO(
  geneList  = gsea_ranks,
  OrgDb     = org.Hs.eg.db,
  keyType   = "ENTREZID",
  ont       = "BP",
  minGSSize = 15,
  maxGSSize = 500
)
```

**Result: 836 enriched GO-BP terms**
(More than ORA because GSEA uses all genes and captures subtle pathway shifts)

### Plots produced

| File | What it shows |
|------|--------------|
| `go_gsea_dotplot.png` | Activated vs suppressed pathways side by side |
| `go_gsea_enrichmentplot.png` | Classic GSEA mountain plot for the top term |

### Reading the enrichment plot

```
Top panel:    Running enrichment score — the "mountain"
              Peak = where pathway genes are concentrated in the ranked list

Middle panel: Black ticks = where pathway genes fall in the ranked list

Bottom panel: Gradient bar — left = high Wald stat (upregulated), right = low
```

A tall mountain peaking in the first third of the ranked list = strong upregulation.
A valley bottoming out in the last third = strong downregulation.

---

## Step 7 — KEGG GSEA

```r
gsea_kegg <- gseKEGG(
  geneList = gsea_ranks,
  organism = "hsa",
  minGSSize = 15,
  maxGSSize = 500
)
```

**Result: 108 enriched KEGG pathways**

Top results include cytokine signalling and steroid biosynthesis pathways —
consistent with dexamethasone's known mechanism.

---

## Results Summary

| Analysis | Result |
|----------|--------|
| GO ORA   | 227 enriched Biological Process terms |
| KEGG ORA | 15 enriched pathways |
| GO GSEA  | 836 enriched Biological Process terms |
| KEGG GSEA | 108 enriched pathways |

---

## Output Files

| File | Contents |
|------|---------|
| `go_ora_dotplot.png` | GO-BP ORA top terms, dot plot |
| `go_ora_barplot.png` | GO-BP ORA top terms, bar plot |
| `go_ora_cnetplot.png` | Gene-concept network for GO-BP ORA |
| `kegg_ora_dotplot.png` | KEGG ORA pathways, dot plot |
| `go_gsea_dotplot.png` | GO-BP GSEA activated/suppressed split |
| `go_gsea_enrichmentplot.png` | GSEA mountain plot for top GO term |
| `kegg_gsea_dotplot.png` | KEGG GSEA pathways split by direction |
| `kegg_gsea_enrichmentplot.png` | GSEA mountain plot for top KEGG pathway |
| `go_ora_results.csv` | Full GO ORA results table |
| `kegg_ora_results.csv` | Full KEGG ORA results table |
| `go_gsea_results.csv` | Full GO GSEA results table |
| `kegg_gsea_results.csv` | Full KEGG GSEA results table |

---

## Key Concepts Reference

| Concept | Meaning |
|---------|---------|
| **ORA** | Tests if known pathway genes appear more than expected in a hit list |
| **GSEA** | Tests if pathway genes cluster at top or bottom of a ranked list |
| **NES** | Normalised Enrichment Score: positive = pathway upregulated, negative = down |
| **GO** | Gene Ontology — hierarchical vocabulary of biological functions |
| **KEGG** | Curated maps of biological pathways with directional interactions |
| **Universe** | Background gene set for ORA — should be all genes you tested |
| **Wald stat** | DESeq2's ranking metric for GSEA: LFC / SE — captures size and confidence |
| **simplify()** | Removes redundant GO terms with overlapping gene sets |
| **GeneRatio** | Fraction of your significant genes that are in the pathway |
| **BgRatio** | Fraction of the universe that are in the pathway |
| **hypergeometric test** | Statistical test behind ORA — like Fisher's exact test |

---

## Connection to Week 1

In Week 1 Day 4 we ran GSEA using **fgsea** on mouse KEGG pathways with Entrez IDs.
Here we use **clusterProfiler**, which is more widely used, supports GO and KEGG,
handles ID conversion internally, and produces publication-quality plots out of the box.

| | Week 1 Day 4 (fgsea) | Week 2 Day 2 (clusterProfiler) |
|-|---|---|
| Package | fgsea | clusterProfiler |
| Organism | Mouse (*Mus musculus*) | Human (*Homo sapiens*) |
| Databases | KEGG only | GO + KEGG |
| Methods | GSEA only | ORA + GSEA |
| ID type | Entrez | Ensembl → Entrez conversion |
| Plots | Manual ggplot | Built-in dotplot, cnetplot, gseaplot2 |

---

*Week 2, Day 2 — GO and KEGG enrichment on airway DESeq2 results*
