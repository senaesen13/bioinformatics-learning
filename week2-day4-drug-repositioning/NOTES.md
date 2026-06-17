# Week 2 — Day 4 — Drug Repositioning

> **Script:** `scripts/drug_repositioning.R`
> **Run from:** `week2-day4-drug-repositioning/` as `Rscript scripts/drug_repositioning.R`
> **Input:** `../week2-day3-kcl-mouse-deseq2/deseq2_results_MI_vs_sham.csv`
> **Database:** MSigDB C2:CGP (Chemical & Genetic Perturbations) via `msigdbr`

---

## What Is Drug Repositioning?

Drug repositioning (also called drug repurposing) means finding **new therapeutic
uses for drugs that already exist** — drugs that are approved, or have at least
gone through safety testing for some other condition.

The core idea is simple: if a disease makes certain genes go up and other genes
go down, then a drug that does the **opposite** — makes the same genes go down and
the same genes go up — might be able to reverse or counteract the disease.

This is attractive for several reasons:
- Approved drugs already have known safety profiles — you skip early-phase toxicity
  testing, which is expensive and slow
- Phase I and II clinical trials may already be done for a different indication
- Development timelines compress from ~15 years to 3–5 years
- Much cheaper than developing a new drug from scratch (~$2 billion average)

A famous example: thalidomide was originally withdrawn for causing birth defects, but
was later repositioned as a treatment for multiple myeloma. Sildenafil (Viagra) was
originally developed for angina; its cardiovascular effects were noticed incidentally.

---

## Where Our Input Data Came From

We used the DESeq2 differential expression results from **Week 2 Day 3**, which
analysed the KCL mouse myocardial infarction (MI) dataset.

| Detail | Value |
|--------|-------|
| Dataset | GSE114134 — Mouse cardiac tissue, MI vs sham surgery |
| Species | *Mus musculus* (mouse) |
| Comparison | MI (heart attack) vs sham-operated controls |
| Total genes tested | 18,074 |
| Upregulated in MI | 1,605 genes (padj < 0.05, |log2FC| > 1) |
| Downregulated in MI | 1,942 genes (padj < 0.05, |log2FC| > 1) |
| File used | `deseq2_results_MI_vs_sham.csv` |

The DESeq2 results give us the "MI disease signature" — a ranked list of genes
that are differentially expressed when the heart has a heart attack. This
signature is what we query against drug databases.

---

## Pipeline Overview

```
Week 2 Day 3 DESeq2 results (mouse, 18,074 genes)
  │
  ▼ STEP 2: Load all genes + log2FC values
  │
  ▼ STEP 3: BioMart — convert mouse Ensembl IDs → human Entrez IDs
  │         (drug databases use human genes)
  │
  ▼ STEP 4: Build gene lists
  │         • GSEA: all 15,098 genes ranked by log2FC (most up → most down)
  │         • ORA: top 200 up genes + top 200 down genes (by padj)
  │
  ▼ STEP 5: MSigDB C2:CGP — load 428 verified drug gene sets
  │         (filtered to those with a PubChem compound ID = real drugs)
  │
  ├─▶ STEP 6: GSEA — score each drug gene set against the ranked MI list
  │           NES < 0 = drug opposes MI (repositioning candidate)
  │           NES > 0 = drug mimics MI (avoid)
  │
  └─▶ STEP 7: ORA — test overlap between MI genes and drug signatures
              MI-up genes vs drug-DN sets (drugs that knock down MI-activated genes)
              MI-dn genes vs drug-UP sets (drugs that restore MI-suppressed genes)
```

---

## Step-by-Step Explanation

### Step 1 — Load Libraries

**What it does:** Loads the R packages needed for the analysis. No new installations
are required — everything was already available.

**Packages used:**
- `clusterProfiler` — runs GSEA and ORA enrichment analyses
- `msigdbr` — provides access to the MSigDB gene set collections
- `biomaRt` — queries the Ensembl database to convert gene IDs between species
- `org.Hs.eg.db` — local human gene annotation database (Ensembl ↔ Entrez ID mapping)
- `AnnotationDbi` — tools for working with annotation databases
- `ggplot2` / `dplyr` — plotting and data manipulation

---

### Step 2 — Load DESeq2 Results

**What it does:** Reads the CSV file from Week 2 Day 3 into R, removing any genes
with missing fold-change or p-values.

**Why we need all genes, not just significant ones:** GSEA (Step 6) works by scoring
the *entire* ranked transcriptome, not just a list of significant genes. It looks at
whether a drug's gene set tends to cluster at the top or bottom of the ranked list.
If we only provided significant genes, we would lose the information about direction
and effect size that makes GSEA powerful.

**What the columns mean:**
- `gene_id` — mouse Ensembl gene ID (ENSMUSG...)
- `log2FoldChange` — how much higher/lower the gene is in MI vs sham (positive = up in MI)
- `padj` — Benjamini-Hochberg adjusted p-value (our confidence that the change is real)
- `sig` — "Up", "Down", or "NS" (not significant)
- `symbol` — the gene's common name (e.g., Rtn4, Asb2)

---

### Step 3 — Mouse → Human Ortholog Conversion

**What it does:** Converts mouse gene IDs to their human equivalents.

**Why this is necessary:** The drug databases (MSigDB C2:CGP) contain gene sets
based on experiments done in *human* cell lines. Our DESeq2 results are from *mouse*
heart tissue. To query the drug database, we need to translate the mouse genes into
the closest human counterparts — called **orthologues**.

**How it works:**
1. We connect to the Ensembl BioMart database (trying multiple mirrors in case one
   is slow or down)
2. For each mouse Ensembl gene ID, we ask: "What is the human gene that is
   evolutionarily equivalent to this mouse gene?"
3. We keep only **high-confidence 1:1 orthologues** — mouse genes that map to
   exactly one human gene with high evolutionary confidence
4. We then convert the human Ensembl IDs to **Entrez IDs** (NCBI's numeric gene
   identifiers), because MSigDB uses Entrez IDs internally

**What's lost:** About 15–20% of mouse genes have no clear 1:1 human orthologue.
These are dropped. The heart is a specialised tissue; some mouse cardiac-specific
genes may not have well-characterised human equivalents.

**Numbers:** Starting from 18,074 mouse genes, we end up with 15,098 genes that
have a high-confidence human orthologue with a valid Entrez ID.

---

### Step 4 — Build Gene Lists

**What it does:** Prepares two different representations of the data for two
different analyses.

**For GSEA (ranked list):**
- All 15,098 human-converted genes, sorted from most upregulated in MI (top) to
  most downregulated (bottom), by log2FoldChange
- Genes most upregulated in MI (like fibrotic, inflammatory genes) sit at the top
- Genes most downregulated in MI (like fatty acid oxidation genes) sit at the bottom

**For ORA (gene sets):**
- Top 200 most significantly upregulated genes in MI (by adjusted p-value)
- Top 200 most significantly downregulated genes in MI
- 200 is a practical balance: enough to overlap with drug signatures, not so many
  that the signal gets diluted

---

### Step 5 — Load and Filter Drug Gene Sets (MSigDB C2:CGP)

**What it does:** Loads the Chemical and Genetic Perturbations (C2:CGP) collection
from MSigDB and filters it down to only actual drug/compound experiments.

**What MSigDB C2:CGP is:**
MSigDB (Molecular Signatures Database) is a curated collection of gene sets maintained
by the Broad Institute. The C2:CGP subcollection contains 3,555 gene sets derived from
experiments where researchers treated cells with drugs, chemicals, or genetic
manipulations and measured which genes changed.

**The filtering problem:** C2:CGP mixes drug experiments with genetic perturbations
(gene knockouts, overexpression) and tissue comparisons. If we used all 3,555 sets,
the "drug names" would actually include things like "OPA1 overexpression" or
"heart atrium vs ventricle comparison" — not drugs.

**How we filter:** We keep only gene sets whose MSigDB description contains a
`[PubChem=...]` identifier. PubChem is NCBI's chemical compound database, and
MSigDB adds these IDs specifically when a gene set comes from a chemical compound
experiment. This is MSigDB's own annotation for verified drug perturbations.

**Result:** 428 gene sets covering 122 distinct compounds — the actual drugs.

**Gene set naming:** Each drug typically has two gene sets:
- `AUTHORNAME_DRUG_CONTEXT_UP` — genes that go UP after drug treatment
- `AUTHORNAME_DRUG_CONTEXT_DN` — genes that go DOWN after drug treatment

---

### Step 6 — GSEA: Find Drug Programs Anti-Correlated with MI

**What it does:** Runs Gene Set Enrichment Analysis comparing the MI gene ranking
against each drug's gene set to compute how similar or opposite each drug's effect
is to MI.

**The concept of GSEA:**
Imagine our 15,098 genes arranged in a line, from most upregulated in MI (left) to
most downregulated in MI (right). For each drug, we ask: "Do the genes that go UP
with this drug tend to cluster at the left end of our MI line, or the right end?"

- If a drug's UP genes cluster on the LEFT (most upregulated in MI): this drug
  mimics MI. NES > 0. Avoid.
- If a drug's UP genes cluster on the RIGHT (most downregulated in MI): this drug
  does the OPPOSITE of MI. NES < 0. Repositioning candidate.
- If the drug's genes are scattered randomly: no relationship. NES ≈ 0.

**The NES (Normalized Enrichment Score):**
- A number between roughly -3 and +3
- NES < 0 = the drug's gene program is ANTI-correlated with MI
- The more negative, the stronger the anti-MI effect
- p.adjust < 0.05 = the enrichment is statistically significant (not by chance)

The analysis tests all 428 drug gene sets and completes in ~13 seconds.

---

### Step 7 — ORA: Connectivity Analysis

**What it does:** A complementary analysis that tests specific directional overlaps
between MI genes and drug signatures.

**The connectivity map logic** (borrowed from the Broad Institute's CMAP approach):
- Genes that are upregulated in MI: a good drug should KNOCK THESE DOWN
  → Test: do our MI-up genes significantly overlap with a drug's DOWN gene set?
- Genes that are downregulated in MI: a good drug should RESTORE THESE
  → Test: do our MI-down genes significantly overlap with a drug's UP gene set?

**How enricher() works:** It uses the hypergeometric test — essentially asking
"given that there are X genes in both our MI list and the drug's gene set,
would this many overlapping genes happen by chance?"

**Universe:** We use all 15,098 human-converted genes as the background universe,
so the test is calibrated to what we could have detected, not all human genes.

**ORA result:** One significant hit — estradiol (the DUTERTRE_ESTRADIOL_RESPONSE_24HR_DN
gene set). 20 of our top 200 MI-upregulated genes are in the set of genes that
estradiol knocks down. This independently confirms the GSEA estradiol finding.

---

## Final Results: The 6 Drug Candidates

All 6 candidates have GSEA p.adjust < 0.05, meaning the anti-MI enrichment is
statistically significant. They are ranked by mean NES (most negative = strongest
anti-MI signal).

| Rank | Compound | Mean NES | adj. p | Gene Sets | What it is | Why it makes biological sense |
|------|----------|----------|--------|-----------|------------|-------------------------------|
| 1 | **WY14643** | -2.05 | 0.0004 | 1 | Synthetic PPARα agonist (research compound) | After MI, the heart loses its ability to oxidise fatty acids and switches to glucose. PPARα is the master regulator of cardiac fatty acid oxidation. WY14643 activates PPARα and restores the metabolic gene programme that MI suppresses. Strong biological rationale. |
| 2 | **Troglitazone** | -1.81 | 0.0006 | 3 | PPARγ agonist (thiazolidinedione class) | First-generation insulin sensitiser (withdrawn for liver toxicity). PPARγ activation has anti-inflammatory and cardioprotective effects. Three independent gene sets agree on the anti-MI signal, increasing confidence. |
| 3 | **CL-387785** | -1.45 | 0.041 | 2 | Irreversible EGFR/ErbB inhibitor | ErbB/EGFR signalling is activated after MI and drives pathological cardiac remodelling and fibrosis. Inhibiting ErbB opposes this remodelling programme. |
| 4 | **TPA** | -1.35 | 0.007 | 4 | Phorbol ester (PKC activator) | A carcinogen and research tool, not therapeutic. Included here because PKC (protein kinase C) activation has complex roles in ischaemic preconditioning — brief PKC activation can paradoxically protect the heart. Four gene sets support this signal. |
| 5 | **Estradiol** | -1.25 | 0.003 | 6 | 17β-oestradiol (female sex hormone) | The best-documented cardioprotective hormone. Pre-menopausal women have lower MI rates than age-matched men; this protection diminishes after menopause. Estradiol reduces inflammation, improves endothelial function, and inhibits pathological hypertrophy. Six gene sets independently confirm this signal. Also confirmed by ORA (20 MI-upregulated genes are knocked down by estradiol). |
| 6 | **Progesterone** | -1.15 | 0.0004 | 11 | Female sex hormone | The most gene-set support of all candidates (11 sets). Progesterone has established cardioprotective effects including anti-inflammatory and anti-fibrotic actions. Its combination with estradiol is well-studied in the context of hormone replacement therapy and cardiovascular risk. |

### Interpreting NES and Confidence

- **WY14643** has the strongest single-set NES (-2.05) and very low p-value. High
  confidence that PPARα agonism opposes MI at the transcriptomic level.
- **Progesterone** has the most independent evidence (11 gene sets, p = 0.0004).
  When many different studies in different cell lines all point to the same direction,
  that's strong convergent evidence.
- **Estradiol** is confirmed by both GSEA and ORA — two independent statistical
  approaches agreeing on the same candidate increases confidence further.
- **TPA** is a carcinogen and not a therapeutic candidate. It appears here because
  PKC signalling overlaps mechanistically with cardioprotective pathways, but this
  illustrates that drug repositioning results need biological interpretation, not
  just statistical ranking.

---

## What This Analysis Does NOT Guarantee

Drug repositioning generates **hypotheses**, not proven therapies. Significant caveats:

| Limitation | Why it matters |
|------------|----------------|
| Cell line context | MSigDB drug signatures come from cancer cell lines (MCF7, PC3, HL60), not cardiac tissue. The same drug may behave differently in cardiomyocytes. |
| Mouse → human conversion | ~15–20% of mouse MI genes have no direct human orthologue and are excluded. Some mouse-specific biology is missed. |
| Gene expression ≠ disease reversal | A drug reversing gene expression patterns does not guarantee it reverses disease outcome. |
| Single concentration | Drug signatures in MSigDB represent one dose in one experimental context. Dose-dependence and pharmacokinetics are ignored. |
| Off-target effects | Low NES means anti-MI transcriptomic effect; it says nothing about toxicity, drug interactions, or off-target effects. |

**Validation steps that would follow a real repositioning hit:**
1. Test in human cardiac cell models (iPSC-derived cardiomyocytes)
2. In vivo mouse MI model with drug treatment
3. Retrospective analysis of clinical databases (patients already on estradiol or PPARγ agonists — do they have better MI outcomes?)
4. Prospective clinical trial

---

## Output Files

| File | What it contains |
|------|-----------------|
| `results/gsea_drug_candidates.csv` | One row per compound: mean NES, best adj. p, number of supporting gene sets |
| `results/gsea_anti_mi_hits.csv` | All anti-MI gene set hits (NES < 0) with full GSEA statistics |
| `results/gsea_all_drug_sets.csv` | Complete GSEA results for all 413 tested drug gene sets |
| `results/drug_gene_set_metadata.csv` | Metadata for all PubChem-verified drug gene sets |
| `results/ora_MI_up_vs_drug_DN.csv` | ORA: MI-upregulated genes vs drug-downregulating sets |
| `plots/gsea_top_drugs.png` | Lollipop chart of the top significant anti-MI drug candidates |
| `plots/gsea_nes_distribution.png` | Histogram of NES scores across all drug gene sets |
| `plots/ora_connectivity.png` | Bubble plot of ORA connectivity hits |

---

## Connection to Previous Weeks

| Week | Day | What was done | Used here |
|------|-----|---------------|-----------|
| Week 1 | Day 4 | KCL workshop: DESeq2 on mouse MI | Context |
| Week 2 | Day 3 | Clean DESeq2 re-analysis of mouse MI (GSE114134) | DESeq2 results CSV |
| Week 2 | Day 4 | This analysis — drug repositioning | — |

The DESeq2 → drug repositioning pipeline represents a real computational drug
discovery workflow used in academic and industry research. Versions of this approach
have identified cardiac drug candidates including repurposed metabolic drugs and
hormonal interventions that are now in clinical investigation.

---

*Week 2, Day 4 — Drug repositioning on mouse MI expression signature using
clusterProfiler + MSigDB C2:CGP (Chemical & Genetic Perturbations)*
