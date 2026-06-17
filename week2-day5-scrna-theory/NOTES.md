# Week 2 — Day 5 — Single-Cell RNA Sequencing: Theory

---

## Table of Contents

1. [What Is scRNA-seq?](#1-what-is-scrna-seq)
2. [Why Does It Exist? The Problem with Bulk RNA-seq](#2-why-does-it-exist-the-problem-with-bulk-rna-seq)
3. [How It Works: Step by Step](#3-how-it-works-step-by-step)
   - 3.1 Tissue Dissociation
   - 3.2 Cell Capture and Barcoding
   - 3.3 UMIs (Unique Molecular Identifiers)
   - 3.4 Library Preparation and Sequencing
4. [Computational Pipeline: From Reads to Biology](#4-computational-pipeline-from-reads-to-biology)
   - 4.1 Alignment and Count Matrix Generation
   - 4.2 Quality Control (QC)
   - 4.3 Normalisation
   - 4.4 Feature Selection and Dimensionality Reduction
   - 4.5 Clustering
   - 4.6 Cell Type Annotation
5. [Key Terms Glossary](#5-key-terms-glossary)
6. [How scRNA-seq Connects to Our Bulk RNA-seq Work](#6-how-scrna-seq-connects-to-our-bulk-rna-seq-work)
7. [Applications in Pharma and Drug Discovery](#7-applications-in-pharma-and-drug-discovery)
8. [Limitations and Caveats](#8-limitations-and-caveats)

---

## 1. What Is scRNA-seq?

**Single-cell RNA sequencing (scRNA-seq)** is a technology that measures gene
expression in **individual cells** — one by one — rather than in bulk tissue.

For every single cell in a sample, scRNA-seq tells you which genes that cell is
expressing and at what level. The result is a huge table: rows are genes, columns
are cells, and every number represents how many transcripts of that gene were
detected in that particular cell.

A typical experiment might profile **5,000 to 50,000 cells** in a single run,
producing a dataset with dimensions like 30,000 genes × 20,000 cells.

### The Central Question It Answers

> "What is every cell in this tissue doing — and is each cell the same, or are
> there distinct populations doing different things?"

This question cannot be answered with bulk RNA-seq. With scRNA-seq, you can:
- Discover rare cell types that represent 0.1% of a tissue
- Show that two cells sitting next to each other can have completely different
  gene expression programmes
- Track how cells change their identity over time (pseudotime / trajectory analysis)
- Compare the composition of a tumour or diseased tissue to healthy tissue
  at cellular resolution

---

## 2. Why Does It Exist? The Problem with Bulk RNA-seq

### What bulk RNA-seq measures

In bulk RNA-seq (like our Week 2 analysis of mouse MI), you grind up a piece of
tissue and sequence everything in it together. The result is an **average**
across all cells in that tissue.

Imagine a heart sample containing:
- 60% cardiomyocytes (the beating muscle cells)
- 20% fibroblasts (structural support, scar tissue)
- 10% endothelial cells (blood vessel lining)
- 5% macrophages (immune cells)
- 5% other

The bulk RNA-seq result you get is the weighted average of all of these. If a gene
is upregulated after MI, you cannot tell whether:
- ALL cell types upregulated it a little, or
- Just macrophages massively upregulated it (dominating the signal), while
  cardiomyocytes did nothing

### The averaging problem in disease

This matters enormously in disease research. In myocardial infarction:
- The **fibrotic response** (scar formation) happens in fibroblasts
- The **inflammatory response** happens in macrophages and immune cells
- The **loss of contractility** is a cardiomyocyte-specific problem
- The **angiogenic response** (forming new blood vessels) happens in
  endothelial cells

Bulk RNA-seq mixes all four signals together. You see that *something* changed,
but not *which cells* are responsible.

scRNA-seq unmixes this. It lets you ask: "In MI vs sham hearts, which specific
cell type drives the fibrotic gene programme — is it the same fibroblasts, or
do macrophages also contribute?"

### The heterogeneity problem

Even within a single cell type, individual cells are not identical. A tumour
contains:
- Cells responding to treatment
- Cells that are treatment-resistant
- Cells in different phases of the cell cycle
- Cancer stem cells vs differentiated cancer cells

Bulk RNA-seq tells you the tumour average. scRNA-seq reveals this internal
heterogeneity — which is why the resistant subpopulation is there even before
treatment starts.

### Brief historical context

| Year | Milestone |
|------|-----------|
| 2009 | First scRNA-seq paper (Tang et al.) — just 1 cell at a time, very low throughput |
| 2013 | CEL-seq and other early platforms — tens to hundreds of cells |
| 2015 | Drop-seq (Macosko et al.) — thousands of cells using droplet microfluidics |
| 2015 | 10x Genomics founded, later commercialising a droplet platform |
| 2017 | 10x Chromium becomes standard — tens of thousands of cells per run |
| 2018 | Human Cell Atlas project begins — map every cell type in the human body |
| 2022+ | Spatial transcriptomics — scRNA-seq with physical location in tissue |

---

## 3. How It Works: Step by Step

### 3.1 Tissue Dissociation

**What happens:** The tissue (e.g., a heart, tumour biopsy, blood sample) is
broken apart into a **single-cell suspension** — a liquid containing individual,
separated cells.

**How it's done:**
- Mechanical dissociation: physically chopping/mincing the tissue
- Enzymatic dissociation: digesting the extracellular matrix proteins that hold
  cells together (collagenase, trypsin, papain, etc.)
- Result: a suspension of single cells in buffer solution

**Why this step is critical and problematic:**

This is arguably the biggest source of **biological artefacts** in scRNA-seq.
The dissociation process is stressful — cells are yanked away from their
neighbours, deprived of their normal signals, and exposed to enzymes. In
response, many cells immediately activate stress response genes (heat shock
proteins, immediate-early genes like *FOS*, *JUN*, *EGR1*).

This means some of what you see in the final data is **dissociation artefact**,
not the biology you wanted to capture. Researchers try to minimise this with:
- Cold dissociation protocols (slower, less stressful)
- Short dissociation times
- Transcriptional inhibitors (actinomycin D) to freeze the transcriptome
- Computational methods to identify and flag stress-activated cells

**Cell viability matters:** Dead cells leak their RNA into the surrounding
solution, which gets captured by other cells' barcodes (ambient RNA contamination).
Libraries like SoupX and CellBender correct for this computationally.

Some cell types don't survive dissociation well:
- **Neurons** are particularly fragile (large, long processes)
- **Adipocytes** (fat cells) are so full of lipid they float and clog equipment
- **Cardiomyocytes** (heart muscle cells) are large and multinucleated — they're
  frequently lost or only partially captured

For these reasons, **single-nucleus RNA-seq (snRNA-seq)** has become popular for
frozen tissue and for cell types that don't survive dissociation: instead of
isolating whole cells, you isolate just their nuclei. Nuclei are more robust,
but you lose cytoplasmic RNA (which includes more mature, processed transcripts).

---

### 3.2 Cell Capture and Barcoding

This is the step that makes scRNA-seq possible — assigning a unique identity tag
to every cell.

**The core problem:** How do you measure gene expression in thousands of cells
simultaneously while keeping each cell's information separate?

**The solution: cell barcodes**

Each cell is tagged with a short DNA sequence called a **cell barcode** — a
unique identifier that gets attached to all RNA from that cell. Every transcript
sequenced later carries this barcode, so the sequencer knows which cell it came from.

#### The 10x Genomics Chromium approach (dominant platform)

The most widely used method uses **droplet microfluidics**:

1. The single-cell suspension is loaded into a microfluidic chip alongside
   thousands of tiny gel beads, each coated with millions of DNA oligonucleotides
2. The chip creates tiny oil droplets (called **GEMs** — Gel Beads in Emulsion)
   at ~2,000 droplets per second
3. The process is calibrated so that, on average, **each droplet captures one
   cell and one bead** — though in practice ~10–20% of droplets capture a cell
4. Inside each droplet:
   - The cell is lysed (broken open), releasing its mRNA
   - The mRNA binds to the oligonucleotides on the bead
   - These oligos contain:
     - A poly-T sequence to capture the poly-A tail on mRNA
     - A **UMI** (unique molecular identifier — see 3.3)
     - A **cell barcode** (unique to that bead/droplet)
     - A sequencing primer

5. The droplets are broken, and reverse transcription converts mRNA → cDNA
6. All cDNA is amplified and sequenced together (they're distinguishable by barcode)

**Result:** For every sequenced read, you know:
- Which cell it came from (cell barcode)
- Which gene it came from (by aligning the sequence to the genome)
- Whether it was a real molecule or a PCR duplicate (UMI — see below)

#### Doublets: when two cells end up in one droplet

If two cells land in the same droplet, their transcripts get the same barcode —
they appear as one "super-cell" in the data. These are called **doublets** (or
multiplets). They are:
- ~1–2% of captured "cells" at typical loading densities
- Higher loading densities → more doublets
- Detected computationally (doublets look like hybrid cell types with unusually
  high gene counts) using tools like DoubletFinder or Scrublet

---

### 3.3 UMIs (Unique Molecular Identifiers)

**The amplification problem:** PCR amplification is necessary to generate enough
material to sequence. But PCR is not perfectly even — some molecules get amplified
50× while others get amplified 5×. If you just count reads, you'll confuse
"this gene had many reads" with "this particular transcript was PCR-amplified many times".

**The UMI solution:**

Before PCR amplification, every individual mRNA molecule gets tagged with a
**UMI** — a short random DNA sequence (typically 10–12 nucleotides) that is
essentially unique to that molecule.

After sequencing, you count **UMIs, not reads**:
- If the same gene in the same cell has 20 reads but only 3 different UMIs,
  those 20 reads all came from 3 original molecules (with PCR generating copies)
- You count: 3 transcripts, not 20

**UMIs are what the count matrix contains.** When you see a number like "47" in a
scRNA-seq count matrix, it means 47 distinct mRNA molecules of that gene were
captured in that cell — not 47 reads.

**Why this matters:**
- Makes quantification more accurate (removes PCR bias)
- Allows proper statistical modelling (UMI counts follow negative binomial distributions)
- Enables meaningful comparison between cells

---

### 3.4 Library Preparation and Sequencing

After droplet breaking and reverse transcription:

1. **cDNA amplification** — PCR amplifies all cDNA (now distinguishable by UMI+barcode)
2. **Fragmentation** — cDNA is fragmented into ~300–500 bp pieces
3. **Sequencing adapters** are ligated to the fragments
4. **Sequencing** — typically on Illumina instruments

**What is sequenced:** In most 10x protocols (3' end capture):
- Read 1: the cell barcode + UMI (28 bp)
- Read 2: the cDNA sequence (used to identify the gene, 90–150 bp)

**Sequencing depth:** Each cell typically gets ~2,000–5,000 reads (often called
"reads per cell" or "sequencing depth"). This is far less than bulk RNA-seq
(which typically uses 20–50 million reads per sample), because depth is spread
across thousands of cells.

The consequence: scRNA-seq data is **sparse**. Most genes in most cells will
have 0 UMI counts — not because the gene isn't expressed, but because the
transcript wasn't captured. This is called **dropout**.

---

## 4. Computational Pipeline: From Reads to Biology

### 4.1 Alignment and Count Matrix Generation

**Tools:** Cell Ranger (10x Genomics), STARsolo, Salmon/Alevin, Kallisto|bustools

**What happens:**
1. Raw sequencing reads (FASTQ files) are aligned to a reference genome
2. Cell barcodes are identified (corrected against a known whitelist of valid barcodes)
3. UMIs are counted per gene per cell
4. Output: a **count matrix** (sparse matrix, typically stored as `.h5` or
   as `barcodes.tsv` + `features.tsv` + `matrix.mtx` — the "10x format")

The count matrix is the starting point for all downstream analysis.
It typically has dimensions: ~30,000 genes × 5,000–20,000 cells, but >95% of
entries are zero (sparse).

**Standard tools for everything downstream:** Seurat (R) and Scanpy (Python)

---

### 4.2 Quality Control (QC)

QC filters out poor-quality cells and technical artefacts before analysis.
Poor-quality cells can dominate clustering results and create spurious "cell types".

**The three key QC metrics per cell:**

#### 1. Number of detected genes (nFeature_RNA)
- Too low → the cell is probably dead, dying, or empty droplet
- Too high → probably a doublet (two cells merged)
- Typical range: 500–5,000 genes per cell (varies by cell type and tissue)

#### 2. Total UMI count (nCount_RNA)
- Too low → dead/empty/poor quality cell
- Too high → doublet
- Correlates with nFeature_RNA (more UMIs usually means more genes detected)

#### 3. Mitochondrial gene percentage (%MT)
- The most important QC metric for dead cell detection
- Mitochondrial genes (MT-CO1, MT-ND1, etc. in human; mt-Co1, mt-Nd1 in mouse)
  are in the cytoplasm, not the nucleus
- When a cell's outer membrane is broken (dying/dead), the cytoplasmic RNA
  leaks out — but the mitochondria (being membrane-enclosed organelles) stay
  intact
- Result: dead cells have disproportionately high mitochondrial RNA
- Threshold: typically flag cells with >10–20% MT reads (varies by tissue;
  heart tissue naturally has high mitochondrial content so thresholds are higher)

**The filtering decision:**
```
Keep a cell if:
  nFeature_RNA > 200      (not an empty droplet)
  nFeature_RNA < 5000     (not a doublet — adjust per experiment)
  percent.mt < 10–20%     (not a dead cell)
```

These thresholds are not universal — they must be **examined per dataset** by
plotting violins and scatter plots. Blindly applying thresholds from a tutorial
to a different tissue type is a common mistake.

**Ambient RNA / soup correction:**
Tools like **SoupX** estimate and remove background RNA that leaked from lysed
cells into the surrounding solution (the "soup"). If uncorrected, every cell
looks like it weakly expresses the most abundant genes from surrounding cells.

---

### 4.3 Normalisation

**The problem:** Different cells have different total UMI counts — not because
they're biologically different, but because of technical variation in capture
efficiency. Cell A might have 3,000 total UMIs and cell B might have 8,000,
even if they're the same cell type with the same biology.

If you compare raw counts, genes in cell B will always look "more expressed"
simply because more molecules were captured from it.

**Standard normalisation: library size normalisation + log transformation**

```
Step 1 — Normalise: divide each cell's counts by its total UMIs × 10,000
         (this is "counts per 10,000" or CP10K — analogous to CPM in bulk)

Step 2 — Log transform: log1p(normalised_count) = log(count + 1)
         (the +1 prevents log(0); compresses dynamic range)
```

After this, all cells are on the same scale, and the log transformation
makes the data more approximately normal (which matters for downstream statistics).

**Why not DESeq2's normalisation?** DESeq2's size factors work by comparing a
cell against the geometric mean across all cells. In scRNA-seq, most genes are 0
in most cells — geometric means collapse to 0. DESeq2 normalisation breaks on
sparse scRNA-seq data. This is why scRNA-seq uses its own normalisation approaches.

**More advanced alternatives:**
- **SCTransform** (Seurat): fits a regularised negative binomial regression per
  gene, accounting for sequencing depth as a covariate. Better than simple
  log-normalisation for highly variable genes.
- **scran pooling-based normalisation**: pools cells to estimate size factors,
  then deconvolutes — handles zeros better than geometric mean approaches.

---

### 4.4 Feature Selection and Dimensionality Reduction

#### Feature selection: highly variable genes (HVGs)

A 30,000-gene matrix is too large to analyse directly, and most genes add noise
rather than signal (housekeeping genes are the same in every cell; they don't
help distinguish cell types).

We select **highly variable genes (HVGs)** — typically the 2,000–5,000 genes
that vary the most across cells while correcting for the mean-variance relationship
(genes with higher expression naturally vary more; we want genes that vary
*more than expected* for their expression level).

HVGs capture the biology that distinguishes different cell types and states.

#### PCA (Principal Component Analysis)

Even with 2,000 HVGs, the data is still too high-dimensional to visualise.
PCA projects the data into a lower-dimensional space by finding linear combinations
of genes that explain the most variance.

In practice, the top **20–50 principal components** capture most of the biologically
meaningful variation. PC1 might correlate with cell cycle, PC2 with cell type,
PC3 with a stress response, etc.

PCA is linear — it can't capture complex non-linear relationships between cells.

#### UMAP and t-SNE: 2D visualisation

**UMAP (Uniform Manifold Approximation and Projection)** and **t-SNE** are
non-linear dimensionality reduction methods that project the data (usually from
the PCA embedding) into **2 dimensions** for visualisation.

They preserve local structure: cells that are similar in high-dimensional gene
expression space end up near each other on the 2D plot.

**Critical warnings about UMAP/t-SNE:**
- The 2D layout is useful for visualisation, but **distances between clusters
  are not meaningful** — two clusters far apart on UMAP may not be more different
  than two clusters nearby
- The shape of clusters can change dramatically with different random seeds and
  hyperparameters (`n_neighbors`, `min_dist`)
- UMAP is not deterministic — rerunning gives a different layout (unless you set
  a seed)
- **Never run statistics directly on UMAP coordinates** — use PCA space or the
  original normalised expression for any quantitative analysis

**PCA is used for clustering; UMAP is only for visualisation.**

---

### 4.5 Clustering

**What happens:** Cells are grouped into clusters based on similarity in their
gene expression profiles (using the PCA embedding, not UMAP coordinates).

**The graph-based approach (most common):**

1. **Nearest-neighbour graph:** For each cell, find its k-nearest neighbours in
   PCA space (typically k = 20). Draw edges between neighbouring cells.
   This creates a graph where connected cells are transcriptionally similar.

2. **Community detection:** Apply a graph clustering algorithm (typically
   **Leiden** or the older **Louvain**) to find communities — groups of cells
   that are more connected to each other than to the rest of the graph.

3. **Resolution parameter:** Controls the number of clusters.
   - High resolution → many small clusters (over-clustering)
   - Low resolution → few large clusters (under-clustering)
   - There is no single "correct" resolution — it must be explored and validated
     against biological knowledge

**What each cluster represents:**
Each cluster should (ideally) correspond to a distinct cell type or cell state.
A cluster of cells all expressing cardiomyocyte markers is probably cardiomyocytes.
A cluster expressing macrophage markers is probably macrophages.

But clusters are **computational constructs**, not guaranteed biological entities.
Two distinct clusters might just be the same cell type in two different cell cycle
states. Or one cluster might contain a mixture of two cell types that the algorithm
didn't separate.

---

### 4.6 Cell Type Annotation

**The hardest step.** After clustering, each cluster must be labelled with a
biological identity. This requires biological knowledge.

#### Manual annotation with marker genes

The most trusted approach: examine known marker genes for each cluster.

```
Cardiomyocytes:  TNNT2, MYH7, MYH6, ACTC1
Fibroblasts:     COL1A1, VIM, DCN, POSTN
Macrophages:     CD68, CD14, LYZ, CSF1R
Endothelial:     PECAM1 (CD31), VWF, CDH5, CLDN5
Smooth muscle:   ACTA2, TAGLN, MYH11
T cells:         CD3D, CD3E, CD8A, CD4
B cells:         CD79A, MS4A1 (CD20)
NK cells:        GNLY, NKG7, KLRD1
```

You look at a **dot plot** or **violin plot** of these marker genes across
clusters, and assign identity based on which markers are enriched.

#### Automated annotation

- **SingleR**: compares each cluster's expression profile against reference
  datasets of known cell types, assigns the closest match
- **CellTypist**: machine learning classifier trained on large atlases
- **AUCell**: scores cells for activity of gene sets (e.g., cell type signatures)

Automated tools are fast but not infallible — they work best when your tissue
has well-characterised reference data. Novel cell types or unusual states may
be mislabelled.

#### Sub-clustering

After annotating major cell types, you can sub-cluster within each type to find
sub-populations. For example, within the macrophage cluster:
- Resident cardiac macrophages (TIMD4+)
- Recruited monocyte-derived macrophages (CCR2+)
- Anti-inflammatory (M2-like) macrophages
- Pro-inflammatory macrophages

This hierarchical exploration is one of the most powerful aspects of scRNA-seq.

#### Pseudotime / trajectory analysis

If you have cells at different stages of a process (e.g., cells differentiating
from stem cells to cardiomyocytes, or macrophages activating in response to injury),
you can order cells along a **pseudotime trajectory** — a reconstructed developmental
or activation path through gene expression space.

Tools: Monocle3, Slingshot, PAGA. Useful for understanding dynamics rather than
just static cell states.

---

## 5. Key Terms Glossary

| Term | Definition |
|------|-----------|
| **scRNA-seq** | Single-cell RNA sequencing — measures gene expression in individual cells |
| **snRNA-seq** | Single-nucleus RNA sequencing — profiles nuclei instead of whole cells; better for frozen tissue |
| **Cell barcode** | A short DNA sequence unique to each cell/droplet, used to assign reads to their cell of origin |
| **UMI** | Unique Molecular Identifier — a random DNA tag on each individual mRNA molecule, used to count original molecules rather than PCR copies |
| **GEM** | Gel Bead in Emulsion — the individual oil droplet in 10x Genomics that contains one cell and one bead |
| **Count matrix** | The output of alignment: rows = genes, columns = cells, values = UMI counts per gene per cell |
| **Dropout** | When a gene that is expressed in a cell is not detected — due to low capture efficiency at low expression levels; causes sparsity |
| **Sparsity** | Most entries in the count matrix are zero (typically >90%); not all zeros mean truly unexpressed genes |
| **Doublet** | A "cell" in the data that is actually two cells co-captured in the same droplet; appears as a hybrid of two cell types |
| **Ambient RNA / soup** | RNA from lysed cells that contaminated the solution; gets captured by other cells' beads; removed with SoupX/CellBender |
| **nFeature_RNA** | Number of distinct genes detected in a cell (a QC metric) |
| **nCount_RNA** | Total UMI count in a cell (a QC metric) |
| **%MT** | Percentage of UMIs from mitochondrial genes; high %MT indicates a dying/dead cell |
| **Normalisation** | Correcting for differences in total UMI counts between cells before comparison |
| **HVG** | Highly Variable Gene — genes that vary more across cells than expected for their expression level; used for dimensionality reduction |
| **PCA** | Principal Component Analysis — linear dimensionality reduction; first step before clustering |
| **UMAP** | Uniform Manifold Approximation and Projection — non-linear 2D visualisation of cell similarity; not used for clustering |
| **t-SNE** | t-Distributed Stochastic Neighbour Embedding — older 2D visualisation method; superseded by UMAP for most purposes |
| **KNN graph** | K-Nearest Neighbour graph — connects each cell to its k most similar cells in PCA space; the basis for clustering |
| **Leiden / Louvain** | Graph community detection algorithms used for clustering cells in scRNA-seq |
| **Resolution** | Parameter controlling how many clusters the algorithm finds; higher = more clusters |
| **Marker gene** | A gene specifically or highly expressed in one cell type, used to identify that cell type |
| **Dot plot** | Visualisation showing expression level (dot size = % cells expressing) and magnitude (colour) of marker genes across clusters |
| **Pseudotime** | A computational ordering of cells along a biological trajectory (e.g., differentiation or activation) |
| **Trajectory analysis** | Reconstructing the path cells take through gene expression space during a dynamic process |
| **Cell type annotation** | The process of assigning biological identities to computational clusters |
| **SingleR** | Automated cell type annotation tool that compares clusters to reference datasets |
| **Seurat** | The dominant R package for scRNA-seq analysis (developed by the Satija Lab) |
| **Scanpy** | The dominant Python package for scRNA-seq analysis (equivalent to Seurat) |
| **Cell Ranger** | 10x Genomics' software for alignment and count matrix generation from raw reads |
| **SCTransform** | Advanced normalisation method in Seurat using regularised negative binomial regression |
| **Human Cell Atlas** | International project aiming to create a reference map of every cell type in the human body |
| **Spatial transcriptomics** | Technology measuring gene expression while preserving the physical location of cells in tissue |
| **10x Chromium** | The dominant commercial platform for droplet-based scRNA-seq (10x Genomics) |

---

## 6. How scRNA-seq Connects to Our Bulk RNA-seq Work

We spent Week 2 (Days 1–4) doing bulk RNA-seq analysis of mouse myocardial
infarction — DESeq2 differential expression, GSEA pathway analysis, and drug
repositioning. Here is how scRNA-seq would extend and deepen every part of
that work.

### 6.1 Our bulk DESeq2 found: "1605 genes upregulated in MI"

**The bulk result:** After MI, 1605 genes are upregulated compared to sham.
This is the average across the whole heart.

**What scRNA-seq would tell us:** Which cell type is responsible for each
of those 1605 genes?

Example: if *Col1a1* (collagen, marker of fibrosis) is upregulated in our
bulk data, scRNA-seq would tell us: is this upregulation happening in
fibroblasts (expected), in macrophages (unexpected but possible), or diffusely
across all cell types?

This is called **deconvolution** — breaking the bulk signal back into its
cellular components.

### 6.2 Our GSEA found: "Inflammatory and fibrotic pathways activated in MI"

**The bulk result:** Pathways like TNF signalling, complement, ECM remodelling
are upregulated.

**What scRNA-seq would tell us:** Which cell types are driving each pathway?

- TNF signalling might be driven specifically by CCR2+ macrophages
- ECM remodelling might be driven by a subpopulation of activated fibroblasts
  (called myofibroblasts)
- Complement might be driven primarily by endothelial cells

In bulk, all these signals are collapsed into one number. In scRNA-seq, each
pathway can be mapped to its cellular source.

### 6.3 Our drug repositioning found: "Estradiol and progesterone as MI candidates"

**The bulk result:** Hormonal cardioprotection as a top signal.

**What scRNA-seq would tell us:**
- Which cell type does estradiol act on to produce its cardioprotective effect?
- Is it cardiomyocytes? Macrophages (reducing inflammation)? Fibroblasts
  (reducing fibrosis)?
- You could look at estrogen receptor (*ESR1*, *ESR2*, *GPER1*) expression
  across cell types in the MI scRNA-seq dataset to predict the relevant target cell

### 6.4 Deconvolution: using scRNA-seq to reinterpret bulk data

Even without running scRNA-seq on your own samples, if a published scRNA-seq
atlas of mouse MI exists, you can **deconvolve** your bulk RNA-seq data using
that atlas as a reference.

Tools: CIBERSORT, MuSiC, BayesPrism, DWLS.

These tools estimate what fraction of your bulk sample came from each cell type,
using the scRNA-seq data as a "cell type signature dictionary". This bridges the
two technologies without needing to re-run experiments.

### 6.5 Summary comparison table

| Feature | Bulk RNA-seq (Week 2) | scRNA-seq |
|---------|----------------------|-----------|
| What is measured | Average across all cells in tissue | Each individual cell separately |
| Resolution | Tissue-level | Single-cell |
| Throughput | 1 measurement per sample | 5,000–50,000 cells per sample |
| Sensitivity | High (deep sequencing per gene) | Lower per cell (sparse) |
| Cost (2024) | ~£300–500 per sample | ~£1,500–3,000 per sample |
| Complexity of analysis | Moderate | High |
| Best for | Comparing conditions (MI vs sham) | Discovering cell types; heterogeneity |
| Missing | Cell-type resolution | Spatial location; lower depth per cell |
| Key output | List of DEGs with effect sizes | UMAP + clusters + cell type identities |

---

## 7. Applications in Pharma and Drug Discovery

scRNA-seq has transformed how pharmaceutical companies approach target discovery,
patient stratification, and drug development.

### 7.1 Target Discovery at Cellular Resolution

Before scRNA-seq, drug targets were identified based on bulk gene expression or
protein data. A gene upregulated in disease was a candidate target — but you
didn't know which cell expressed it, or whether it was the disease-driving cell
or a bystander.

**scRNA-seq enables:** "This target is upregulated specifically in CCR2+ macrophages
in the infarcted heart, and those same macrophages are the ones secreting pro-fibrotic
signals. Targeting CCR2 in macrophages should reduce fibrosis without affecting
cardiomyocytes."

This cell-type-specific target identification is now standard in cardiac, oncology,
and inflammatory disease pipelines.

### 7.2 Patient Stratification and Precision Medicine

Patients diagnosed with the same disease are often biologically heterogeneous.
Two patients with "heart failure" may have completely different cellular
compositions in their failing hearts:
- Patient A: macrophage-dominated inflammatory phenotype
- Patient B: fibroblast-dominated fibrotic phenotype

These patients likely need different treatments. scRNA-seq of biopsy material or
circulating immune cells can stratify patients based on their cellular phenotype,
not just their symptoms.

**Example in oncology:** In cancer immunotherapy, whether a tumour responds to
checkpoint inhibitors (anti-PD-1) depends on the immune cell composition of the
tumour microenvironment. scRNA-seq of tumour biopsies reveals:
- Exhausted T cells (poor responders to immunotherapy)
- Activated cytotoxic T cells (good responders)
- Immunosuppressive regulatory T cells or myeloid cells (predict resistance)

Companies like Genentech, BMS, and AstraZeneca use scRNA-seq at scale to
stratify patients in clinical trials.

### 7.3 Drug Mechanism of Action

After a drug is developed, scRNA-seq can reveal what it actually does at
cellular resolution:

- "We gave this anti-inflammatory drug. Bulk RNA-seq says inflammation is
  reduced. But which immune cell subset changed?"
- scRNA-seq before and after treatment shows exactly which populations expand,
  contract, or change state

This is used to:
- Confirm on-target effects
- Identify off-target effects early
- Understand why some patients respond and others don't

### 7.4 The Human Cell Atlas and Reference Maps

The **Human Cell Atlas (HCA)** is profiling every cell type in the human body
across development, health, and disease. By 2024, it has catalogued:
- Millions of cells from hundreds of tissues
- Cell types in foetal development
- Age-related changes in cellular composition
- Disease-specific alterations

For pharma, this creates a **reference baseline**: "What do cells in a healthy
heart look like? What changes in heart failure? Which changes are consistent
across patients (core disease signature) vs. variable (patient heterogeneity)?"

### 7.5 Spatial Transcriptomics: The Next Frontier

scRNA-seq's limitation is that it destroys spatial information — you don't know
where in the tissue each cell was sitting.

**Spatial transcriptomics** (10x Visium, Slide-seq, Stereo-seq, MERFISH, seqFISH)
measures gene expression while preserving spatial location. You can see:
- "The fibrotic scar zone has these cells"
- "The border zone between scar and healthy tissue has a unique cell state"
- "The macrophages nearest to dead cardiomyocytes have a different programme
  than macrophages in the healthy remote zone"

**Pharma applications:**
- Understanding which cells in which tissue zones respond to therapy
- Identifying spatial niches of drug-resistant cells in tumours
- Drug target validation in spatial context

### 7.6 Drug Safety and Toxicology

Off-target effects of drugs can be predicted by looking at which unexpected cell
types express the drug's target. scRNA-seq of multiple tissues can flag:
- "This anti-cancer drug targets receptor X in the tumour, but receptor X is
  also highly expressed in cardiac pacemaker cells — cardiac toxicity risk"

### 7.7 Timeline of scRNA-seq in Pharma

| Era | What pharma used scRNA-seq for |
|-----|-------------------------------|
| 2017–2019 | Proof-of-concept studies; understanding tumour microenvironments |
| 2019–2021 | Standard tool in oncology target discovery; immune cell profiling |
2021–2023 | Integrated with spatial transcriptomics; applied to non-oncology indications |
| 2023–present | Routine patient stratification in clinical trials; used in regulatory submissions |

---

## 8. Limitations and Caveats

Understanding what scRNA-seq *cannot* do is as important as knowing what it can.

| Limitation | Explanation |
|------------|-------------|
| **Sparsity / dropout** | Low capture efficiency means most gene-cell combinations show 0, even if the gene is expressed. This complicates statistical testing. |
| **Snapshot in time** | Each cell is measured once — the RNA content at the moment of lysis. You see a snapshot, not a movie. Pseudotime is a reconstruction, not a direct measurement. |
| **No protein data** | mRNA levels don't always predict protein levels (translational regulation, protein stability). CITE-seq addresses this by measuring surface proteins alongside RNA. |
| **Dissociation bias** | Some cell types (neurons, adipocytes, cardiomyocytes) are under-represented because they don't survive dissociation well. snRNA-seq mitigates this. |
| **Transcriptional noise** | At the single-cell level, gene expression is "bursty" — stochastic. Two identical cells will show different counts for the same gene by chance. Biological signal must be separated from this noise. |
| **Cost and complexity** | 3–10× more expensive than bulk RNA-seq. Requires specialised equipment (10x Chromium controller). Computational analysis requires more expertise. |
| **Ambient RNA** | Leaked RNA from lysed cells contaminates droplets. Must be corrected computationally. |
| **Doublets** | Two cells in one droplet create fake hybrid cell types. Must be detected and removed. |
| **No spatial information** | Dissociation destroys physical context. A macrophage from the scar zone and one from healthy tissue look the same unless tissue was sub-dissected first. |
| **Clustering subjectivity** | The resolution parameter and annotation step both require human judgement. Different analysts can arrive at different numbers of cell types from the same data. |

---

*Week 2, Day 5 — scRNA-seq theory notes*
*Follows from: Week 2 Days 1–4 (bulk RNA-seq, DESeq2, GSEA, drug repositioning on mouse MI)*
