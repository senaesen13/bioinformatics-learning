# Week 2 Day 5 — scRNA-seq Theory

---

## What is scRNA-seq and why does it exist

Bulk RNA-seq (everything we did in Week 2) measures the **average** gene expression
across all cells in a sample. If a heart has 60% cardiomyocytes and 20% macrophages,
you get a blended signal — you can't tell which cell type is driving a given change.

**scRNA-seq measures gene expression in individual cells.** Each cell gets its own
profile. A typical experiment captures 5,000–20,000 cells and tells you: what is
every cell doing, and are there distinct populations?

This matters because disease isn't uniform. In a heart attack:
- Fibroblasts form the scar
- Macrophages drive inflammation
- Cardiomyocytes lose contractility

Bulk RNA-seq sees all three signals merged. scRNA-seq separates them.

---

## How it works

1. **Tissue dissociation** — the sample is broken into a single-cell suspension using enzymes. Stressful for cells; some activate stress genes (FOS, JUN) as an artefact.
2. **Droplet capture (10x Chromium)** — cells flow through a microfluidic chip. Each cell is captured in an oil droplet with a gel bead.
3. **Cell barcoding** — the bead carries a unique DNA tag (barcode) that gets attached to all RNA from that cell. This is how you later know which read came from which cell.
4. **UMI tagging** — each individual mRNA molecule also gets a unique random tag (UMI) before PCR. This lets you count original molecules, not PCR copies.
5. **Sequencing** — all cDNA is sequenced together. Each read carries a cell barcode + UMI + gene sequence.
6. **Count matrix** — software (Cell Ranger) aligns reads, counts UMIs per gene per cell. Output: a matrix of genes × cells.

---

## The analysis pipeline

| Step | What it does |
|------|-------------|
| **QC** | Remove dead cells (high % mitochondrial reads), empty droplets (too few genes), and doublets (two cells in one droplet). Thresholds vary by tissue. |
| **Normalisation** | Divide each cell's counts by its total UMIs × 10,000, then log-transform. Corrects for differences in how much RNA was captured per cell. |
| **Highly variable genes** | Keep the ~2,000–5,000 genes that vary most across cells. Removes uninformative housekeeping genes and reduces noise. |
| **PCA** | Compress 2,000 genes into ~30 principal components. Captures most of the biological variation in a manageable form. |
| **Clustering** | Build a nearest-neighbour graph in PCA space, then find communities (Leiden algorithm). Each community = a candidate cell type or state. |
| **UMAP** | Project cells into 2D for visualisation only. Nearby cells = similar expression. Do not use UMAP coordinates for statistics. |
| **Cell type annotation** | Label each cluster using known marker genes (e.g. TNNT2 = cardiomyocytes, CD68 = macrophages). The hardest step — requires biology knowledge. |

---

## Key terms

| Term | Plain English |
|------|--------------|
| Cell barcode | A short DNA tag unique to each cell, labels all RNA from that cell |
| UMI | A random tag on each mRNA molecule; used to count real molecules, not PCR duplicates |
| Dropout | A gene that's expressed but not detected — common at low expression levels; makes data sparse |
| Doublet | Two cells captured in one droplet; appears as a fake hybrid cell type; must be removed |
| Ambient RNA | RNA from dead cells floating in solution; contaminates other cells' counts |
| Count matrix | Genes × cells table of UMI counts; the starting point for all analysis |
| HVG | Highly variable gene — varies more across cells than expected; used for PCA |
| Resolution | A clustering parameter; higher = more clusters. No single correct value. |
| Pseudotime | A reconstructed ordering of cells along a biological process (e.g. differentiation) |
| snRNA-seq | Profiles nuclei instead of whole cells — better for frozen tissue or fragile cells |

---

## Why it matters for pharma careers

- **Target discovery** — scRNA-seq tells you *which cell type* expresses your drug target, so you can predict on-target and off-target effects before going into the clinic.
- **Patient stratification** — two patients with the same diagnosis can have completely different cellular compositions in their tissue. scRNA-seq can predict who will respond to a drug.
- **Tumour microenvironment** — cancer immunotherapy response depends heavily on the immune cells surrounding the tumour. scRNA-seq is now standard for predicting checkpoint inhibitor response.
- **Drug mechanism** — after treatment, scRNA-seq shows which cell populations expanded, contracted, or changed state — revealing what the drug actually did at cellular resolution.
- **It is now expected** — most target discovery groups at pharma companies (Roche, AZ, Pfizer, BMS) run scRNA-seq routinely. Knowing the pipeline makes you hireable.

---

## What's coming next — Week 3: Seurat PBMC 3k

Week 3 will be the hands-on practical: running a full scRNA-seq analysis in R
using **Seurat** on the classic PBMC 3k dataset (3,000 human peripheral blood
mononuclear cells from 10x Genomics).

We will go from count matrix → QC → normalisation → PCA → UMAP → clusters →
cell type annotation, producing the labelled UMAP that is the standard output
of every scRNA-seq paper.

Everything in this theory note will become concrete code.
