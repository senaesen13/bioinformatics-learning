# Week 6 Day 1 — NAFLD scRNA-seq: GSE136103 Setup

Human liver and PBMC single-cell RNA-seq data from Ramachandran et al. 2019 (Nature).
This folder covers dataset download, human-only filtering, and pre-processing setup
(Seurat QC / normalisation / clustering will be in a subsequent script).

---

## Dataset Overview

**GEO Accession:** GSE136103  
**Paper:** Ramachandran et al. 2019, *Nature* — "Resolving the fibrotic niche of human
liver cirrhosis at single-cell level"  
**Organism:** Homo sapiens only (see note below on mouse samples)

### Biological Design

Cells from human liver and blood were sorted into two fractions before sequencing:
- **CD45+** (Leucocytes): immune / inflammatory cells
- **CD45−** (Other NPC): hepatic stellate cells, endothelial cells, cholangiocytes, etc.

This sorting strategy means each biological donor contributes 1–3 separate sequencing
libraries (one per sort fraction), so the library count (24) is higher than the donor
count (14).

---

## Sample Inventory

### Biological donors: 14 human subjects

| Group | Donors | Library files |
|---|---|---|
| Healthy liver | 5 (Healthy1–5) | 11 |
| Cirrhotic liver | 5 (Cirrhotic1–5) | 9 |
| PBMC (cirrhotic) | 4 (Blood1–4) | 4 |
| **Total** | **14** | **24** |

### Why 24 library files from 14 donors?

Healthy1 → 3 libraries (cd45+, cd45-A, cd45-B)  
Healthy2 → 2 libraries (cd45+, cd45-)  
Healthy3 → 3 libraries (cd45+, cd45-A, cd45-B)  
Healthy4 → 2 libraries (cd45+, cd45-)  
Healthy5 → 1 library  (cd45+ only)  
Cirrhotic1 → 3 libraries (cd45+, cd45-A, cd45-B)  
Cirrhotic2 → 2 libraries (cd45+, cd45-)  
Cirrhotic3 → 2 libraries (cd45+, cd45-)  
Cirrhotic4 → 1 library  (cd45+ only)  
Cirrhotic5 → 1 library  (cd45+ only)  
Blood1–4   → 1 library each (PBMC, no sort fraction split)

All 24 libraries must be read in and merged/integrated into a single Seurat object
for downstream analysis.

### Full GSM list

| GSM | Title | Group |
|---|---|---|
| GSM4041150 | Healthy1_Cd45+ | Healthy liver |
| GSM4041151 | Healthy1_Cd45-A | Healthy liver |
| GSM4041152 | Healthy1_Cd45-B | Healthy liver |
| GSM4041153 | Healthy2_Cd45+ | Healthy liver |
| GSM4041154 | Healthy2_Cd45- | Healthy liver |
| GSM4041155 | Healthy3_Cd45+ | Healthy liver |
| GSM4041156 | Healthy3_Cd45-A | Healthy liver |
| GSM4041157 | Healthy3_Cd45-B | Healthy liver |
| GSM4041158 | Healthy4_Cd45+ | Healthy liver |
| GSM4041159 | Healthy4_Cd45- | Healthy liver |
| GSM4041160 | Healthy5_Cd45+ | Healthy liver |
| GSM4041161 | Cirrhotic1_Cd45+ | Cirrhotic liver |
| GSM4041162 | Cirrhotic1_Cd45-A | Cirrhotic liver |
| GSM4041163 | Cirrhotic1_Cd45-B | Cirrhotic liver |
| GSM4041164 | Cirrhotic2_Cd45+ | Cirrhotic liver |
| GSM4041165 | Cirrhotic2_CD45- | Cirrhotic liver |
| GSM4041166 | Cirrhotic3_CD45+ | Cirrhotic liver |
| GSM4041167 | Cirrhotic3_Cd45- | Cirrhotic liver |
| GSM4041168 | Cirrhotic4_Cd45+ | Cirrhotic liver |
| GSM4041169 | Cirrhotic5_Cd45+ | Cirrhotic liver |
| GSM4041170 | Blood1 | PBMC (cirrhotic) |
| GSM4041171 | Blood2 | PBMC (cirrhotic) |
| GSM4041172 | Blood3 | PBMC (cirrhotic) |
| GSM4041173 | Blood4 | PBMC (cirrhotic) |

---

## Important Note: Mouse Samples

**GSE136103 contains zero Mus musculus samples.** The `organism_ch1` field for all
24 samples is "Homo sapiens". The CCl4-treated mouse macrophage data referenced in
the Ramachandran 2019 paper is submitted under a separate GEO accession (not
GSE136103). No mouse filtering was needed here, but the organism check in the script
confirms this explicitly.

---

## Cirrhotic Liver Etiology (Cause of Disease)

The cirrhotic group is heterogeneous by cause:

| GSM | Donor | Cause |
|---|---|---|
| GSM4041161–163 | Cirrhotic1 | NAFLD |
| GSM4041164–165 | Cirrhotic2 | Alcohol |
| GSM4041166–167 | Cirrhotic3 | Alcohol |
| GSM4041168 | Cirrhotic4 | NAFLD |
| GSM4041169 | Cirrhotic5 | PBC |

This matters for interpretation: if studying NAFLD-specific biology, Cirrhotic2, 3,
and 5 are non-NAFLD cirrhosis and may need to be handled carefully in differential
analyses.

Similarly, PBMC donors:
- Blood1 / Blood3 / Blood4: NAFLD cirrhosis
- Blood2: Hereditary Haemochromatosis

---

## File Format

Each GSM folder under `data/` contains three 10x Genomics Market Exchange (MEX) files:
- `*_barcodes.tsv.gz` — cell barcodes
- `*_genes.tsv.gz` — gene list (Ensembl ID + symbol)
- `*_matrix.mtx.gz` — sparse count matrix

These are loaded with `Read10X()` in Seurat (note: file names are not the standard
`barcodes.tsv.gz / features.tsv.gz / matrix.mtx.gz` so `Read10X()` needs the
directory path and the non-standard filenames must be renamed or read manually).

---

## Output Files (this step)

| File | Description |
|---|---|
| `scripts/01_download_setup.R` | Metadata download, human filtering, count data download |
| `results/human_pdata.rds` | Full pData for 24 human libraries (RDS) |
| `results/human_sample_metadata.csv` | GSM, title, organism, group label (CSV) |
| `data/GSM404XXXX/` | 24 directories, each with barcodes / genes / matrix files |

---

## Next Steps (Week 6 Day 2)

1. Read all 24 10x directories into R with `Read10X()` / `CreateSeuratObject()`
2. Add metadata (donor, disease status, CD45 fraction, disease cause) to each object
3. Merge into a single Seurat object
4. QC: nFeature_RNA, nCount_RNA, percent mitochondrial genes
5. Filter low-quality cells and doublets
6. Normalise + find variable features
7. Integrate across donors (Harmony or Seurat integration) to remove batch effects
8. Cluster + UMAP, annotate cell types
