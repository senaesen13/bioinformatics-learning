# Week 1 — Day 3 — RNA Sequencing (RNA-seq)

> Theory notes covering the full RNA-seq workflow, cDNA synthesis, fragmentation, adapter anatomy, FASTQ files, and tool installation.

---

## What is RNA Sequencing?

RNA sequencing is a technique that takes a snapshot of all the genes that are **active** in a cell at a given moment. The more RNA a gene produced, the more copies you find in your data — this tells you how active (expressed) that gene was.

---

## The 7 Core Steps of RNA-seq

### Step 1 — RNA Extraction
- Break open cells using chemicals to release RNA
- RNA degrades very fast (within minutes at room temperature) → must work on ice
- Quality is measured using a **RIN score** (RNA Integrity Number) — scored 1 to 10
- Generally need RIN > 7 for reliable results
- **Analogy:** Cracking open a walnut. Cell wall = shell, RNA = nut inside.

### Step 2 — Quality Check
- Run RNA on a machine (e.g. Agilent Bioanalyzer) to measure integrity
- RIN 1 = completely degraded, RIN 10 = perfect
- **Analogy:** Checking if a cassette tape is still playable or snapped in the middle.

### Step 3 — RNA Enrichment / Selection
- Over 80% of RNA in a cell is ribosomal RNA (rRNA) — not useful for gene expression
- Two methods to remove it:
  - **Poly-A selection** → keeps only mRNA (which has a poly-A tail at the end)
  - **Ribo-depletion** → destroys rRNA, keeps everything else
- **Analogy:** Tuning out 50 loud radio stations to hear the one you actually want.

### Step 4 — cDNA Synthesis (Reverse Transcription)
- Sequencing machines cannot read RNA directly — they read **DNA only**
- Enzyme **reverse transcriptase** converts RNA → cDNA (complementary DNA)
- An **oligo-dT primer** (string of T's) sticks to the poly-A tail and acts as a starting point

```
mRNA:  5'— AAAUUGCCAGUUCGAUGGC...AAAAAAAAAA —3'  (poly-A tail)
cDNA:  3'— TTTAACGGTCAAGCTACCG...TTTTTTTTTT —5'
```

- **Analogy:** Converting a voice memo (RNA) into a text transcript (DNA) so a computer can process it.

### Step 5 — Library Preparation
- cDNA is **fragmented** into ~150 letter pieces (sequencers can only read short pieces)
- Fragments overlap so the computer can stitch them back together later
- **Adapters** are glued to both ends of every fragment (see Adapter section below)
- The full collection of fragments = **sequencing library**
- **Analogy:** Stuffing millions of letters into envelopes with barcodes so the post office knows how to handle each one.

### Step 6 — Sequencing (NGS)
- Library is loaded onto a sequencing machine (e.g. Illumina)
- Machine reads each fragment one letter at a time
- Produces tens of millions of short **reads** (50–150 letters each)
- Output is a **FASTQ file** (see FASTQ section below)

### Step 7 — Data Analysis (Bioinformatics)
- Software **aligns** each short read back to the reference genome
- Counts how many reads came from each gene:
  - Many reads = gene is **highly active / upregulated**
  - Few reads = gene is **switched off / downregulated**
- Compare two samples (e.g. healthy vs diseased) to find which genes changed
- **Analogy:** Reassembling shredded book strips, then counting how many times each word appears.

---

## DNA → mRNA → Counts — The Biological Flow

RNA-seq intercepts mRNA before it becomes a protein:

```
DNA (in nucleus)
  └─ Transcription
       └─ pre-mRNA  (raw copy — has introns + exons)
            └─ RNA splicing
                 └─ mature mRNA  (introns removed — has poly-A tail)
                      ├─ normally → Protein  (cell function)
                      └─ RNA-seq intercepts here
                           └─ Extract + reverse transcribe
                                └─ cDNA library
                                     └─ Sequence (NGS)
                                          └─ FASTQ file
                                               └─ Align to genome
                                                    └─ Count matrix
```

> **Key point:** RNA-seq measures mRNA activity — not proteins.

### What is a count matrix?

The final output — a table of genes x samples:

```
Gene       Sample_A   Sample_B   Sample_C
BRCA1        4521       203        4890
TP53          890      8821         754
GAPDH        9200      9100        9050
```

High number = gene was active. Low number = gene was quiet.
This is what you load into R and analyse with DESeq2.

---

## cDNA Synthesis — Step by Step

### Step 1 — Mature mRNA
Single-stranded RNA with a poly-A tail at the end:
```
5'— AUGCCAGUUCGAUGGCUAGC...AAAAAAAAAA —3'
```

### Step 2 — Oligo-dT primer binds
A string of T's sticks to the poly-A tail (A pairs with T):
```
mRNA:   5'—...AAAAAAAAAA—3'
                ||||||||||
primer: 3'—TTTTTTTTTT—5'
```

### Step 3 — Reverse transcriptase copies
Enzyme reads the mRNA and builds a complementary DNA strand:
```
mRNA:  5'— AUGCCAGUUCGAUGGC...AAAAAAAAAA —3'
cDNA:  3'— TACGGTCAAGCTACCG...TTTTTTTTTT —5'
```

### Step 4 — mRNA is removed
Original mRNA destroyed by RNase H. Only the cDNA strand remains.

### Step 5 — Second strand synthesis
DNA polymerase builds a second strand → double-stranded cDNA (dsDNA):
```
1st strand: 3'— TACGGTCAAGCTACCG... —5'
2nd strand: 5'— ATGCCAGTTCGATGGC... —3'
```

### Step 6 — Fragmentation
dsDNA cut into ~150bp overlapping pieces:
```
Full cDNA:   ATGCCAGTTCGATGGCTAGCAATGGCT...  (2000+ letters)

Fragments:   ATGCCAGTTCGATGGC        (fragment 1)
                 CAGTTCGATGGCTAG      (fragment 2, overlaps)
                       GATGGCTAGCAAT  (fragment 3, overlaps)
```

---

## Adapters — Anatomy

Every fragment gets a synthetic DNA tag (adapter) on both ends:

```
5'—[ P5 ]—[ Barcode ]—[ SP1 ]—— cDNA fragment ——[ SP2 ]—[ Barcode ]—[ P7 ]—3'
       LEFT ADAPTER                                        RIGHT ADAPTER
```

| Zone | Name | Job |
|------|------|-----|
| P5 / P7 | Flow cell binding | Sticks fragment to sequencer surface |
| Barcode | Sample index | Identifies which sample this read came from |
| SP1 / SP2 | Sequencing primer | Tells sequencer where to start reading |

> **Common mistake:** The barcode tells you the **sample origin**, NOT the gene.
> Gene identity is determined later by alignment.

### Parcel analogy

| Adapter part | Analogy |
|---|---|
| cDNA fragment | The letter inside |
| Full adapter | The envelope |
| P5 / P7 | Sticky grip on the conveyor belt |
| Barcode | Postcode — which sample pile it goes to |
| SP1 / SP2 | "Open here" arrow |

---

## FASTQ Files — Deep Dive

Raw output from the sequencer. Every read = exactly 4 lines:

```
@SRR1234567.1 read_001        <- Line 1: read name (starts with @)
ATGCCAGTTCGATGGCTAGCAATGGCT   <- Line 2: DNA sequence (the real data)
+                              <- Line 3: separator (always just +)
IIIIIIIHHHHIIIIIIHHIIIIIIHH   <- Line 4: quality scores per base
```

### Quality scores (Phred+33 encoding)

Each character in line 4 scores the matching base in line 2:

| Character | Score | Error rate | Meaning |
|-----------|-------|------------|---------|
| I | 40 | 1 in 10,000 | Excellent |
| H | 39 | 1 in 8,000 | Excellent |
| ? | 30 | 1 in 1,000 | Good (Q30 threshold) |
| 5 | 20 | 1 in 100 | Poor |
| ! | 0 | — | Failed |

- **Q30 or above** = reliable base call (aim for >80% of bases)
- **Q20 or below** = too many errors, may need trimming
- Scores normally drop toward the end of reads — this is expected

### Scale

- One RNA-seq sample = **20–100 million reads**
- File size = **5–20 GB** compressed (.fastq.gz)
- This is why FastQC (Week 2) checks quality before any analysis

---

## Tools Installed — Week 1

| Tool | Type | Purpose |
|------|------|---------|
| Homebrew | Mac package manager | Installs terminal tools |
| Kallisto | Terminal (command line) | FASTQ → gene counts |
| DESeq2 | R / Bioconductor | Differential expression analysis |
| tximport | R / Bioconductor | Import Kallisto output into R |
| ggplot2 | R / CRAN | Plotting |
| pheatmap | R / CRAN | Heatmaps (Week 3) |
| ggrepel | R / CRAN | Gene labels on volcano plots (Week 3) |

### The pipeline these tools build

```
FASTQ file
  └─ Kallisto (terminal)       → count table per gene
       └─ tximport (R)         → load counts into R
            └─ DESeq2 (R)      → find significant genes
                 └─ ggplot2 + ggrepel  → volcano plot
                 └─ pheatmap           → heatmap
```

### Installation commands (Mac)

```bash
# 1. Homebrew
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

# 2. Kallisto
brew install kallisto
kallisto version   # verify

# 3. DESeq2
Rscript -e 'install.packages("BiocManager", repos="https://cran.r-project.org"); BiocManager::install("DESeq2")'

# 4. Supporting R packages
Rscript -e 'install.packages(c("ggplot2","dplyr","pheatmap","ggrepel"), repos="https://cran.r-project.org")'

# 5. tximport (Bioconductor)
Rscript -e 'BiocManager::install("tximport")'

# 6. Verify everything
Rscript -e 'library(DESeq2); library(tximport); library(ggplot2); library(pheatmap); library(ggrepel); cat("All packages working!\n")'
```

---

## Key Terms

| Term | Meaning |
|---|---|
| mRNA | Messenger RNA — carries gene instructions from DNA |
| rRNA | Ribosomal RNA — structural RNA, makes up 80%+ of total RNA |
| pre-mRNA | Raw RNA copy before splicing — contains introns |
| Intron | Non-coding section of pre-mRNA — removed during splicing |
| Exon | Coding section of pre-mRNA — kept after splicing |
| cDNA | Complementary DNA — a DNA copy of an RNA molecule |
| Reverse transcriptase | Enzyme that converts RNA → DNA |
| Poly-A tail | String of A's at the end of mRNA |
| RIN score | RNA Integrity Number — quality measure (1–10) |
| Adapter | Synthetic DNA tag added to fragments for sequencing |
| Barcode / Index | Short unique sequence identifying sample origin |
| Multiplexing | Running multiple samples in one sequencing run |
| Read | A short DNA sequence produced by the sequencer |
| FASTQ | File format storing reads and quality scores |
| Q30 | Quality threshold — less than 1 error per 1000 bases |
| Alignment | Mapping reads back to reference genome positions |
| Count matrix | Table of genes x samples showing expression levels |
| Flow cell | Glass surface inside the sequencer where fragments bind |
| Library | Complete collection of adapter-ligated cDNA fragments |
| Kallisto | Tool that aligns reads to transcriptome and counts them |
| DESeq2 | R package for differential gene expression analysis |

---

## The Big Picture

```
Cells
  └─ Extract RNA → quality check (RIN > 7)
       └─ Remove rRNA (poly-A selection or ribo-depletion)
            └─ Reverse transcribe → cDNA
                 └─ Fragment + add adapters → Library
                      └─ Sequence (NGS) → FASTQ file
                           └─ Kallisto → count table
                                └─ DESeq2 in R → significant genes
                                     └─ Volcano plot / heatmap
```

**Core idea:** More RNA from a gene → more reads → gene was more active.

---

*Bioinformatics learning journey — Week 1, Day 3*
*Topics: RNA-seq theory · DNA→mRNA→counts · cDNA synthesis · FASTQ · adapters · tool installation*
