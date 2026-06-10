# Week 1 — Day 3 — RNA Sequencing (RNA-seq)

> Theory notes covering the full RNA-seq workflow, cDNA synthesis, fragmentation, and adapter anatomy.

---

## What is RNA Sequencing?

RNA sequencing is a technique that takes a snapshot of all the genes that are **active** in a cell at a given moment. The more RNA a gene produced, the more copies you find in your data — this tells you how active (expressed) that gene was.

---

## The 7 Core Steps of RNA-seq

### Step 1 — RNA Extraction
- Break open cells using chemicals to release RNA
- RNA degrades very fast (within minutes at room temperature) → must work on ice
- Quality is measured using a RIN score (RNA Integrity Number) — scored 1 to 10
- Generally need RIN > 7 for reliable results
- Analogy: Cracking open a walnut. Cell wall = shell, RNA = nut inside.

### Step 2 — Quality Check
- Run RNA on a machine (e.g. Agilent Bioanalyzer) to measure integrity
- RIN 1 = completely degraded, RIN 10 = perfect
- Analogy: Checking if a cassette tape is still playable or snapped in the middle.

### Step 3 — RNA Enrichment / Selection
- Over 80% of RNA in a cell is ribosomal RNA (rRNA) — not useful for gene expression
- Two methods to remove it:
  - Poly-A selection  → keeps only mRNA (which has a poly-A tail at the end)
  - Ribo-depletion   → destroys rRNA, keeps everything else
- Analogy: Tuning out 50 loud radio stations to hear the one you actually want.

### Step 4 — cDNA Synthesis (Reverse Transcription)
- Sequencing machines cannot read RNA directly — they read DNA only
- Enzyme reverse transcriptase converts RNA into cDNA (complementary DNA)
- An oligo-dT primer (string of T's) sticks to the poly-A tail and acts as a starting point

  mRNA:  5'— AAAUUGCCAGUUCGAUGGC...AAAAAAAAAA —3'  (poly-A tail)
  cDNA:  3'— TTTAACGGTCAAGCTACCG...TTTTTTTTTT —5'

- Analogy: Converting a voice memo (RNA) into a text transcript (DNA) so a computer can process it.

### Step 5 — Library Preparation
- cDNA is fragmented into ~150 letter pieces (sequencers can only read short pieces)
- Fragments overlap so the computer can stitch them back together later
- Adapters are glued to both ends of every fragment (see Adapter section below)
- The full collection of fragments = sequencing library
- Analogy: Stuffing millions of letters into envelopes with barcodes so the post office knows how to handle each one.

### Step 6 — Sequencing (NGS)
- Library is loaded onto a sequencing machine (e.g. Illumina)
- Machine reads each fragment one letter at a time
- Produces tens of millions of short reads (50–150 letters each)
- Output is a FASTQ file:

  @read_001
  TTTAACGGTCAAGCTACCGATCG
  +
  IIIIIIHHHHIIIIIHHIIIIII

  Line 1 = read name
  Line 2 = DNA sequence
  Line 3 = separator
  Line 4 = quality scores per base

- Analogy: Shredding a book into millions of tiny strips, then photographing each strip.

### Step 7 — Data Analysis (Bioinformatics)
- Software aligns each short read back to the reference genome
- Counts how many reads came from each gene:
  - Many reads = gene is highly active / upregulated
  - Few reads  = gene is switched off / downregulated
- Compare two samples (e.g. healthy vs diseased) to find which genes changed
- Analogy: Reassembling the shredded book strips, then counting how many times each word appears.

---

## Deep Dive — Adapters

An adapter is a short synthetic DNA tag glued to both ends of every cDNA fragment.
Made of three distinct zones:

  5'—[ P5 ]—[ Barcode ]—[ SP1 ]—— cDNA fragment ——[ SP2 ]—[ Barcode ]—[ P7 ]—3'
         LEFT ADAPTER                                        RIGHT ADAPTER

### Zone 1 — P5 / P7 (Flow Cell Binding)
- Complementary to tiny DNA hooks on the sequencer's glass surface (flow cell)
- Fragment physically sticks to the surface and stays in place while being read
- Without these the fragment floats away and is never sequenced

### Zone 2 — Barcode (Sample Index)
- A short unique sequence of ~6–8 letters (e.g. ATCGGT)
- Identifies which SAMPLE the fragment came from — NOT which gene
- Enables multiplexing — mixing up to 96 samples in one run to save cost

  Barcode ATCGGT → Patient A tumour sample
  Barcode TTGCCA → Patient B healthy tissue

  NOTE: The barcode tells you the sample origin, not the gene.
        Gene identity is determined later by aligning the cDNA to the reference genome.

### Zone 3 — SP1 / SP2 (Sequencing Primer Binding Sites)
- Known DNA sequences the sequencer uses as a starting point
- SP1 = start reading from the left end
- SP2 = enables paired-end sequencing (reads fragment from both ends for higher accuracy)

### Adapter Summary — The Parcel Analogy

  cDNA fragment  →  the letter inside
  Full adapter   →  the envelope
  P5 / P7        →  sticky part gripping the conveyor belt
  Barcode        →  postcode telling the sorter which pile it goes to
  SP1 / SP2      →  "open here" arrow showing where to start reading

---

## Key Terms

  mRNA                 Messenger RNA — carries gene instructions from DNA
  rRNA                 Ribosomal RNA — structural RNA, makes up 80%+ of total RNA
  cDNA                 Complementary DNA — a DNA copy of an RNA molecule
  Reverse transcriptase  Enzyme that converts RNA into DNA
  Poly-A tail          String of A's at the end of mRNA, used to select mRNA
  RIN score            RNA Integrity Number — quality measure (1 to 10)
  Adapter              Synthetic DNA tag added to fragments for sequencing
  Barcode / Index      Short unique sequence identifying sample origin
  Multiplexing         Running multiple samples in one sequencing run
  Read                 A short DNA sequence produced by the sequencer
  FASTQ                File format storing reads and quality scores
  Alignment            Mapping reads back to reference genome positions
  Flow cell            Glass surface inside the sequencer where fragments bind
  Library              Complete collection of adapter-ligated cDNA fragments

---

## The Big Picture

  Cells
    └─ Extract RNA
         └─ Quality check (RIN > 7)
              └─ Remove rRNA (poly-A selection or ribo-depletion)
                   └─ Reverse transcribe → cDNA
                        └─ Fragment + add adapters → Library
                             └─ Sequence (NGS) → FASTQ file
                                  └─ Align to genome → Count reads per gene
                                       └─ Find active / inactive genes

  Core idea: more RNA from a gene → more reads → more active that gene was.

