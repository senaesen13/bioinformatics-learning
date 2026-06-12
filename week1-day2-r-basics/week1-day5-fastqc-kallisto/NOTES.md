# Week 1 — Day 5 — FastQC + Kallisto Pipeline

> Hands-on: learning the first half of the bulk RNA-seq pipeline from scratch.
> Sample used: SRR6068409 — healthy sham mouse heart, day 1 post-surgery

---

## Why This Day Matters

In Week 1 Day 4, the Kallisto output files were already pre-made for us.
Today we learned how to produce those files ourselves from raw data.

The complete pipeline looks like this:

```
RAW FASTQ FILE
    ↓
FastQC          → check data quality before doing anything
    ↓
Kallisto index  → build a searchable reference (done once per organism)
    ↓
Kallisto quant  → match reads to transcripts, count them
    ↓
abundance.h5    → output file with counts per transcript
    ↓
tximport (R)    → collapse transcripts to genes  [Week 1 Day 4]
    ↓
DESeq2 (R)      → find significantly changed genes  [Week 1 Day 4]
    ↓
Results
```

Week 1 Day 4 covered the bottom half.
Today we covered the top half — the part that produces the raw data.

---

## Key Concepts

### What is a FASTQ file?
The raw output from a sequencing machine. Contains millions of short DNA reads.
Every single read has exactly 4 lines:

```
@SRR6068409.1 read name        ← Line 1: unique ID for this read
ATGCCAGTTCGATGGCTAGCAATGGCT    ← Line 2: the actual DNA sequence
+                               ← Line 3: separator (always just +)
IIIIIIIHHHHIIIIIIHHIIIIIIHH    ← Line 4: quality score per base
```

Quality scores use Phred+33 encoding:
- I = score 40 = 99.99% confidence = excellent
- ? = score 30 = 99.9% confidence = good (Q30 is the minimum threshold)
- 5 = score 20 = 99% confidence = poor
- ! = score 0 = completely unreliable

Rule: you want >80% of bases above Q30.

### What is Paired-End Sequencing?
The sequencer reads each DNA fragment from BOTH ends simultaneously.

```
DNA fragment (300 letters):
ATGCCAGTTCGATGGCTAGCAATGGCTTAGCGATCGTAGC...

Read 1 (forward): ATGCCAGTTCGATGGC...   ← reads from left end
Read 2 (reverse): ...GCTAGCAATGGCTTAG   ← reads from right end
```

This produces TWO files per sample:
- SRR6068409_1.fastq = all forward reads
- SRR6068409_2.fastq = all reverse reads

Why paired-end?
1. More accurate — two anchors instead of one
2. Better at detecting splice junctions
3. More information about each fragment

### What is SRA?
SRA = Sequence Read Archive.
A public database run by NCBI (USA) where ALL published RNA-seq data is stored.
Every sample gets a unique accession number starting with SRR.

The heart attack experiment used in this course:
- SRR6068402 = MI mouse, day 3, replicate 1
- SRR6068403 = MI mouse, day 3, replicate 2
- SRR6068404 = sham mouse, day 3, replicate 1
- SRR6068405 = sham mouse, day 3, replicate 2
- SRR6068406 = MI mouse, day 1, replicate 1
- SRR6068407 = MI mouse, day 1, replicate 2
- SRR6068408 = sham mouse, day 1, replicate 1
- SRR6068409 = sham mouse, day 1, replicate 2  ← used today

### What is a Reference Transcriptome?
The reference transcriptome is a FASTA file containing all known transcript sequences
for an organism. It contains only exon sequences (introns already removed).

Think of it like a dictionary of all known genes:
```
>ENSMUST00000001234 Brca1
ATGCCAGTTCGATGGCTAGCAATGGCTTAGCGATCGTAGCATG...

>ENSMUST00000005678 Tp53
GCTAGCAATGGCTTAGCGATCGTAGCATGCCAGTTCGATGGC...
```

Without the reference, Kallisto has no idea what gene each read came from.
With the reference, it can match every read to a transcript and count them.

### What is Pseudoalignment?
Traditional aligners (STAR, HISAT2) align every base of every read to the genome.
This is accurate but very slow — hours per sample.

Kallisto uses pseudoalignment — instead of aligning every base, it breaks reads
into short overlapping 31-letter chunks (k-mers) and checks which transcripts
contain those same k-mers. Much faster — minutes per sample.

```
Read:       ATGCCAGTTCGATGGCTAGC
k-mers:     ATGCCAGTTCGATGGCTAG  (letters 1-31)
             TGCCAGTTCGATGGCTAGC  (letters 2-32)

If both k-mers are found in transcript Brca1 → this read came from Brca1!
Count +1 for Brca1.
```

Speed comparison:
- Traditional alignment: 4-8 hours per sample
- Kallisto pseudoalignment: 2-5 minutes per sample

---

## Tools Installed

```bash
brew install fastqc      # quality control for FASTQ files
brew install sratoolkit  # download data from SRA database
```

Verify:
```bash
fastqc --version    # should show FastQC v0.12.1
kallisto version    # should show kallisto, version 0.52.0
```

---

## Step 1 — Download FASTQ Data from SRA

```bash
# Create working folder
mkdir week1-day5-fastqc-kallisto
cd week1-day5-fastqc-kallisto

# Download first 100,000 reads from sample SRR6068409
fastq-dump --split-files -X 100000 SRR6068409
```

What each flag means:
- --split-files = split paired-end reads into two files (_1 and _2)
- -X 100000 = only download first 100,000 reads (full files are 5-20GB)

Output:
```
SRR6068409_1.fastq   26MB   forward reads
SRR6068409_2.fastq   26MB   reverse reads
```

For future reference — how to download full files:
```bash
fasterq-dump SRR6068409   # faster download tool
gzip SRR6068409.fastq     # compress to save space
```

---

## Step 2 — FastQC Quality Control

FastQC reads your FASTQ file and produces an HTML quality report.
Always run FastQC FIRST before any analysis.

```bash
# Run FastQC on both files
fastqc SRR6068409_1.fastq SRR6068409_2.fastq

# Open the HTML report in browser
open SRR6068409_1_fastqc.html
```

### FastQC traffic light system
- Green tick = PASS — no problem
- Orange exclamation = WARNING — worth checking but usually okay
- Red cross = FAIL — potential problem

### The 10 quality checks

| Check | What it tests | Our result | Notes |
|-------|--------------|-----------|-------|
| Basic Statistics | Summary info | PASS | |
| Per base sequence quality | Quality score at each read position | PASS | Most important! |
| Per sequence quality scores | Average quality per read | PASS | |
| Per base sequence content | A/T/G/C balance | FAIL | Normal for RNA-seq! |
| Per sequence GC content | GC% vs expected | WARNING | Minor, acceptable |
| Per base N content | Unreadable bases | PASS | |
| Sequence length distribution | All reads same length | PASS | |
| Sequence duplication | Duplicate reads | WARNING | Normal for RNA-seq! |
| Overrepresented sequences | Any sequence too frequent | WARNING | Usually highly expressed genes |
| Adapter content | Adapter contamination | PASS | Critical check! |

### Why RNA-seq always has some fails/warnings
- Per base sequence content FAIL = normal! First 10-12 bases are biased
  due to random hexamer priming during library preparation. Expected and acceptable.
- Sequence duplication WARNING = normal! Highly expressed genes naturally
  produce many identical reads.

### Overall verdict for our data
The two most critical checks PASSED:
1. Per base sequence quality PASS = bases are reliably called
2. Adapter content PASS = no adapter contamination

Conclusion: GOOD QUALITY DATA — safe to proceed with analysis.

---

## Step 3 — Download Mouse Reference Transcriptome

```bash
curl -O https://ftp.ensembl.org/pub/release-109/fasta/mus_musculus/cdna/Mus_musculus.GRCm39.cdna.all.fa.gz
```

What's in this file:
- Mus_musculus = mouse (Latin name)
- GRCm39 = genome build version 39 (most recent mouse genome)
- cdna = only transcribed sequences (exons, no introns)
- all = all transcripts including all splice variants
- ~140,000 transcripts from ~22,000 genes
- Size: ~48MB compressed

Reference transcriptomes for other organisms:
```bash
# Human
curl -O https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz

# Zebrafish
curl -O https://ftp.ensembl.org/pub/release-109/fasta/danio_rerio/cdna/Danio_rerio.GRCz11.cdna.all.fa.gz
```

Important: always match the reference to your organism!

---

## Step 4 — Build Kallisto Index

The index is a compressed, fast-searchable version of the reference transcriptome.
Build it ONCE per organism — you can reuse it for all samples forever.

```bash
kallisto index -i mouse_index.idx Mus_musculus.GRCm39.cdna.all.fa.gz
```

- -i mouse_index.idx = name of the index file to create
- Takes 2-3 minutes
- Output: mouse_index.idx file

Analogy: like building the index at the back of a textbook.
You do it once, then looking things up is instant.

---

## Step 5 — Run Kallisto Quantification

Now align your reads to the index and count how many came from each transcript.

```bash
kallisto quant \
  -i mouse_index.idx \
  -o SRR6068409_output \
  --bias \
  SRR6068409_1.fastq \
  SRR6068409_2.fastq
```

What each flag means:
- -i mouse_index.idx = the index we built
- -o SRR6068409_output = folder to save results in
- --bias = correct for sequence-specific bias in library prep
- last two = paired-end FASTQ files (read 1 and read 2)

Kallisto output:
```
[quant] processed 100,000 reads
[quant] 49,586 reads pseudoaligned   = 49.6% alignment rate
```

Note on alignment rate:
- 50% for a 100k subset = acceptable (small subset behaves differently)
- Full datasets typically give 70-85% alignment rate
- Below 50% on a full dataset = investigate (wrong reference? bad data?)

### Output files produced
```
SRR6068409_output/
    abundance.h5      main output file (tximport reads this in R)
    abundance.tsv     same data in readable text format
    run_info.json     summary statistics of the run
```

---

## Step 6 — Inspect the Output

```bash
# List the output files
ls SRR6068409_output/

# View the abundance table
head SRR6068409_output/abundance.tsv

# View only transcripts with counts above 0
awk '$4 > 0' SRR6068409_output/abundance.tsv | head -20
```

### abundance.tsv columns explained

| Column | Meaning |
|--------|---------|
| target_id | Transcript ID (Ensembl format) |
| length | Full transcript length in bases |
| eff_length | Effective length (adjusted for fragment size distribution) |
| est_counts | Estimated number of reads from this transcript |
| tpm | TPM normalised expression value |

TPM = Transcripts Per Million
Normalises for transcript length AND sequencing depth so samples can be compared.

---

## Step 7 — Biological Validation

How do we know the analysis is correct?
Check if the results make biological sense!

Most highly expressed transcript in our data:
- ID: ENSMUST00000082402
- TPM: 13,026 (extremely high)

Check gene name in R:
```r
load("proteinCodingGenes.Rda")
proteinCodingGenes[proteinCodingGenes$ensembl_transcript_id == "ENSMUST00000082402", ]
```

Result:
- Gene name: mt-Co1
- Full name: mitochondrially encoded cytochrome c oxidase I

Why this makes perfect sense:
- Heart muscle cells beat 70 times per minute, 24 hours a day, 7 days a week
- This requires enormous amounts of energy
- Energy is produced by mitochondria
- Cytochrome c oxidase is the key enzyme in the mitochondrial energy chain
- Finding mt-Co1 as the most expressed gene in heart tissue = correct biology

This is how scientists validate bioinformatics results:
does the biology make sense?
In our case: YES.

---

## Complete Pipeline (copy-paste ready)

```bash
# 0. Setup
mkdir week1-day5-fastqc-kallisto
cd week1-day5-fastqc-kallisto

# 1. Download 100k reads from SRA
fastq-dump --split-files -X 100000 SRR6068409

# 2. Quality check
fastqc SRR6068409_1.fastq SRR6068409_2.fastq
open SRR6068409_1_fastqc.html

# 3. Download mouse reference transcriptome
curl -O https://ftp.ensembl.org/pub/release-109/fasta/mus_musculus/cdna/Mus_musculus.GRCm39.cdna.all.fa.gz

# 4. Build Kallisto index (do once, reuse forever)
kallisto index -i mouse_index.idx Mus_musculus.GRCm39.cdna.all.fa.gz

# 5. Run Kallisto quantification
kallisto quant -i mouse_index.idx -o SRR6068409_output --bias SRR6068409_1.fastq SRR6068409_2.fastq

# 6. Check output
ls SRR6068409_output/
awk '$4 > 0' SRR6068409_output/abundance.tsv | head -20
```

---

## Key Terms Glossary

| Term | Meaning |
|------|---------|
| FASTQ | Raw sequencing file containing reads and quality scores |
| Read | A single short DNA sequence produced by the sequencer |
| Paired-end | Sequencing both ends of each DNA fragment |
| Q30 | Quality score 30 = 99.9% base accuracy (minimum acceptable) |
| SRA | Sequence Read Archive — public database of all sequencing data |
| SRR | Unique accession number for one sample in SRA |
| Reference transcriptome | All known transcript sequences for an organism |
| Kallisto index | Compressed searchable version of the reference |
| Pseudoalignment | Fast k-mer based matching instead of base-by-base alignment |
| k-mer | Short overlapping sequence chunk (31 letters by default) |
| abundance.h5 | Kallisto output — counts per transcript |
| TPM | Transcripts Per Million — normalised expression value |
| Alignment rate | Percentage of reads that matched the reference |

---

## Connection to the Rest of the Course

```
Week 1 Day 5 (today):
FASTQ → FastQC → Kallisto → abundance.h5

Week 1 Day 4 (workshop):
abundance.h5 → tximport → count matrix → DESeq2 → results

Together = the complete bulk RNA-seq pipeline.
```

Note: data files (FASTQ, reference, index) are NOT stored on GitHub.
They are too large (26MB-48MB each). Only notes and code go on GitHub.
Data lives on your local computer or on SRA.

---

*Week 1, Day 5 — FastQC and Kallisto pipeline*
*Sample: SRR6068409 — healthy mouse heart tissue, day 1*
*Validation: most expressed gene = mt-Co1 (mitochondrial energy gene — correct for heart)*
