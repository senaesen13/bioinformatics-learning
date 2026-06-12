# Week 1 – Day 5: SRA, FastQC & Kallisto
## Complete Study Notes

---

## PART 1: What is SRA and SRR?

### The big picture

When scientists do an RNA-seq experiment, they generate millions of sequencing reads. These reads are uploaded to a public database so other researchers can use them. That database is called the **SRA**.

```
Scientist does experiment
        ↓
Uploads raw reads to SRA (public database)
        ↓
You download them using an SRR ID
        ↓
You analyse them yourself
```

### SRA – Sequence Read Archive

- SRA = **Sequence Read Archive**
- It is a public database run by NCBI (National Center for Biotechnology Information)
- It stores raw sequencing data from experiments all over the world
- It is free to access
- Think of it like a giant library of raw sequencing files

### SRR – what is it?

- SRR = **SRA Run accession number**
- Every single sequencing run in the database has a unique ID
- It always starts with **SRR** followed by numbers (e.g. `SRR12345678`)
- This ID is how you find and download a specific dataset

### The ID hierarchy (you mostly only need SRR)

```
PRJNA123456     ← BioProject ID (the whole study / paper)
    └── SRX123456     ← SRX = one experiment (one sample)
            └── SRR123456     ← SRR = one sequencing run (the actual file)
```

> You use the **SRR ID** to download the actual data. The others are just for navigation on the website.

### How to find SRR IDs

1. Go to a paper you are interested in
2. Find the "Data Availability" section
3. They will list a BioProject ID (e.g. PRJNA12345)
4. Go to https://www.ncbi.nlm.nih.gov/sra
5. Search that ID → find individual SRR numbers for each sample

---

## PART 2: Downloading SRA data

### Step 1 – Download with prefetch

`prefetch` downloads the raw SRA file (`.sra` format) to your computer.

```bash
prefetch SRR12345678
```

This creates a folder with a `.sra` file inside it.

### Step 2 – Convert to FASTQ with fasterq-dump

The `.sra` file is not readable yet. You need to convert it to **FASTQ format**.

```bash
# For paired-end data (most RNA-seq is paired-end)
fasterq-dump --split-files SRR12345678
```

The `--split-files` flag splits paired-end reads into two files:
- `SRR12345678_1.fastq` → Read 1 (forward)
- `SRR12345678_2.fastq` → Read 2 (reverse)

### Step 3 – Compress the files (good practice)

```bash
gzip SRR12345678_1.fastq
gzip SRR12345678_2.fastq
```

Now you have `.fastq.gz` files, which are compressed and faster to work with.

### Why do we do all this?

Because raw public data comes in SRA format. You need FASTQ files to run FastQC and Kallisto. This pipeline converts public data into a format you can actually analyse.

---

## PART 3: What is a FASTQ file?

Before FastQC makes sense, you need to understand what a FASTQ file is.

Each read in a FASTQ file has exactly 4 lines:

```
@SRR12345678.1               ← Line 1: read name (always starts with @)
ATCGGTTACGATCGGATTCGATCG     ← Line 2: the actual DNA sequence
+                             ← Line 3: just a + separator (ignore this)
IIIIIIHHHHGGGGFFFFEEEEDD     ← Line 4: quality score for each base
```

The quality scores on line 4 are encoded as ASCII characters. Each character represents a **Phred quality score** for the base at that position.

---

## PART 4: FastQC

### What is FastQC?

FastQC reads your FASTQ file and produces an HTML report telling you whether your data is good and where the problems are (if any).

**Why do we run it?**
Before doing any analysis, you need to check if your data is trustworthy. Bad quality data = wrong results. FastQC is always the very first thing you run after downloading data.

### Running FastQC

```bash
# On one file
fastqc sample_1.fastq.gz

# On multiple files at once
fastqc *.fastq.gz

# Save results to a specific folder
fastqc *.fastq.gz -o fastqc_results/
```

FastQC produces two output files per sample:
- `sample_fastqc.html` → open this in your browser to see the full report
- `sample_fastqc.zip` → contains the raw numbers

---

### Understanding Phred quality scores

A Phred score tells you how confident the sequencer was when it called a base.

| Phred score | Confidence | Error chance |
|---|---|---|
| Q10 | 90% | 1 in 10 bases is wrong |
| Q20 | 99% | 1 in 100 bases is wrong |
| Q30 | 99.9% | 1 in 1000 bases is wrong |
| Q40 | 99.99% | 1 in 10,000 bases is wrong |

> **Q30 is the gold standard.** You want most of your bases to be Q30 or above.

FastQC uses a traffic light system for each module:
- ✅ PASS – looks good
- ⚠️ WARN – worth checking but may be fine
- ❌ FAIL – potential problem

---

### FastQC modules explained

#### 1. Per base sequence quality
**What it shows:** A box plot of Phred scores at every position along the read (position 1, 2, 3... up to end of read).

**What you want:** All boxes in the green zone (Q > 28) across the whole read. It is completely normal for quality to drop slightly at the very end of reads.

**What FAIL means:** The end of your reads are low quality. Fix: trim those bases with Trim Galore or Trimmomatic.

---

#### 2. Per sequence quality scores
**What it shows:** Distribution of the average quality score across all your reads.

**What you want:** A single sharp peak at a high Q score (Q30+). Most reads should have a high average quality.

**What FAIL means:** Many reads have low overall quality throughout.

---

#### 3. Per base sequence content
**What it shows:** The percentage of A, T, G, C at each position across all reads.

**What you want:** Four roughly parallel lines (each base ~25% on average). The first 10–15 bases often look wavy – this is normal for RNA-seq due to random hexamer priming bias and is not a problem.

**What FAIL means:** Strong systematic sequence bias. Can indicate adapter contamination or library preparation issues.

---

#### 4. Per sequence GC content
**What it shows:** The GC content distribution across all reads, overlaid on a theoretical bell curve for your organism.

**What you want:** Your data follows a smooth bell curve that roughly matches the expected GC content.

**What FAIL means:** A strange shape usually means contamination (DNA from another organism mixed in) or a library preparation problem.

---

#### 5. Per base N content
**What it shows:** How many positions have an N (unknown base – the sequencer could not decide) instead of A/T/G/C.

**What you want:** Near zero everywhere.

**What FAIL means:** The sequencer frequently could not call bases – usually a sign of quality issues with the run itself.

---

#### 6. Sequence length distribution
**What it shows:** The distribution of read lengths.

**What you want:** A single sharp peak (e.g. all reads exactly 150 bp).

**What FAIL means:** Variable read lengths – expected if trimming has been applied, unexpected otherwise.

---

#### 7. Sequence duplication levels
**What it shows:** How many sequences appear more than once.

**What you want for RNA-seq:** High duplication is NORMAL and EXPECTED. Highly expressed genes produce thousands of identical reads. FastQC will flag this as FAIL for RNA-seq – ignore it.

**What to worry about:** Unexpectedly high duplication in whole-genome or ChIP-seq data can indicate a problem.

---

#### 8. Overrepresented sequences
**What it shows:** Sequences that appear far more often than statistically expected.

**What you want:** No hits, or hits identified as known adapters.

**What FAIL means:** Adapter contamination is the most common cause. If adapters are found, trim them before alignment.

---

#### 9. Adapter content
**What it shows:** Whether Illumina adapter sequences are present at the ends of your reads.

**What you want:** Flat lines at 0% across all positions.

**What FAIL means:** Adapter sequences got included in your reads. This happens when the DNA fragment is shorter than the read length (the sequencer reads off the end of the fragment and into the adapter). Fix: trim adapters with Trim Galore before proceeding.

---

## PART 5: Kallisto

### What is Kallisto?

Kallisto quantifies gene and transcript expression from RNA-seq reads. It answers the question: **how much is each gene/transcript expressed in this sample?**

### Why Kallisto and not a traditional aligner?

Traditional RNA-seq aligners (STAR, HISAT2) work like this:
1. Take each read
2. Find the exact genomic position it maps to
3. Count reads per gene

This is accurate but slow (hours per sample) and uses a lot of memory (40–100 GB RAM).

Kallisto uses **pseudoalignment**:
1. Take each read
2. Ask: *which transcripts is this read compatible with?* (no exact position needed)
3. Estimate how many reads came from each transcript
4. Output expression values

This makes Kallisto extremely fast (minutes per sample) with low memory requirements (~4–8 GB RAM). Accuracy is comparable for most experiments.

---

### The Kallisto workflow

```
Reference transcriptome (FASTA file from Ensembl)
        ↓
[Step 1] kallisto index  →  creates index file (.idx)
        ↓
Your FASTQ reads (downloaded from SRA)
        ↓
[Step 2] kallisto quant  →  abundance.tsv (expression values per transcript)
```

---

### Step 1: Build the index

The index is a pre-processed, compressed version of the reference transcriptome. Kallisto uses it to quickly match reads against known transcripts without reading the full FASTA file every time.

```bash
kallisto index -i mouse_transcriptome.idx Mus_musculus.GRCm38.cdna.all.fa.gz
```

- `Mus_musculus.GRCm38.cdna.all.fa.gz` → the reference transcriptome (download from Ensembl)
- `-i mouse_transcriptome.idx` → name of the index file Kallisto will create

> This step only needs to be done **once** per reference genome/transcriptome. The index file is large (1–2 GB). Always add it to `.gitignore` – never push it to GitHub.

---

### Step 2: Quantify expression

This is where Kallisto compares your reads to the index and estimates expression levels.

```bash
kallisto quant \
  -i mouse_transcriptome.idx \
  -o sample1_results \
  -b 100 \
  SRR12345678_1.fastq.gz SRR12345678_2.fastq.gz
```

Flags explained:

| Flag | Meaning |
|---|---|
| `-i` | Path to the index file |
| `-o` | Folder where results will be saved |
| `-b 100` | Bootstrap 100 times – estimates uncertainty in the quantification |
| Last two arguments | Your paired FASTQ files (R1 then R2) |

---

### Step 3: Understanding the output

Kallisto creates a results folder containing:

| File | What it is |
|---|---|
| `abundance.tsv` | Main results – expression values per transcript |
| `abundance.h5` | Same data in HDF5 format (used by sleuth for DE analysis) |
| `run_info.json` | Run summary: how many reads mapped, what parameters were used |

### Reading the abundance.tsv file

```
target_id           length   eff_length   est_counts   tpm
ENSMUST00000001     1500     1380.2       245.0        18.4
ENSMUST00000002     800      680.1        12.0         1.8
ENSMUST00000003     2100     1980.5       0.0          0.0
```

| Column | What it means |
|---|---|
| target_id | Transcript ID from Ensembl |
| length | Full length of the transcript in base pairs |
| eff_length | Effective length – adjusted for average fragment size |
| est_counts | Estimated number of reads that came from this transcript |
| tpm | TPM – the normalized expression value you actually use |

---

### What is TPM?

TPM = **Transcripts Per Million**

You cannot directly compare raw read counts between genes or samples because:
- Longer genes naturally get more reads (more sequence to catch reads)
- Samples sequenced more deeply get more reads everywhere

TPM corrects for both problems:
1. Divide counts by transcript length → removes gene length bias
2. Scale so everything sums to 1,000,000 → removes sequencing depth bias

The sum of all TPM values in one sample always equals exactly 1,000,000.

> TPM is what you use when comparing expression levels across genes or across samples.

---

### Checking the run_info.json

Always check this after running Kallisto:

```bash
cat sample1_results/run_info.json
```

Key things to look at:
- `n_processed` – how many reads were processed
- `p_pseudoaligned` – percentage of reads that mapped (aim for >60%, ideally >80%)

If p_pseudoaligned is very low (e.g. 10%), something went wrong – wrong reference transcriptome, wrong organism, or quality issues.

---

## PART 6: The full Day 5 pipeline from start to finish

```
STEP 1 – Find dataset
  → Go to NCBI SRA, search for a BioProject ID from a paper
  → Note down the SRR IDs for each sample

STEP 2 – Download
  prefetch SRR12345678

STEP 3 – Convert to FASTQ
  fasterq-dump --split-files SRR12345678

STEP 4 – Compress
  gzip SRR12345678_1.fastq SRR12345678_2.fastq

STEP 5 – Quality check
  fastqc SRR12345678_1.fastq.gz SRR12345678_2.fastq.gz -o qc_results/

STEP 6 – Read FastQC report
  → Open the .html file in your browser
  → Check: per base quality, adapter content, GC content

STEP 7 – Trim if needed (if adapters found or quality is poor)
  trim_galore --paired SRR12345678_1.fastq.gz SRR12345678_2.fastq.gz

STEP 8 – Build Kallisto index (only once per reference)
  kallisto index -i transcriptome.idx transcriptome.fasta

STEP 9 – Quantify expression
  kallisto quant -i transcriptome.idx -o sample1_out -b 100 \
  SRR12345678_1.fastq.gz SRR12345678_2.fastq.gz

STEP 10 – Check results
  cat sample1_out/run_info.json
  head sample1_out/abundance.tsv
```

---

## Quick reference glossary

| Term | What it means |
|---|---|
| SRA | Sequence Read Archive – public database of all sequencing data |
| SRR | Unique ID for one sequencing run in the SRA database |
| prefetch | Tool to download a .sra file from NCBI |
| fasterq-dump | Tool to convert .sra to .fastq files |
| FASTQ | File format storing sequencing reads + quality scores |
| Phred score | Confidence score for each base call (Q30 = 99.9% accurate) |
| FastQC | Quality control tool – checks if your reads are good quality |
| Pseudoalignment | Kallisto's method – asks which transcript fits, not where exactly |
| Index (.idx) | Pre-processed reference file that makes Kallisto fast |
| TPM | Normalized expression value – corrects for length and depth |
| est_counts | Raw estimated read counts per transcript |
| bootstrap (-b) | Repeated resampling to estimate uncertainty in quantification |
| p_pseudoaligned | % of reads Kallisto successfully mapped (aim for >60%) |

---

*Week 1 – Day 5 | Bioinformatics Learning*
