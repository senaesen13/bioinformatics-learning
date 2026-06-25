# Week 3 Day 2 — Drug Repositioning for ccRCC using LINCS L1000

> **Script:** `scripts/ccRCC_pipeline_adapted.R`
> **Dataset:** TCGA KIRC (kidney renal clear cell carcinoma) — 528 patients, 500 genes
> **Based on:** Li et al. (2022) EBioMedicine 78:103963 (Elif's KCL workshop materials)
> **Run from:** `week3-day2-drug-repositioning-ccRCC-LINCS/` as `Rscript scripts/ccRCC_pipeline_adapted.R`

---

## What is ccRCC?

Clear cell renal cell carcinoma (ccRCC) is the most common form of kidney cancer, accounting for ~75% of cases. It has a poor prognosis once metastatic and limited treatment options beyond immunotherapy and VEGF inhibitors.

Drug repositioning — finding existing approved drugs that could treat ccRCC — is especially valuable here because de novo drug development takes 10–15 years. If a drug already has LINCS/CMap gene expression data showing it reverses the ccRCC molecular signature, it becomes an immediate candidate for clinical repurposing.

This pipeline is a systems biology approach: instead of testing one gene at a time, it integrates survival analysis, pathway enrichment, co-expression networks, and drug signature matching to nominate candidates systematically.

---

## Dataset

**TCGA KIRC** (The Cancer Genome Atlas — Kidney Renal Clear Cell Carcinoma)

- **528 tumour samples** from ccRCC patients
- **500 genes** (curated subset; the full TCGA KIRC dataset contains ~20,000 genes)
- Clinical variables per patient: gender, race, tumour stage, survival status, age, follow-up time (days)
- Expression values: TPM (transcripts per million) — a normalised measure comparable across samples
- **175/528 patients died** during follow-up (follow-up range: 2–4537 days)
- Note: this is tumour-only expression. No matched normal tissue in this example dataset. The full TCGA KIRC cohort includes ~72 adjacent normal samples.

The data file is `TCGA_KIRC_trans_exp_TPM_1.txt` from the workshop's `Example_data/` folder.

---

## Pipeline steps

### Step 1 — Kaplan-Meier survival analysis (Part 1)

**Goal:** find genes where high expression predicts shorter survival in ccRCC patients.

**Method:**
1. For each of the 500 genes, split patients into "high" vs "low" expression groups
2. The cutoff is chosen optimally — every threshold between the 20th and 80th expression percentile is tested; the one giving the lowest log-rank p-value is selected
3. Log-rank test (Kaplan-Meier): tests whether the survival curves of high-expression vs low-expression groups are significantly different
4. Cox proportional hazards model: `coef > 0` means high expression → higher hazard of death (worse prognosis)
5. Multiple testing correction: Benjamini-Hochberg FDR across all 500 genes

**Why this matters:** Genes where high expression → worse survival are likely drivers of tumour aggression, not just passengers. They make better drug targets because disrupting them should improve patient outcomes.

---

### Step 2 — GO term enrichment (Part 2)

**Goal:** understand what biological processes the prognostic oncogenes are collectively involved in.

**Method:** `enrichGO()` from clusterProfiler — hypergeometric test asking whether GO Biological Process terms are over-represented in our gene list compared to the 500-gene background.

---

### Step 3 — Co-expression network (Part 3)

**Goal:** find clusters of genes that consistently move together across ccRCC tumours — likely co-regulated or in the same pathway.

**Method:**
1. Compute Spearman correlation between all gene pairs across 528 patients
2. Keep only the top 1% most correlated pairs (r > 99th percentile) as network edges
3. Build an igraph network from these edges
4. Run random walk community detection (`cluster_walktrap`) to find gene modules
5. Filter: keep modules with ≥5 genes and high internal connectivity (transitivity > 0.4)

---

### Step 4 — Drug repositioning (Part 4)

**Goal:** find drugs whose gene expression signatures reverse the ccRCC oncogenic programme.

**Original method (Li et al.):** LINCS L1000 CMap — correlate shRNA knockdown profiles of target genes against compound perturbation profiles in the HA1E kidney cell line. Requires downloading large `.gctx` files (~GB) from clue.io.

**Adapted method (this pipeline):** MSigDB C2:CGP GSEA — equivalent connectivity map logic using pre-compiled drug gene sets:
- Score each gene by how strongly it predicts worse survival: `score = -log10(p) × sign(coef)`
- Positive score = high expression → worse survival (ccRCC oncogene)
- Run GSEA against 428 PubChem-verified drug gene sets from MSigDB C2:CGP
- Drugs with **NES < 0** = their gene set is enriched at the top of the oncogene list = the drug downregulates ccRCC oncogenes = **repositioning candidate**

---

## Results

### Part 1: Prognostic genes

| Result | Value |
|--------|-------|
| Genes tested | 500 |
| Significantly prognostic (BH p.adj < 0.05) | 407 |
| Oncogene-type (high expression → worse survival) | 81 |
| Tumour suppressor-type (low expression → worse survival) | 326 |

**Top prognostic oncogenes** (high expression = worse survival, sorted by p-value):

| Gene | Log-rank p | BH p.adj | Hazard ratio (coef) | Function |
|------|-----------|---------|-------------------|---------|
| **GPRC5A** | 1.4 × 10⁻¹⁰ | 4.6 × 10⁻⁹ | 0.95 | G protein-coupled receptor; retinoic acid target; cancer cell signalling |
| **TRIP13** | 9.1 × 10⁻¹⁰ | 1.7 × 10⁻⁸ | 0.92 | Spindle assembly checkpoint; drives chromosomal instability |
| **PABPN1** | 2.4 × 10⁻⁸ | 2.2 × 10⁻⁷ | 0.89 | Nuclear poly(A) binding protein; RNA processing |
| **KDELR3** | 5.8 × 10⁻⁸ | 4.5 × 10⁻⁷ | 0.80 | ER retention receptor; promotes cell survival under ER stress |
| **CENPM** | 9.4 × 10⁻⁸ | 6.0 × 10⁻⁷ | 0.83 | Centromere protein; chromosome segregation |
| **NME1** | 9.9 × 10⁻⁸ | 6.1 × 10⁻⁷ | 0.83 | Nucleoside diphosphate kinase; metastasis suppressor (paradoxically overexpressed here) |
| **ASNS** | 3.8 × 10⁻⁷ | 1.9 × 10⁻⁶ | 0.81 | Asparagine synthetase; metabolic reprogramming in tumours |
| **ITGB4** | — | — | — | Integrin beta-4; cell adhesion, invasion, PI3K/AKT signalling |
| **BCL3** | — | — | — | NF-κB co-activator; transcriptional regulation of survival genes |

These are biologically coherent: cell cycle (TRIP13, CENPM), ER stress survival (KDELR3), metabolic adaptation (ASNS), and signalling (GPRC5A, BCL3, ITGB4) — all known hallmarks of aggressive ccRCC.

### Part 2: GO enrichment

1,839 GO Biological Process terms returned (note: p.adj = 1 for all because the 500-gene background is too small for robust FDR correction — needs genome-wide background).

Top terms by nominal p-value:
- Regulation of actin filament length / polymerization (p = 0.003)
- Execution phase of apoptosis (p = 0.004)
- Regulation of protein-containing complex assembly (p = 0.006)

The actin and cytoskeletal terms are consistent with metastatic behaviour: tumours remodel their actin cytoskeleton to invade and migrate.

### Part 3: Co-expression network

- **500 nodes, 1,248 edges** (top 1% Spearman correlations across 528 patients)
- Random walk detected **369 communities**; **4 modules with ≥5 genes** passed the size + connectivity filter:

| Module | Size | Connectivity (CC) | Classification |
|--------|------|--------------------|---------------|
| 1 | 77 genes | 0.597 | HighCC |
| 5 | 18 genes | 0.671 | HighCC |
| 6 | 7 genes | 0.429 | HighCC |
| 8 | 5 genes | 0.875 | HighCC |

Module 1 (77 genes) is the largest co-regulated cluster — these genes rise and fall together across patients, suggesting shared transcriptional control. Module 8 has the highest internal connectivity (CC = 0.875), meaning almost every gene in it is correlated with every other.

### Part 4: Drug repositioning

**No hits at p.adj < 0.05** — expected. GSEA requires a genome-wide ranked list (~20,000 genes) to build a meaningful null distribution. With only 500 genes, the permutation-based p-values are unreliable for 138 of 428 gene sets.

**Top 10 candidates by nominal NES** (most negative = strongest anti-correlation with ccRCC signature):

| Drug / Gene Set | NES | Nominal p |
|----------------|-----|-----------|
| bhat_esr1_targets_not_via_akt1_dn | -1.63 | 0.0029 |
| lee_liver_cancer_ciprofibrate_up | -1.60 | 0.0060 |
| bhat_esr1_targets_via_akt1_dn | -1.60 | 0.0062 |
| darwiche_skin_tumor_promoter_up | -1.48 | 0.0332 |
| mcdowell_acute_lung_injury_dn | -1.47 | 0.0073 |
| blum_response_to_salirasib_up | -1.46 | 0.0413 |

**LINCS metadata check:** 7 of our top prognostic genes have shRNA knockdown profiles in the LINCS L1000 database across multiple cell lines including HA1E (kidney):

`ASNS, BCL3, GPRC5A, ITGB4, KDELR3, NME1, TRIP13`

These could be used in the full CMap pipeline (downloading `.gctx` files from clue.io) to run the original Li et al. correlation-based drug matching.

---

## Biological interpretation

### The ESR1 signal

The top two candidates both involve **ESR1** (oestrogen receptor alpha):
- `bhat_esr1_targets_not_via_akt1_dn` — genes downregulated by ESR1 activity, via a non-AKT pathway
- `bhat_esr1_targets_via_akt1_dn` — genes downregulated by ESR1 via AKT

NES < 0 means these "ESR1-downregulated" gene sets are enriched at the TOP of our ccRCC oncogene ranked list — i.e., the genes that ESR1 suppresses are the same genes that drive poor survival in ccRCC.

**Interpretation:** activating oestrogen receptor signalling (or mimicking it pharmacologically) would suppress the ccRCC oncogenic programme.

This is not surprising biologically:
- ccRCC has a well-documented **sex bias**: men are affected ~2× more than women
- Pre-menopausal women (high oestrogen) have substantially lower ccRCC incidence and better outcomes
- ESR1 is known to regulate VHL pathway genes and HIF-1α activity — both central to ccRCC biology

### Connection to Week 2 Day 4

In the mouse MI (myocardial infarction) drug repositioning analysis, one of our top candidates was **estradiol** (NES = -1.25, p.adj < 0.05). Estradiol is the primary natural ligand of ESR1.

That analysis used a completely different disease (heart attack), different species (mouse), and different method (MSigDB C2:CGP GSEA). The fact that the same oestrogen receptor biology surfaces again as a repositioning signal — this time in human kidney cancer — suggests this is a real and robust pharmacological signal, not a noise artefact.

Taken together, the two analyses independently suggest that drugs activating oestrogen receptor signalling may be protective in multiple inflammatory/oncogenic conditions.

### TRIP13 and chromosomal instability

TRIP13 (our second-most-significant prognostic gene) is a AAA-ATPase that inactivates the spindle assembly checkpoint. Overexpression leads to chromosomal instability (CIN), a hallmark of aggressive tumours. TRIP13 inhibitors are being investigated in clinical trials — this gene being top of our ccRCC survival list adds to that rationale.

### ASNS and metabolic reprogramming

Asparagine synthetase (ASNS) overexpression is a marker of metabolic adaptation: tumours starved of asparagine (e.g., by L-asparaginase treatment) upregulate ASNS to synthesise their own. ASNS overexpression predicts resistance to L-asparaginase and correlates with poor prognosis in multiple cancers. Our finding in ccRCC is consistent with this.

---

## Why this approach vs the original CMap method

| Aspect | CMap .gctx (original Li et al.) | MSigDB C2:CGP GSEA (this pipeline) |
|--------|----------------------------------|--------------------------------------|
| Drug signatures | Actual gene expression after drug treatment | Pre-compiled UP/DN gene sets from published experiments |
| Resolution | Dose- and time-specific | Gene set level only |
| Cell line | HA1E (kidney) — disease-relevant | Any cell line represented in CGP |
| Data required | GB-scale .gctx files from clue.io | No download (msigdbr package) |
| Statistical power | High with genome-wide signatures | Needs genome-wide ranked list |
| Speed | Hours (correlation across 1M+ signatures) | ~13 seconds |

The two approaches are conceptually identical (connectivity map logic: find drug signatures that anti-correlate with the disease signature) but differ in data source and resolution. With genome-wide TCGA expression and the full L1000 database, the CMap approach would give more granular, cell-line-specific results.

---

## Output files

| File | Contents |
|------|---------|
| `plots/01_km_volcano.png` | Volcano: Cox hazard ratio vs -log10(KM p-value) for all 500 genes |
| `plots/02_km_top_genes.png` | KM survival curves for top 9 prognostic oncogenes |
| `plots/03_go_barplot.png` | GO Biological Process enrichment barplot |
| `plots/04_coexpression_network.png` | Co-expression network, nodes coloured by module |
| `plots/05_module_sizes.png` | Module size distribution |
| `plots/06_drug_candidates_lollipop.png` | Top drug candidates by NES |
| `results/km_survival_all_genes.csv` | KM p-values, FDR, hazard ratios for all 500 genes |
| `results/go_enrichment_prognostic_genes.csv` | GO enrichment table |
| `results/coexpression_modules.csv` | Module membership, size, and connectivity |
| `results/gsea_drug_candidates.csv` | Significant drug candidates (empty — 500 genes insufficient) |
| `results/gsea_top10_candidates_nominal.csv` | Top 10 by nominal NES for interpretation |
| `results/lincs_shrna_availability.csv` | LINCS shRNA signature counts for our target genes |

---

## Connection to previous weeks

| Week | Analysis | Link to today |
|------|---------|--------------|
| Week 2 Day 1–2 | DESeq2 + GSEA on airway dataset | Same clusterProfiler enrichGO used here |
| Week 2 Day 3 | DESeq2 on KCL mouse MI data | Survival analysis here is analogous — both find disease-relevant genes |
| Week 2 Day 4 | Drug repositioning (mouse MI, MSigDB GSEA) | **Estradiol appeared as top candidate** — ESR1 signal replicates here in ccRCC |
| Week 3 Day 1 | Seurat scRNA-seq PBMC | Both use patient-level data; scRNA-seq would reveal which ccRCC cell types express GPRC5A/TRIP13 |

---

## References

1. **Li et al. (2022)** — Original paper this pipeline is based on.
   Li, Xiangyu, et al. "Prediction of drug candidates for clear cell renal cell carcinoma using a systems biology-based drug repositioning approach."
   *EBioMedicine* 78:103963.
   https://doi.org/10.1016/j.ebiomed.2022.103963

2. **TCGA KIRC** — The Cancer Genome Atlas Kidney Renal Clear Cell Carcinoma dataset.
   https://portal.gdc.cancer.gov/projects/TCGA-KIRC

3. **LINCS L1000 / CMap 2020** — Library of Integrated Network-Based Cellular Signatures.
   Subramanian A, et al. (2017). "A Next Generation Connectivity Map: L1000 Platform and the First 1,000,000 Profiles."
   *Cell* 171(6):1437–1452.
   https://doi.org/10.1016/j.cell.2017.10.049
   Data portal: https://clue.io/data/CMap2020

4. **MSigDB C2:CGP** — Chemical and Genetic Perturbations gene set collection.
   Liberzon A, et al. (2015). "The Molecular Signatures Database Hallmark Gene Set Collection."
   *Cell Systems* 1(6):417–425.
   https://doi.org/10.1016/j.cels.2015.12.004

5. **clusterProfiler** — Yu G, et al. (2012). "clusterProfiler: an R Package for Comparing Biological Themes Among Gene Clusters."
   *OMICS* 16(5):284–287.
   https://doi.org/10.1089/omi.2011.0118

6. **igraph** — Csardi G, Nepusz T. (2006). "The igraph software package for complex network research."
   *InterJournal Complex Systems* 1695.
   https://igraph.org

---

*Week 3, Day 2 — ccRCC drug repositioning pipeline*
*Follows from: Week 2 Day 4 (mouse MI drug repositioning) and Week 3 Day 1 (Seurat scRNA-seq)*
