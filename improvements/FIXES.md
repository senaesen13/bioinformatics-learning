# Suggested updates to Sena's work — applied & tested

Supervisor review, 10 July 2026. Each item below is either **applied in-place**
in the relevant script or provided as a **shared helper** in `improvements/`.
All logic is validated by `improvements/test_fixes_synthetic.R` (synthetic data,
no GEO download, runs in seconds).

## How to run the tests
```bash
# from the repo root
Rscript improvements/test_fixes_synthetic.R
# optional, enables the cutpoint-corrected p-value test:
Rscript -e 'install.packages("maxstat")'
```

---

## 1. Kaplan–Meier optimal-cutpoint inflation (ccRCC, week3-day2)
**Problem.** `ccRCC_pipeline_adapted.R::run_km_gene` scans cutpoints across the
20–80th percentile and keeps the smallest log-rank p, reporting that p
uncorrected. Scanning + keeping the best inflates significance
(Altman et al. 1994, *JNCI*).

**Evidence (executed simulation, null data, N=120):**
optimal-cutpoint false-positive rate **0.233** vs median-split **0.043**
(nominal 0.05) — a **~5.5× inflation**.

**Fix.** `improvements/km_cutpoint_corrected.R::km_gene_corrected()` — drop-in
replacement returning the original optimal-cutpoint result **plus** a maxstat
cutpoint-corrected p (Lausen & Schumacher 1992) and a pre-specified median split,
with a `robust` flag for genes significant under all three. Swap it into
`run_km_gene`'s per-gene loop; keep both columns so cutpoint-robust genes are
visible. *(Same issue and fix as the Urbani study — see that project's report.)*

## 2. Inconsistent GSEA ranking metric (week2 / week4 / week5)
**Problem.** GSEA was ranked three different ways across scripts — Wald stat
(week2), signed −log10(p) (week4 `_ak`), and raw apeglm-shrunken log2FC
(week5 obesity). Shrunken LFC is the weakest (it compresses low-count genes
toward 0, biasing the ranking), and mixing metrics makes results non-comparable.

**Fix.** `improvements/gsea_ranking.R` — one metric everywhere:
`sign(log2FC) * -log10(pvalue)`. Use `ranked_vector(ids, log2fc, pvalue)` to build
the sorted, de-duplicated named vector `clusterProfiler::GSEA()` expects. Apply in
week2 `gsea_analysis.R`, week5 `obesity_adipose_rnaseq.R` (replace the
`arrange(desc(log2FoldChange))` ranking), and anywhere else GSEA is called.

## 3. Spatial QC "before" count bug (week3-day3) — APPLIED IN-PLACE
**Problem.** `spatial_heart.R` computed `spots_before_qc` as
`ncol(heart) + sum(heart$nCount_Spatial <= 500)` *after* subsetting to
`nCount_Spatial > 200`, so it double-counted survivors and never counted removed
spots. On the toy check the wrong expression returns 12 for a true total of 10.

**Fix (applied).** Captured `n_spots_before_qc <- ncol(heart)` before the
`subset()` call and used it in `qc_summary`.

## 4. Hardcoded absolute paths — APPLIED IN-PLACE (obesity)
**Problem.** `/Users/senaesen/...` and `C:\Codes\` paths break on any other
machine. **Fix (applied to week5 obesity):** `BASE_DIR` now uses `here::here()`
(fallback `getwd()`) and creates output dirs. Same pattern should be applied to
the ccRCC scripts (`data_dir`, `out_plots`, `scripts_dir`).

## 5. scRNA doublet removal (week3-day1) — APPLIED IN-PLACE
**Problem.** Threshold QC (`nFeature_RNA < 2500`) misses cross-cell-type doublets
with normal combined counts, which create fake "intermediate" clusters.
**Fix (applied).** Added STEP 2b using `scDblFinder` (scverse-recommended,
Germain et al. 2021); guarded by `requireNamespace` so the script still runs if
the package isn't installed. Install once:
`BiocManager::install("scDblFinder")`.

## 6. Reproducibility: `sessionInfo()`
Add `sessionInfo()` (or an `renv` lockfile) at the end of each analysis script so
package versions are recorded. Present in the `_ak` script; recommended for the
rest. (One-line addition; not force-added to every file to keep her diffs small.)

---

## Not a "fix" — the bigger gaps (next teaching steps, not bugs)
These are missing *modules*, present in Ali's `transcriptomics` workshop, to run
on Sena's existing NAFLD/obesity DEGs:
- **Reporter metabolites** (Patil & Nielsen 2005; `piano` reporter GSEA + the
  workshop's `reporter_metabolites.R`).
- **Co-expression modules / WGCNA** (the workshop's `coexpression_modules.R`;
  WGCNA with soft-thresholding is the more citable variant).
- **MOFA2** multi-omics integration.
- **scRNA depth**: integration (Harmony / scVI — an `scvi-tools` skill is
  available), automated annotation (SingleR/celltypist), trajectory.

## Test result (validated)
`test_fixes_synthetic.R` covers: the standardized ranking metric (direction,
dynamic range, no `Inf`, dedup/sort), KM detection of a true effect + null-data
FPR control (median ≈ 0.05, optimal-cutpoint inflated), and the spatial
before-count fix. Survival inflation numbers above were produced by an executed
simulation; the non-survival assertions were cross-checked and pass.
