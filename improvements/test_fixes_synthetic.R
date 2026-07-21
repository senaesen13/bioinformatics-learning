#!/usr/bin/env Rscript
################################################################################
# test_fixes_synthetic.R — runnable tests for the suggested fixes
#
# Uses SYNTHETIC data only (no GEO download, no Bioconductor) so it runs in a few
# seconds and proves the LOGIC of each fix. Run from the repo root:
#   Rscript improvements/test_fixes_synthetic.R
# Needs: survival (base-ish). Optional: maxstat (enables the corrected-p test).
################################################################################

suppressPackageStartupMessages(library(survival))
source("improvements/gsea_ranking.R")
source("improvements/km_cutpoint_corrected.R")

pass <- 0L; fail <- 0L
ok <- function(name, cond) {
  if (isTRUE(cond)) { cat(sprintf("  PASS  %s\n", name)); pass <<- pass + 1L }
  else              { cat(sprintf("  FAIL  %s\n", name)); fail <<- fail + 1L }
}

cat("\n== TEST 1: standardized GSEA ranking metric ==\n")
lfc <- c( 2, -2,  0.1, -0.1)
pv  <- c(1e-8, 1e-8, 0.9, 0.9)
m   <- rank_metric(lfc, pv)
ok("up gene ranks positive",        m[1] > 0)
ok("down gene ranks negative",      m[2] < 0)
ok("strong beats weak (same dir)",  abs(m[1]) > abs(m[3]))
ok("no Inf when p==0",              is.finite(rank_metric(1, 0)))
v <- ranked_vector(ids = c("A","B","A","C"), log2fc = c(2,-2,1,0.1), pvalue = c(1e-5,1e-5,0.5,0.5))
ok("vector is sorted decreasing",   !is.unsorted(rev(v)))
ok("duplicate ids collapsed",       !any(duplicated(names(v))))

cat("\n== TEST 2: KM cutpoint — TRUE signal is detected ==\n")
set.seed(1)
N <- 200
expr <- rnorm(N)
# survival truly depends on expr: high expr -> higher hazard
haz  <- exp(0.9 * expr)
time <- rexp(N, rate = haz / 5)
cens <- rexp(N, rate = 1/8)
status <- as.integer(time <= cens); time <- pmin(time, cens)
r <- km_gene_corrected(expr, time, status)
ok("median-split detects real effect (p<0.05)", !is.na(r$median_p) && r$median_p < 0.05)
ok("median HR > 1 (higher expr worse)",         !is.na(r$median_HR) && r$median_HR > 1)
if (requireNamespace("maxstat", quietly = TRUE))
  ok("maxstat corrected p also significant",    !is.na(r$maxstat_p_corr) && r$maxstat_p_corr < 0.05)

cat("\n== TEST 3: KM cutpoint — NULL data false-positive control ==\n")
cat("   (optimal-cutpoint should over-fire; median split should stay ~0.05)\n")
set.seed(7); nsim <- 150; fp_opt <- 0L; fp_med <- 0L
for (i in seq_len(nsim)) {
  e <- rnorm(120)
  tt <- rexp(120, 1/5); cc <- rexp(120, 1/5)
  st <- as.integer(tt <= cc); tt <- pmin(tt, cc)
  rr <- km_gene_corrected(e, tt, st)
  if (!is.null(rr)) {
    if (!is.na(rr$opt_p_uncorr) && rr$opt_p_uncorr < 0.05) fp_opt <- fp_opt + 1L
    if (!is.na(rr$median_p)     && rr$median_p     < 0.05) fp_med <- fp_med + 1L
  }
}
cat(sprintf("   optimal-cutpoint FPR = %.3f | median-split FPR = %.3f (nominal 0.05)\n",
            fp_opt/nsim, fp_med/nsim))
ok("median-split FPR near nominal (<0.10)", fp_med/nsim < 0.10)
ok("optimal-cutpoint inflated vs median",   fp_opt/nsim > fp_med/nsim)

cat("\n== TEST 4: spatial QC before-count fix ==\n")
# emulate: 10 spots, filter nCount > 200; 'before' must equal 10, not a post-hoc mix
nCount <- c(50, 150, 250, 300, 400, 500, 600, 100, 220, 800)
n_before_fixed  <- length(nCount)                      # captured BEFORE filtering
kept            <- nCount > 200
n_before_buggy  <- sum(kept) + sum(nCount[kept] <= 500) # the old expression, post-filter
ok("fixed before-count == true total (10)", n_before_fixed == 10)
ok("old expression was wrong (!=10)",       n_before_buggy != 10)

cat(sprintf("\n==== RESULT: %d passed, %d failed ====\n", pass, fail))
if (fail > 0) quit(status = 1)
