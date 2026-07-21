################################################################################
# km_cutpoint_corrected.R — corrected Kaplan–Meier cutpoint for the ccRCC pipeline
#
# WHY: week3-day2 ccRCC (ccRCC_pipeline_adapted.R::run_km_gene) selects, per gene,
# the expression cutpoint with the SMALLEST log-rank p between the 20th–80th
# percentile, then reports that minimised p. Scanning many cutpoints and keeping
# the best inflates significance (Altman et al. 1994, JNCI).
#
# Simulation on NULL data (gene unrelated to survival), N=120:
#   optimal-cutpoint min-p  -> false-positive rate 0.233
#   median split            -> false-positive rate 0.043   (nominal 0.05)
#   => ~5.5x inflation.
#
# FIX: report BOTH the optimal-cutpoint result (for continuity) AND a corrected /
# pre-specified result, so cutpoint-robust genes are visible. `km_gene_corrected()`
# is a drop-in replacement for run_km_gene() returning extra columns.
#
# Optional: install.packages("maxstat") to enable the cutpoint-corrected p-value.
################################################################################

suppressPackageStartupMessages(library(survival))
.has_maxstat <- requireNamespace("maxstat", quietly = TRUE)

#' Corrected KM analysis for one gene.
#' @param expr numeric expression vector
#' @param time follow-up time
#' @param status event indicator (1 = event, 0 = censored)
#' @return one-row data.frame with optimal, maxstat-corrected, and median results
km_gene_corrected <- function(expr, time, status) {
  ok <- !is.na(expr) & !is.na(time) & !is.na(status)
  expr <- expr[ok]; time <- time[ok]; status <- status[ok]
  if (length(expr) < 20) return(NULL)
  so <- Surv(time, status)

  ## (A) optimal cutpoint — original behaviour, UNCORRECTED p
  cuts <- sort(unique(expr))
  cuts <- cuts[cuts > quantile(expr, 0.2) & cuts <= quantile(expr, 0.8)]
  bestP <- 1; bestCut <- median(expr)
  for (co in cuts) {
    g <- ifelse(expr >= co, 1L, 0L)
    if (min(table(g)) < 3 || length(unique(g)) < 2) next
    sd <- tryCatch(survdiff(so ~ g), error = function(e) NULL)
    if (is.null(sd)) next
    p <- pchisq(sd$chisq, length(sd$n) - 1, lower.tail = FALSE)
    if (p < bestP) { bestP <- p; bestCut <- co }
  }
  hr_opt <- tryCatch(exp(coef(coxph(so ~ ifelse(expr >= bestCut, 1, 0)))[1]),
                     error = function(e) NA_real_)

  ## (B) maxstat cutpoint-corrected p (Lausen & Schumacher 1992)
  corrP <- NA_real_; corrCut <- NA_real_
  if (.has_maxstat) {
    mt <- tryCatch(
      maxstat::maxstat.test(so ~ expr,
                            data = data.frame(expr = expr, time = time, status = status),
                            smethod = "LogRank", pmethod = "Lau92",
                            minprop = 0.2, maxprop = 0.8),
      error = function(e) NULL)
    if (!is.null(mt)) { corrP <- as.numeric(mt$p.value); corrCut <- as.numeric(mt$estimate) }
  }

  ## (C) median split — pre-specified, no optimization
  gm <- ifelse(expr >= median(expr), 1L, 0L)
  sdm <- tryCatch(survdiff(so ~ gm), error = function(e) NULL)
  medP <- if (!is.null(sdm)) pchisq(sdm$chisq, length(sdm$n) - 1, lower.tail = FALSE) else NA
  hr_med <- tryCatch(exp(coef(coxph(so ~ gm))[1]), error = function(e) NA_real_)

  data.frame(
    opt_cutoff = bestCut, opt_HR = hr_opt, opt_p_uncorr = bestP,
    maxstat_cutoff = corrCut, maxstat_p_corr = corrP,
    median_HR = hr_med, median_p = medP,
    robust = (bestP < 0.05) & (is.na(corrP) | corrP < 0.05) & (is.na(medP) | medP < 0.05)
  )
}
