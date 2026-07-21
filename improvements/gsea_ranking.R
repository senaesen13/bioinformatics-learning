################################################################################
# gsea_ranking.R — ONE standardized GSEA ranking metric for the whole repo
#
# WHY: across the learning scripts the GSEA input was ranked three different
# ways — Wald statistic (week2 gsea), raw shrunken log2FC (week5 obesity),
# and signed -log10(p) (week4 _ak). For results to be comparable and
# publishable, use ONE defensible metric everywhere.
#
# CHOICE: signed -log10(p-value)  =  sign(log2FoldChange) * -log10(pvalue)
#   • preserves direction (up vs down)
#   • gives dynamic range to genes with small but highly consistent changes
#   • not distorted by apeglm shrinkage (shrunken LFC compresses low-count
#     genes toward 0, which biases the ranking)
#   • standard in clinical/NAFLD RNA-seq GSEA
#
# Use `rank_metric()` to build the metric column, then `ranked_vector()` to make
# the named, sorted, de-duplicated vector clusterProfiler::GSEA() expects.
################################################################################

#' Signed -log10(p) ranking metric.
#' @param log2fc numeric vector of (preferably UNSHRUNKEN / MLE) log2 fold changes
#' @param pvalue numeric vector of raw p-values (same length)
#' @return numeric vector; genes with tiny p and consistent direction rank extreme
rank_metric <- function(log2fc, pvalue) {
  # guard p == 0 (would give Inf): floor at the smallest representable double
  pvalue <- pmax(pvalue, .Machine$double.xmin)
  sign(log2fc) * -log10(pvalue)
}

#' Build the named, sorted, de-duplicated ranked vector for clusterProfiler::GSEA.
#' @param ids gene identifiers (symbols or Entrez) matching your gene set DB
#' @param log2fc log2 fold changes
#' @param pvalue raw p-values
#' @return named numeric vector, decreasing, first occurrence kept for duplicates
ranked_vector <- function(ids, log2fc, pvalue) {
  keep <- !is.na(ids) & ids != "" & !is.na(log2fc) & !is.na(pvalue)
  ids <- ids[keep]; log2fc <- log2fc[keep]; pvalue <- pvalue[keep]
  v <- rank_metric(log2fc, pvalue)
  v <- v[order(v, decreasing = TRUE)]
  ids <- ids[order(rank_metric(log2fc, pvalue), decreasing = TRUE)]
  names(v) <- ids
  v[!duplicated(names(v))]
}
