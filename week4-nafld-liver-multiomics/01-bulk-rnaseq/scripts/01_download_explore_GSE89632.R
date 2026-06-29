library(GEOquery)
library(Biobase)

# ── Download ──────────────────────────────────────────────────────────────────
gse <- getGEO("GSE89632", GSEMatrix = TRUE, getGPL = FALSE)

# getGEO can return a list when multiple platforms exist
if (is.list(gse)) {
  message("Multiple ExpressionSets returned — picking first")
  gse <- gse[[1]]
}

# ── Basic dimensions ──────────────────────────────────────────────────────────
cat("\n=== Dataset dimensions ===\n")
cat("Features (probes/genes):", nrow(exprs(gse)), "\n")
cat("Samples:               ", ncol(exprs(gse)), "\n")

# ── Sample metadata ───────────────────────────────────────────────────────────
pd <- pData(gse)
cat("\n=== Column names in sample metadata ===\n")
print(colnames(pd))

cat("\n=== First few rows of key columns ===\n")
key_cols <- grep("title|source|characteristics|disease|status|condition|tissue|grade|stage",
                 colnames(pd), ignore.case = TRUE, value = TRUE)
print(head(pd[, key_cols, drop = FALSE]))

# ── Condition breakdown ───────────────────────────────────────────────────────
cat("\n=== Sample characteristics ===\n")
char_cols <- grep("characteristics", colnames(pd), ignore.case = TRUE, value = TRUE)
for (col in char_cols) {
  cat("\n--", col, "--\n")
  print(table(pd[[col]]))
}

# ── Expression matrix preview ─────────────────────────────────────────────────
cat("\n=== Expression matrix (first 5 genes × 5 samples) ===\n")
print(exprs(gse)[1:5, 1:min(5, ncol(gse))])

cat("\n=== Expression value range ===\n")
cat("Min:", round(min(exprs(gse), na.rm = TRUE), 3), "\n")
cat("Max:", round(max(exprs(gse), na.rm = TRUE), 3), "\n")
cat("Any NAs:", anyNA(exprs(gse)), "\n")

# ── Feature metadata ──────────────────────────────────────────────────────────
cat("\n=== Feature (gene) metadata columns ===\n")
print(colnames(fData(gse)))

cat("\n=== First few feature rows ===\n")
print(head(fData(gse)))

# ── Save metadata tables ──────────────────────────────────────────────────────
out_dir <- here::here("01-bulk-rnaseq/results")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

write.csv(pd,       file.path(out_dir, "GSE89632_sample_metadata.csv"), row.names = TRUE)
write.csv(fData(gse), file.path(out_dir, "GSE89632_feature_metadata.csv"), row.names = TRUE)

cat("\nSaved sample and feature metadata to", out_dir, "\n")
cat("\nDone — ready for differential expression.\n")
