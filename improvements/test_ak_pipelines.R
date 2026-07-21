#!/usr/bin/env Rscript
## ==============================================================================
## Script: test_ak_pipelines.R
## Author: Dr. Ali Kaynar (King's College London)
## Description: Master Test Suite for AK Standardized Bioinformatics Pipelines (01 to 06)
## ==============================================================================

cat("============================================================\n")
cat("RUNNING AK PIPELINES TEST SUITE (PIPELINES 01 - 06)\n")
cat("============================================================\n\n")

get_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", cmd_args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("--file=", "", file_arg[1]))))
  }
  if (exists(".rs.getScriptPath", mode = "function")) {
    sp <- .rs.getScriptPath()
    if (!is.null(sp) && nchar(sp) > 0) return(dirname(normalizePath(sp)))
  }
  return(getwd())
}

imp_dir <- get_script_dir()
cat("Improvements Directory:", imp_dir, "\n")

pipeline_files <- c(
  "01_deg_overlap_jaccard_ak.R",
  "02_coexpression_wgcna_ak.R",
  "03_gsea_enrichment_ak.R",
  "04_cross_cohort_concordance_ak.R",
  "05_drug_repositioning_lincs_ak.R",
  "06_reporter_metabolites_ak.R"
)

results_status <- list()

for (p_file in pipeline_files) {
  full_path <- file.path(imp_dir, p_file)
  cat("\n------------------------------------------------------------\n")
  cat("Testing:", p_file, "\n")
  cat("------------------------------------------------------------\n")
  
  if (!file.exists(full_path)) {
    cat("[FAIL] File not found:", full_path, "\n")
    results_status[[p_file]] <- "MISSING"
    next
  }
  
  status <- tryCatch({
    source(full_path, local = new.env())
    "PASS"
  }, error = function(e) {
    cat("[ERROR]:", conditionMessage(e), "\n")
    "FAIL"
  })
  
  results_status[[p_file]] <- status
}

cat("\n============================================================\n")
cat("SUMMARY OF TEST RESULTS\n")
cat("============================================================\n")
for (p_file in names(results_status)) {
  cat(sprintf("  %-35s : %s\n", p_file, results_status[[p_file]]))
}
cat("============================================================\n")

if (all(unlist(results_status) == "PASS")) {
  cat("\n🎉 ALL 6 AK PIPELINES EXECUTED SUCCESSFULLY WITH ZERO ERRORS!\n")
} else {
  stop("Some pipeline tests failed. Please review output trace above.")
}
