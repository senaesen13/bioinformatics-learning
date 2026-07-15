## Week 6 Day 1 — GSE136103 scRNA-seq: metadata download & human-only filtering
## Stop after confirming correct human sample set and downloading count data.
## Do NOT run Seurat pipeline yet.

library(GEOquery)
library(dplyr)

set.seed(42)

# ── 1. Download metadata ──────────────────────────────────────────────────────

message("\n=== Downloading GSE136103 metadata via GEOquery ===\n")
gse <- getGEO("GSE136103", GSEMatrix = TRUE, getGPL = FALSE)

# getGEO returns a list when multiple platforms; take the first element
if (is.list(gse)) gse <- gse[[1]]

pdata <- pData(gse)

message("All columns in pData:\n")
print(colnames(pdata))

message("\nFull pData (selected informative columns):\n")
info_cols <- c("title", "geo_accession", "organism_ch1",
               "source_name_ch1", "characteristics_ch1",
               "characteristics_ch1.1", "characteristics_ch1.2")
info_cols <- info_cols[info_cols %in% colnames(pdata)]
print(pdata[, info_cols])

# ── 2. Organism breakdown ─────────────────────────────────────────────────────

message("\n=== Organism breakdown ===\n")
print(table(pdata$organism_ch1))

# ── 3. Filter to human only ───────────────────────────────────────────────────

human_pdata <- pdata[pdata$organism_ch1 == "Homo sapiens", ]

message("\n=== Human-only samples (", nrow(human_pdata), " total) ===\n")
print(human_pdata[, info_cols])

# Assign group labels from source / title / characteristics
# Inspect what labels are present
message("\nUnique source_name_ch1 values in human samples:\n")
print(unique(human_pdata$source_name_ch1))

message("\nUnique characteristics_ch1 values in human samples:\n")
print(unique(human_pdata$characteristics_ch1))

# Build a tidy group column
human_pdata <- human_pdata %>%
  mutate(
    group = case_when(
      grepl("PBMC|pbmc", source_name_ch1, ignore.case = TRUE) ~ "PBMC (cirrhotic)",
      grepl("cirrhot|cirrhos|fibrosis|F4", source_name_ch1, ignore.case = TRUE) ~ "Cirrhotic liver",
      grepl("normal|healthy|non.cirrhot|non.fibrotic", source_name_ch1, ignore.case = TRUE) ~ "Healthy liver",
      grepl("PBMC|pbmc", title, ignore.case = TRUE) ~ "PBMC (cirrhotic)",
      grepl("cirrhot|cirrhos", title, ignore.case = TRUE) ~ "Cirrhotic liver",
      grepl("normal|healthy", title, ignore.case = TRUE) ~ "Healthy liver",
      TRUE ~ "Unknown"
    )
  )

message("\n=== Group breakdown of human samples ===\n")
print(table(human_pdata$group))

message("\n=== Confirmed human-only sample list ===\n")
print(human_pdata[, c("geo_accession", "title", "organism_ch1", "group")])

# ── 4. Expected: 5 healthy liver, 5 cirrhotic liver, 4 PBMC = 14 total ───────

n_healthy  <- sum(human_pdata$group == "Healthy liver")
n_cirrh    <- sum(human_pdata$group == "Cirrhotic liver")
n_pbmc     <- sum(human_pdata$group == "PBMC (cirrhotic)")
n_unknown  <- sum(human_pdata$group == "Unknown")
n_total    <- nrow(human_pdata)
n_mouse    <- nrow(pdata) - n_total

message("\n========================================")
message("FILTERING SUMMARY")
message("========================================")
message("Total samples in GSE136103 : ", nrow(pdata))
message("Mouse (excluded)           : ", n_mouse)
message("Human (retained)           : ", n_total)
message("  Healthy liver            : ", n_healthy)
message("  Cirrhotic liver          : ", n_cirrh)
message("  PBMC (cirrhotic)         : ", n_pbmc)
message("  Unknown (check labels)   : ", n_unknown)
message("========================================\n")

if (n_unknown > 0) {
  message("WARNING: ", n_unknown, " samples have unrecognised labels. ",
          "Review their source_name_ch1 / title and fix the case_when() above.\n")
  print(human_pdata[human_pdata$group == "Unknown",
                    c("geo_accession", "title", "source_name_ch1",
                      "characteristics_ch1")])
}

# Save the filtered metadata for downstream use
saveRDS(human_pdata, "../results/human_pdata.rds")
write.csv(human_pdata[, c("geo_accession", "title", "organism_ch1", "group",
                           "source_name_ch1")],
          "../results/human_sample_metadata.csv", row.names = FALSE)
message("Saved: results/human_pdata.rds and results/human_sample_metadata.csv\n")

# ── 5. Download supplementary files for human samples only ───────────────────

message("=== Downloading supplementary count files for human samples only ===\n")

human_gsm_ids <- human_pdata$geo_accession

# getGEOSuppFiles downloads to a subfolder named by GSM accession
# under destdir. We point it at data/.
destdir <- "../data"

for (gsm_id in human_gsm_ids) {
  gsm_label <- human_pdata[human_pdata$geo_accession == gsm_id, "group"]
  message("Downloading: ", gsm_id, " [", gsm_label, "]")
  tryCatch(
    getGEOSuppFiles(gsm_id, makeDirectory = TRUE, baseDir = destdir),
    error = function(e) message("  ERROR for ", gsm_id, ": ", e$message)
  )
}

message("\n=== Download complete. Files written to data/<GSM>/ ===\n")
message("Listing downloaded files:\n")
downloaded <- list.files(destdir, recursive = TRUE, full.names = FALSE)
print(downloaded)

# ── 6. Final confirmation report ─────────────────────────────────────────────

message("\n========================================")
message("FINAL CONFIRMATION")
message("========================================")
message("Human samples downloaded: ", length(human_gsm_ids))
message("  Healthy liver  : ", n_healthy)
message("  Cirrhotic liver: ", n_cirrh)
message("  PBMC (cirrh.)  : ", n_pbmc)
message("Mouse samples    : 0 (explicitly excluded)")
message("Ready for: Week 6 Day 1 — Seurat QC pipeline")
message("========================================\n")
