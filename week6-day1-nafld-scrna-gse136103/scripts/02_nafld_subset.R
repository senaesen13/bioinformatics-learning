## Week 6 Day 1 — GSE136103: NAFLD-specific subset definition
## Keep: 5 healthy liver donors + Cirrhotic1 (NAFLD) + Cirrhotic4 (NAFLD)
## Exclude: Cirrhotic2/3 (alcohol), Cirrhotic5 (PBC), all PBMC
## Output: nafld_subset_metadata.csv for use in the Seurat pipeline

meta <- read.csv("../results/human_sample_metadata.csv", stringsAsFactors = FALSE)

# ── Parse donor and CD45 fraction from title ─────────────────────────────────

meta$donor <- sub("_(Cd45|CD45).*", "", meta$title, ignore.case = TRUE)
meta$donor <- ifelse(meta$donor == meta$title, meta$title, meta$donor)  # PBMC: Blood1-4

meta$fraction <- dplyr::case_when(
  grepl("Cd45\\+|CD45\\+", meta$title, ignore.case = TRUE) ~ "CD45pos",
  grepl("Cd45-|CD45-",     meta$title, ignore.case = TRUE) ~ "CD45neg",
  TRUE                                                       ~ "unsorted"
)

# ── Assign etiology ───────────────────────────────────────────────────────────

meta$etiology <- dplyr::case_when(
  meta$donor %in% c("Healthy1","Healthy2","Healthy3","Healthy4","Healthy5") ~ "None (healthy)",
  meta$donor == "Cirrhotic1" ~ "NAFLD",
  meta$donor == "Cirrhotic2" ~ "Alcohol",
  meta$donor == "Cirrhotic3" ~ "Alcohol",
  meta$donor == "Cirrhotic4" ~ "NAFLD",
  meta$donor == "Cirrhotic5" ~ "PBC",
  TRUE ~ "PBMC"
)

# ── Define NAFLD subset ───────────────────────────────────────────────────────

keep_donors <- c("Healthy1","Healthy2","Healthy3","Healthy4","Healthy5",
                 "Cirrhotic1","Cirrhotic4")

nafld_subset <- meta[meta$donor %in% keep_donors, ]
excluded     <- meta[!meta$donor %in% keep_donors, ]

# ── Report ────────────────────────────────────────────────────────────────────

message("\n========================================")
message("NAFLD SUBSET — INCLUDED (", nrow(nafld_subset), " libraries, ",
        length(unique(nafld_subset$donor)), " donors)")
message("========================================\n")
print(nafld_subset[, c("geo_accession","title","donor","fraction",
                        "group","etiology")])

message("\n========================================")
message("EXCLUDED (", nrow(excluded), " libraries)")
message("========================================\n")
print(excluded[, c("geo_accession","title","donor","group","etiology")])

message("\n--- Summary by group ---")
print(table(nafld_subset$group))

message("\n--- Summary by donor ---")
print(table(nafld_subset$donor))

# ── Check data files exist for every kept GSM ─────────────────────────────────

message("\n--- Data directory check ---")
for (gsm in nafld_subset$geo_accession) {
  d <- file.path("../data", gsm)
  files <- list.files(d)
  if (length(files) == 3) {
    message("OK  ", gsm, " (", 3, " files)")
  } else {
    message("MISSING  ", gsm, " — found ", length(files), " file(s) in ", d)
  }
}

# ── Save ──────────────────────────────────────────────────────────────────────

write.csv(nafld_subset, "../results/nafld_subset_metadata.csv", row.names = FALSE)
message("\nSaved: results/nafld_subset_metadata.csv\n")

message("========================================")
message("FINAL CONFIRMATION")
message("========================================")
message("Total GSE136103 human libraries : 24")
message("NAFLD subset (included)         : ", nrow(nafld_subset))
message("  Healthy liver libraries       : ",
        sum(nafld_subset$group == "Healthy liver"))
message("  NAFLD cirrhotic libraries     : ",
        sum(nafld_subset$group == "Cirrhotic liver"))
message("Excluded                        : ", nrow(excluded))
message("  Alcohol cirrhosis             : ",
        sum(excluded$etiology == "Alcohol"))
message("  PBC cirrhosis                 : ",
        sum(excluded$etiology == "PBC"))
message("  PBMC                          : ",
        sum(excluded$etiology == "PBMC"))
message("========================================")
message("Ready for: Seurat pipeline (script 03)")
