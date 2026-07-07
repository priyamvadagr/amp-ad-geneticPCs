################################################################################
# build_master_droplist.R
# Consolidate all sample exclusions into one master drop-list with reasons.
#
# Exclusion categories:
#   - contamination : sequence contamination > 2%
#   - relatedness   : one member dropped per related pair/cluster + duplicate/swap pair
#   - sex_mismatch  : genetic sex != recorded sex (probable swaps)
#
# A sample dropped for >1 reason gets all reasons listed (comma-separated).
# Produces:
#   - master_drop_list.csv : one row per dropped sample, with reason(s)
#   - a console summary by reason
################################################################################

library(dplyr)
library(tidyr)

out_dir <- "~/data/pca_results/QC"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ---- 1. CONTAMINATION drops (>2%) ------------------------------------------
amp_ad_wgs_df <- read.table('~/data/processed/metadata/DivCO/DivCo_assay_WGS_metadata.csv',
                            sep = ',', header = TRUE)

contamination_drops <- amp_ad_wgs_df %>%
  filter(!is.na(Percent.Sequence.Contamination) &
         Percent.Sequence.Contamination > 2) %>%
  transmute(sample.id = specimenID,
            reason = "contamination",
            detail = paste0("contamination=",
                            round(Percent.Sequence.Contamination, 2), "%"))
cat("Contamination drops (>2%):", nrow(contamination_drops), "\n")

## ---- 2. RELATEDNESS drops --------------------------------------------------
# Built earlier: duplicate/swap pair + one member per resolved pair/cluster.
# Pull from the relatedness drop outputs. Adjust paths/objects to your session.

# (a) duplicate/swap pair (both members removed at the start)
#     from the kinship table: relationship == 'dup/MZ'
kin <- read.table('~/data/pca_results/QC/Sample_relatedness_flag.csv',
                  sep = ',', header = TRUE)
dup_ids <- kin %>% filter(relationship == 'dup/MZ') %>%
  { unique(c(.$ID1, .$ID2)) }

duplicate_drops <- tibble(sample.id = dup_ids,
                          reason = "duplicate/swap",
                          detail = "dup/MZ-level kinship with discordant metadata")

# (b) relatedness drops (one per disjoint pair + triad removals)
#     these were written to CSVs by the resolution script; read whichever you have.
#     Expected: a column 'sample.id' for the dropped members.
relatedness_drops <- read.csv('~/data/pca_results/QC/Sample_relatedness_drop_flag.csv')

relatedness_drops <- relatedness_drops %>% filter(kinship_related_drop) 

relatedness_drops <- relatedness_drops[!(relatedness_drops$sample.id %in% dup_ids), ]

relatedness_drops <- relatedness_drops %>% mutate("reason" = drop_reason) %>%
  mutate(detail = case_when(reason == "relatedness" ~ "one member dropped per related pair/cluster",
                    TRUE ~ "dropped due to contamination also flagged in kinship analysis")) %>%
                    select(sample.id, reason, detail)

relatedness_ids <- relatedness_drops$sample.id
cat("Relatedness drops (incl. duplicate/swap):",
    length(dup_ids) + length(relatedness_ids), "\n")

## ---- 3. SEX MISMATCH drops -------------------------------------------------
sex_mismatch <- read.table('~/data/pca_results/QC/suggested_drop_sex_mismatch.txt',
                           header = TRUE)
# expects a column 'sample.id'
sex_drops <- tibble(sample.id = sex_mismatch$sample.id,
                    reason = "sex_mismatch",
                    detail = "genetic sex != recorded sex")
cat("Sex-mismatch drops:", nrow(sex_drops), "\n")

## ---- 4. COMBINE all drops, collapsing multi-reason samples -----------------
all_drops <- bind_rows(contamination_drops,
                       duplicate_drops,
                       relatedness_drop,
                       sex_drops) %>%
  distinct()

# collapse to one row per sample, listing all reasons + details
master_drop_list <- all_drops %>%
  group_by(sample.id) %>%
  summarise(
    n_reasons   = n_distinct(reason),
    drop_reason = paste(sort(unique(reason)), collapse = ", "),
    detail      = paste(unique(detail), collapse = " | "),
    .groups = "drop"
  ) %>%
  arrange(desc(n_reasons), drop_reason, sample.id)

## ---- 5. Report + write -----------------------------------------------------
cat("\n=== MASTER DROP LIST SUMMARY ===\n")
cat("Total unique samples to drop:", nrow(master_drop_list), "\n\n")

cat("By individual reason (samples may appear in >1):\n")
print(all_drops %>% count(reason, name = "n_samples"))

cat("\nSamples flagged for MULTIPLE reasons:\n")
print(master_drop_list %>% filter(n_reasons > 1) %>% select(sample.id, drop_reason))

write.csv(master_drop_list,
          file.path(out_dir, "master_sample_drop_list.csv"), row.names = FALSE)
cat("\nMaster drop list written to",
    file.path(out_dir, "master_sample_drop_list.csv"), "\n")

## ---- 6. (optional) full-cohort flag table for keep-list assembly -----------
# If you have the full sample list, mark every sample keep/drop:
all_samples <- read.table('~/data/processed/sample_ids.txt')$V1
cohort_flags <- tibble(sample.id = as.character(all_samples)) %>%
  left_join(master_drop_list %>% select(sample.id, drop_reason), by = "sample.id") %>%
  mutate(drop = !is.na(drop_reason))
cat("\nFull cohort:", nrow(cohort_flags), "samples;",
    sum(cohort_flags$drop), "to drop;",
    sum(!cohort_flags$drop), "to keep (PCA set)\n")
write.csv(cohort_flags,
          file.path(out_dir, "cohort_keep_drop_flags.csv"), row.names = FALSE)

table(master_drop_list$ drop_reason)
