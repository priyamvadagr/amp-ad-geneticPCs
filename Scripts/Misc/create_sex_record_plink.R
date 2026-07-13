################################################################################
# Create a recorded_sex file required to check against the estimated sex 
# needed to check sex using plink 
# Also check the mismatch between estimated and recorded sex
################################################################################

# add specimenID to individual df 
combined_ind_df <- read.table('~/data/processed/metadata/Combined/DivCo_AMPAD_1.0_individual_metadata.csv',
                             sep = ',', header = T)
combined_biospec_df <- read.table('~/data/processed/metadata/Combined/DivCo_AMPAD_1.0_biospecimen_metadata.csv',
                                   sep = ',', header = T)

combined_ind_df <- combined_ind_df |>
  mutate(sample.id = combined_biospec_df$specimenID[
    match(individualID, combined_biospec_df$individualID)
  ])

recorded_sex <- combined_ind_df %>%
  transmute(
    IID = sample.id,
    sex_code = case_when(
      sex == "male"   ~ 1L,
     sex == "female" ~ 2L,
      TRUE ~ 0L      # unknown/missing
    )
  ) %>%
  filter(!is.na(IID))

# PLINK format: FID IID SEX  (with --double-id, FID == IID)
sex_file <- recorded_sex %>%
  transmute(FID = IID, IID = IID, SEX = sex_code)

write.table(sex_file, "~/data/processed/recorded_sex.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)




