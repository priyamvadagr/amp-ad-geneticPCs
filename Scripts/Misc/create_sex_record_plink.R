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


sc2 <- read.table("~/data/processed/sexcheck_filt_vs_recorded.sexcheck", header = TRUE)
table(recorded = sc2$PEDSEX, inferred = sc2$SNPSEX)
sc2 %>% filter(STATUS == "PROBLEM" & PEDSEX != 0)

sc %>% filter(SNPSEX == 0) %>% pull(F) %>% summary()
hist(sc$F[sc$SNPSEX == 0], breaks = 30)

# join ancestry and see if undetermined females concentrate in certain groups
sc %>%
  filter(SNPSEX == 0) %>%
  left_join(combined_ind_df, by = c("IID" = "sample.id")) %>%
  count(race)

hist(sc$F, breaks = 100)   # the WHOLE distribution, not just SNPSEX==0


