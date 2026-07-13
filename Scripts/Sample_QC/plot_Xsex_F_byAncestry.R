################################################################################
# plot_Xsex_F_byAncestry.R
# Plot X-chromosome inbreeding F (from plink --check-sex) colored by
# self-reported race, to show that the female F bimodality is ancestry-driven.
#
# Requires:
#   sexcheck file from plink --check-sex (cols: FID IID PEDSEX SNPSEX STATUS F)
#   combined_ind_df : metadata with sample.id + race (+ optional isHispanic)
################################################################################


library(dplyr)
library(ggplot2)

out_dir <- "~/AMP-AD_genetic_PCs/Results/Filtering/Sex_check"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
F_FEMALE_MAX <- 0.5   # female if F < 0.5
F_MALE_MIN   <- 0.8   # male   if F > 0.8


## ---- 1. Load sex-check output and join ancestry + Hispanic -----------------
sc <- read.table("~/data/processed/sexcheck_filt_vs_recorded.sexcheck", header = TRUE)


# combine recorded sex across your metadata sources, keyed by the sample IDs
# in the genotype data. Adjust object/column names to yours.
combined_ind_df <- read.table('~/data/processed/metadata/Combined/DivCo_AMPAD_1.0_individual_metadata.csv',
                             sep = ',', header = T)
combined_biospec_df <- read.table('~/data/processed/metadata/Combined/DivCo_AMPAD_1.0_biospecimen_metadata.csv',
                                   sep = ',', header = T)


unique(combined_ind_df$isHispanic)
combined_ind_df <- combined_ind_df |>
  mutate(sample.id = combined_biospec_df$specimenID[
    match(individualID, combined_biospec_df$individualID)
  ])

dat <- sc %>%
  left_join(combined_ind_df %>% select(sample.id, race, isHispanic, sex),
            by = c("IID" = "sample.id")) %>%
  mutate(
    race = ifelse(is.na(race) | race == "" | race == "missing or unknown", "Missing/Unknown", 
                  ifelse(race == "Other", "Not specified", race)),
    hisp = case_when(
      isHispanic == "True"  ~ TRUE,
      isHispanic == "False"  ~ FALSE,
      TRUE ~ NA
    ),
    # genetic sex call from F using our thresholds
    genetic_sex_plink = case_when(
      SNPSEX == 2 ~ "Female",
      SNPSEX == 1  ~ "Male",
      SNPSEX == 0  ~ "Undetermined"
    ),
    # genetic sex per PLINK threshold
    genetic_sex_new_threshold = case_when(
      F < F_FEMALE_MAX ~ "Female",
      F > F_MALE_MIN   ~ "Male",
      TRUE             ~ "Undetermined"
    ),
    recorded_sex = case_when(
      sex == "male" ~ "Male",
      sex == "female" ~ "Female"
    ),
    race_hisp = case_when(
      is.na(hisp)  ~ race,
      hisp == TRUE  ~ paste0(race, " (Hispanic)"),
      hisp == FALSE ~ paste0(race, " (non-Hispanic)")
    ),
    # mismatch = recorded sex known AND disagrees with genetic call
    mismatch_new = !is.na(recorded_sex) &
                   recorded_sex != genetic_sex_new_threshold,

    mismatch_plink = !is.na(recorded_sex) &
                     recorded_sex != genetic_sex_plink
  )

cat("Mismatches (recorded != genetic plink):", sum(dat$mismatch_plink, na.rm = TRUE), "\n")
cat("Mismatches (recorded != genetic new threshold):", sum(dat$mismatch_new, na.rm = TRUE), "\n")


race_colors <-  c("White" = "#CC79A7",  # orange
               "Asian" = "#0039A6",  # dark blue
               "Black or African American" = "#009E73",  # green
               "Not specified" = "#85144b",  # pink/magenta
               "American Indian or Alaska Native" = "#D55E00")  # vermillion/red-orange  


race_colors["Missing/Unknown"] <- "#374057"   # mute the unknowns

sex_colors <- c(
  Male   = "#5566AA",  # muted indigo
  Female = "#CC8800"   # ochre/gold
)



## ---- 2. Plots with standard PLINK thersholds -----------------------------
pos <- position_jitter(width = 0.2, seed = 42)

# swarm plot faceted by race
p_strip <- ggplot(dat, aes(x = recorded_sex, y = F, colour = recorded_sex)) +
  facet_wrap(~race) +
  geom_jitter(width = 0.2, size = 1.6, alpha = 0.7) +
  scale_color_manual(values = sex_colors) +
  geom_hline(yintercept = 0.2, linetype = "dotted", colour = "grey50", , linewidth = 0.7) +
  geom_hline(yintercept = 0.8, linetype = "dashed", colour = "red", linewidth = 0.7) +
  annotate("text", x = 0.7, y = 0.25, label = "Female cutoff (F < 0.2)",
           colour = "black", hjust = 0, size = 3.5, fontface = "italic") +
  annotate("text", x = 0.7, y = 0.85, label = "Male cutoff (F > 0.8)",
           colour = "black", hjust = 0, size = 3.5, fontface = "italic") +         
  labs(title = "X-chromosome F by self-reported sex",
       subtitle = "Default F cutoff",
       x = NULL, y = "X-chromosome F") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "none")

ggsave(file.path(out_dir, "Default_F_distribution_by_reported_sex_and_race.png"),
       p_strip, width = 9, height = 9, dpi = 150)



p_strip_hisp <- ggplot(dat, aes(x = recorded_sex, y = F, colour = recorded_sex)) +
  facet_wrap(~race_hisp) +
  geom_jitter(width = 0.2, size = 1.6, alpha = 0.7) +
  scale_color_manual(values = sex_colors) +
  geom_hline(yintercept = 0.2, linetype = "dotted", colour = "grey50", , linewidth = 0.7) +
  geom_hline(yintercept = 0.8, linetype = "dashed", colour = "red", linewidth = 0.7) +
  annotate("text", x = 0.7, y = 0.25, label = "Female cutoff (F < 0.2)",
           colour = "black", hjust = 0, size = 3.5, fontface = "italic") +
  annotate("text", x = 0.7, y = 0.85, label = "Male cutoff (F > 0.8)",
           colour = "black", hjust = 0, size = 3.5, fontface = "italic") +         
  labs(title = "X-chromosome F by self-reported sex",
       subtitle = "Default F cutoff",
       x = NULL, y = "X-chromosome F") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "none")

ggsave(file.path(out_dir, "Default_F_distribution_by_reported_sex_and_race+hisp.png"),
       p_strip_hisp, width = 9, height = 9, dpi = 150)

# Stacked histogram: shows which race groups occupy which F clusters.
p_hist <- ggplot(dat, aes(x = F, fill = race)) +
  geom_histogram(binwidth = 0.02, colour = "grey30", linewidth = 0.1) +
  geom_vline(xintercept = c(0.2, 0.8), linetype = "dashed", colour = "grey40") +
  scale_fill_manual(values = race_colors) +
  annotate("text", x = -0.25, y = 100, label = "Female cutoff (F < 0.2)",
           colour = "black", hjust = 0, size = 3.5, fontface = "italic") +
  annotate("text", x = 0.4, y = 100, label = "Male cutoff (F > 0.8)",
           colour = "black", hjust = 0, size = 3.5, fontface = "italic") +
  labs(title = "X-chromosome inbreeding F by self-reported race",
       subtitle = "Female clusters (low F) split by ancestry; males at F~1. Dashed = sex-call thresholds",
       x = "X-chromosome F", y = "Count",
       fill = "Self-reported race") +
  theme_bw(base_size = 12)

ggsave(file.path(out_dir, "Default_Xsex_F_hist_by_reported_race.png"),
       p_hist, width = 10, height = 6, dpi = 150)

## ---- 3. Plots with new threshold -----------------------------
# swarm plot faceted by race
dat <- dat %>%
  mutate(plot_group = ifelse(mismatch_new, "mismatch", recorded_sex))

p_strip <- ggplot(dat, aes(x = recorded_sex, y = F, colour = plot_group)) +
  facet_wrap(~race) +
  geom_jitter(width = 0.2, size = 1.6, alpha = 0.7) +
  scale_color_manual(values = c(sex_colors, mismatch = "red")) +
  geom_hline(yintercept = 0.5, linetype = "dotted", colour = "grey50", linewidth = 0.7) +
  geom_hline(yintercept = 0.8, linetype = "dashed", colour = "red", linewidth = 0.7) +
  annotate("text", x = 0.7, y = 0.25, label = "Female cutoff (F < 0.5)",
           colour = "black", hjust = 0, size = 3.5, fontface = "italic") +
  annotate("text", x = 0.7, y = 0.85, label = "Male cutoff (F > 0.8)",
           colour = "black", hjust = 0, size = 3.5, fontface = "italic") +
  labs(title = "X-chromosome F by self-reported sex",
       subtitle = "Default F cutoff",
       x = NULL, y = "X-chromosome F") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "none")

ggsave(file.path(out_dir, "New_thresh_F_distribution_by_reported_sex_and_race.png"),
       p_strip, width = 9, height = 9, dpi = 150)



p_strip_hisp <- ggplot(dat, aes(x = recorded_sex, y = F, colour = plot_group)) +
  facet_wrap(~race_hisp) +
  geom_jitter(width = 0.2, size = 1.6, alpha = 0.7) +
  scale_color_manual(values = c(sex_colors, mismatch = "red")) +
  geom_hline(yintercept = 0.5, linetype = "dotted", colour = "grey50", , linewidth = 0.7) +
  geom_hline(yintercept = 0.8, linetype = "dashed", colour = "red", linewidth = 0.7) +
  annotate("text", x = 0.7, y = 0.25, label = "Female cutoff (F < 0.5)",
           colour = "black", hjust = 0, size = 3.5, fontface = "italic") +
  annotate("text", x = 0.7, y = 0.85, label = "Male cutoff (F > 0.8)",
           colour = "black", hjust = 0, size = 3.5, fontface = "italic") +         
  labs(title = "X-chromosome F by self-reported sex",
       subtitle = "Default F cutoff",
       x = NULL, y = "X-chromosome F") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "none")

ggsave(file.path(out_dir, "New_thresh_F_distribution_by_reported_sex_and_race+hisp.png"),
       p_strip_hisp, width = 9, height = 9, dpi = 150)

# Stacked histogram: shows which race groups occupy which F clusters.
p_hist <- ggplot(dat, aes(x = F, fill = race)) +
  geom_histogram(binwidth = 0.02, colour = "grey30", linewidth = 0.1) +
  geom_vline(xintercept = c(0.5, 0.8), linetype = "dashed", colour = "grey40") +
  scale_fill_manual(values = race_colors) +
  annotate("text", x = 0.0, y = 150, label = "Female cutoff (F < 0.5)",
           colour = "black", hjust = 0, size = 3.5, fontface = "italic") +
  annotate("text", x = 0.4, y = 100, label = "Male cutoff (F > 0.8)",
           colour = "black", hjust = 0, size = 3.5, fontface = "italic") +
  labs(title = "X-chromosome inbreeding F by self-reported race",
       subtitle = "Female clusters (low F) split by ancestry; males at F~1. Dashed = sex-call thresholds",
       x = "X-chromosome F", y = "Count",
       fill = "Self-reported race") +
  theme_bw(base_size = 12)

ggsave(file.path(out_dir, "New_thresh_Xsex_F_hist_by_reported_race.png"),
       p_hist, width = 10, height = 6, dpi = 150)


## Genuine mismatches to drop
samples_to_drop_sex <- dat %>% filter(mismatch_new) %>% mutate(sample.id = FID)

write.table(samples_to_drop_sex, '~/data/pca_results/QC/suggested_drop_sex_mismatch.txt')

samples_to_drop_sex_flag <- dat %>% 
                        mutate(sample.id = FID, sex_mismatch_flag = mismatch_new,
                         reported_sex = sex, genetic_sex =  tolower(genetic_sex_new_threshold)) %>%
                         select(sample.id, reported_sex, genetic_sex, sex_mismatch_flag)

write.table(samples_to_drop_sex_flag, '~/data/pca_results/QC/Sex_mismatch_flag.txt')

summary(dat[dat$F < 0.5, 'F'])

###Checking the mismatches 

library(xlsx)
sex_check_excel <- read.xlsx('/home/ec2-user/data/metadata/AMP-AD_DivCO-WGS/AMP_AD_WGS_sex_check_updated_260706.xlsx',  sheetIndex = 1)
sample_map <- read.table('/home/ec2-user/data/metadata/AMP-AD_DivCO-WGS/DivCo_SampleMapping.csv', sep = ',', header = T)
sample_map$wgs_specimenID <- sample_map$final_specimenid
sample_map$wgs_specimenID[grep('DLPFC', sample_map$wgs_specimenID)] <- paste0(sample_map$wgs_specimenID[grep('DLPFC', sample_map$wgs_specimenID)], '_WGS')
sex_check_excel$wgs_specimenID <- sample_map$wgs_specimenID[match(sex_check_excel$Sample, sample_map$filename.sample)]

View(sex_check_excel[sex_check_excel$wgs_specimenID %in% samples_to_drop_sex$FID,])
View(samples_to_drop_sex[samples_to_drop_sex$FID %in% sex_check_excel$wgs_specimenID,])

