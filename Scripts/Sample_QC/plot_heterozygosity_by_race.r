################################################################################
# plot_het_raw_byRace_hispanic.R
# Raw heterozygosity (inbreeding F) vs missing rate, colored by self-reported
# race AND distinguishing Hispanic ethnicity, with het-flagged samples ringed.
#
# Requires in session / on disk:
#   het_df, samp_miss_df  : sample-level F, flag_het, missingness
#   ind_metadata_df       : individual metadata incl. race + isHispanic
#   biospec_metadata_df   : maps individualID <-> specimenID (sample.id)
################################################################################

library(dplyr)
library(ggplot2)

out_dir <- "/home/ec2-user/AMP-AD_genetic_PCs/Results/Filtering/Heterozygosity"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

het_df       <- read.table('~/data/pca_results/QC/Heterozygosity_flag.txt')
samp_miss_df <- read.table('~/data/pca_results/QC/Sample_missingness.txt')

ind_metadata_df <- read.table(
  '~/data/processed/metadata/Combined/DivCo_AMPAD_1.0_individual_metadata.csv',
  sep = ',', header = TRUE)
biospec_metadata_df <- read.table(
  '~/data/processed/metadata/Combined/DivCo_AMPAD_1.0_biospecimen_metadata.csv',
  sep = ',', header = TRUE)



# map individualID -> specimenID (sample.id used in the genotype data)
ind_metadata_df <- ind_metadata_df |>
  mutate(sample.id = biospec_metadata_df$specimenID[
    match(individualID, biospec_metadata_df$individualID)
  ])

## ---- 0. Build a combined race + Hispanic ethnicity label -------------------
# Normalize isHispanic to a clean logical-ish flag
ind_metadata_df <- ind_metadata_df |>
  mutate(
    hisp = case_when(
      tolower(as.character(isHispanic)) %in% c("true","t","yes","1")  ~ TRUE,
      tolower(as.character(isHispanic)) %in% c("false","f","no","0")  ~ FALSE,
      TRUE                                                            ~ NA
    ),
    race_clean = ifelse(is.na(race) | race == "missing or unknown", "Missing/Unknown", race),
    # combined category: append "(Hispanic)" where isHispanic is TRUE
    ethnicity = case_when(
      is.na(hisp)       ~ race_clean,
      hisp              ~ paste0(race_clean, " (Hispanic)"),
      TRUE              ~ race_clean
    )
  )

## ---- 1. Assemble per-sample data -------------------------------------------
plot_dat <- het_df %>%
  left_join(samp_miss_df %>% select(sample.id, miss.rate), by = "sample.id") %>%
  left_join(ind_metadata_df %>% select(sample.id, race_clean, hisp, ethnicity),
            by = "sample.id") %>%
  mutate(
    ethnicity = ifelse(is.na(ethnicity), "Missing/Unknown", ethnicity),
    hisp      = ifelse(is.na(hisp), FALSE, hisp),
    het_raw   = 1 - inbreeding
  )
## ---- 2. Plot: colour = race, SHAPE encodes Hispanic vs not -----------------
# Colour by race (the base category); shape distinguishes Hispanic (triangle)
# from non-Hispanic (circle). Flagged samples get a black ring on top.

race_colors <- c("White" = "#E69F00",  # orange
               "Asian" = "#0039A6",  # dark blue
               "Black or African American" = "#009E73",  # green
               "Other" = "#85144b",  # pink/magenta
               "American Indian or Alaska Native" = "#D55E00")  # vermillion/red-orange,  

race_colors["Missing/Unknown"] <- "#374057"   # mute the unknowns

sex_colors <- c(
  Male   = "#5566AA",  # muted indigo
  Female = "#CC8800"   # ochre/gold
)
hisp_shapes <- c("Hispanic" = 17, "Non-Hispanic" = 16, "Unknown" = 4)
miss_breaks <- c(1e-4, 4e-4, 1.6e-3, 6.36e-3, 2.496e-2, 9.29e-2)

p_raw <- ggplot(plot_dat,
                aes(x = miss.rate, y = inbreeding,
                    colour = race_clean, shape = hisp)) +
  scale_colour_manual(values = race_colors, name = "Self-reported Race") +
  geom_point(size = 1.8, alpha = 0.75) +
  # ring the het-flagged samples so they stand out regardless of race/ethnicity
  geom_point(data = subset(plot_dat, flag_het),
             shape = 21, size = 3.4, stroke = 1.1,
             colour = "black", fill = NA) +
  scale_x_log10(breaks = miss_breaks,
                labels = formatC(miss_breaks, format = "g")) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17),
                     labels = c(`FALSE` = "Non-Hispanic", `TRUE` = "Hispanic"),
                     name = "Hispanic ethnicity") +
  labs(title = "Raw heterozygosity vs missing rate",
       subtitle = paste0("Colour = self-reported race; triangle = Hispanic; ",
                         sum(plot_dat$flag_het, na.rm = TRUE),
                         " het-flagged samples circled in black"),
       x = "Missing rate", y = "Inbreeding coefficient (F)",
       colour = "Self-reported race") +
  theme_bw(base_size = 12)

ggsave(file.path(out_dir, "het_raw_byRace_hispanic.png"),
       p_raw, width = 10, height = 6.5, dpi = 150)

## ---- 3. Quick tally so you can see the Hispanic breakdown ------------------
cat("Race x Hispanic counts:\n")
print(table(race = plot_dat$race_clean, hispanic = plot_dat$hisp, useNA = "ifany"))
cat("\nFlagged samples by ethnicity:\n")
print(table(plot_dat$ethnicity[plot_dat$flag_het]))

cat("\nPlot written to", file.path(out_dir, "het_raw_byRace_hispanic.png"), "\n")

#-------Plot by Percent.Sequence.Contamination for samples with this info---------
wgs_assay_metadata_df <- read.table(
  '~/data/processed/metadata/DivCO/DivCo_assay_WGS_metadata.csv', sep = ',', header = TRUE)


# Subset plot data to samples with contamination info 
plot_dat <- plot_dat %>% filter(sample.id %in% wgs_assay_metadata_df$specimenID) 

## ---- Attach sequence contamination -------------------------------------
# wgs_assay_df is keyed on specimenID, which matches sample.id here
plot_dat <- plot_dat %>%
  mutate(
    contamination = wgs_assay_metadata_df$Percent.Sequence.Contamination[
      match(sample.id, wgs_assay_metadata_df$specimenID)
    ]
  )

cat("Samples with contamination value:",
    sum(!is.na(plot_dat$contamination)), "of", nrow(plot_dat), "\n")


## ---- 4. Plot: raw heterozygosity vs percent contamination ------------------
p_contam <- ggplot(subset(plot_dat, !is.na(contamination)),
                   aes(x = contamination, y = inbreeding,
                       colour = race_clean, shape = hisp)) +
    scale_colour_manual(values = race_colors, name = "Self-reported Race") +                     
  geom_point(size = 1.8, alpha = 0.75) +
  geom_point(data = subset(plot_dat, flag_het & !is.na(contamination)),
             shape = 21, size = 3.4, stroke = 1.1,
             colour = "black", fill = NA) +
  geom_vline(xintercept = 2, linetype = "dashed", colour = "#fa2929") +
  annotate("text", x = 2.1, y = 0.1, label = "Contamination cutoff (> 2%)",
           colour = "black", hjust = 0, size = 3.5, fontface = "italic") +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17),
                     labels = c(`FALSE` = "Non-Hispanic", `TRUE` = "Hispanic"),
                     name = "Hispanic ethnicity") +
  labs(title = "Raw heterozygosity vs percent contamination",
       subtitle = paste0("Colour = self-reported race; triangle = Hispanic; ",
                         sum(plot_dat$flag_het & !is.na(plot_dat$contamination), na.rm = TRUE),
                         " het-flagged samples circled in black; dashed = 2% cutoff"),
       x = "Percent Sequence Contamination", y = "Inbreeding coefficient (F)",
       colour = "Self-reported race") +
  theme_bw(base_size = 12)

ggsave(file.path(out_dir, "het_raw_byContamination_hispanic.png"),
       p_contam, width = 15, height = 9, dpi = 150)

cat("\nPlot written to", file.path(out_dir, "het_raw_byContamination_hispanic.png"), "\n")
