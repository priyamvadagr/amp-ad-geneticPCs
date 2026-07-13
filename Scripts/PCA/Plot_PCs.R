################################################################################
# plot_pca_by_ancestry.R
# Plot principal components coloured by self-reported race/ethnicity to confirm
# the PCs capture ancestry (not batch), and that the top ~4 PCs (from the scree
# elbow) separate the ancestry groups.
#
# Inputs:
#   pca_eigenvectors.csv : per-sample PC scores (from run_pca.R)
#   individual metadata  : sample.id -> race, isHispanic, cohort/batch
################################################################################

library(dplyr)
library(ggplot2)
library(tidyr)
library(pals)

fig_dir <- "/home/ec2-user/AMP-AD_genetic_PCs/Results/PCA"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

## ---- 1. Load PC scores -----------------------------------------------------
pc <- read.csv("~/data/pca_results/PCA/pca_eigenvectors.csv",
               stringsAsFactors = FALSE)
pc$sample.id <- as.character(pc$sample.id)

## ---- 2. Attach ancestry metadata -------------------------------------------
# Map specimen (sample.id) -> individualID -> race / isHispanic.
biospec <- read.table('~/data/processed/metadata/Combined/DivCo_AMPAD_1.0_biospecimen_metadata.csv',
                      sep = ',', header = TRUE)
ind     <- read.table('~/data/processed/metadata/Combined/DivCo_AMPAD_1.0_individual_metadata.csv',
                      sep = ',', header = TRUE)

meta <- biospec %>% select(specimenID, individualID) %>% distinct() %>%
  left_join(ind %>% select(individualID, race, isHispanic, cohort, dataContributionGroup), by = "individualID")

pc <- pc %>%
  left_join(meta, by = c("sample.id" = "specimenID")) %>% 
  mutate(
    race = ifelse(is.na(race) | race == "" | race == "missing or unknown", "Missing/Unknown",
                  ifelse(race == "Other", "Not specified", race)),
    hisp = case_when(
      as.character(isHispanic) %in% c("True") ~ "Hispanic",
      as.character(isHispanic) %in% c("False") ~ "Non-Hispanic",
      TRUE ~ "unknown")
  )

## ---- 3. Colour palette for race --------------------------------------------
races <- unique(pc$race)
race_cols <- c("White" = "#CC79A7",  # orange
               "Asian" = "#0039A6",  # dark blue
               "Black or African American" = "#009E73",  # green
               "Not specified" = "#85144b",  # pink/magenta
               "American Indian or Alaska Native" = "#D55E00")  # vermillion/red-orange  

race_cols["Missing/Unknown"] <- "#374057"   # mute the unknowns
hisp_shapes <- c("Hispanic" = 17, "Non-Hispanic" = 16, "Unknown" = 4)

## ---- 4. Helper to plot a PC pair -------------------------------------------
plot_pc_pair <- function(df, xpc, ypc, colour_var = "race",
                         colours = race_cols, shape_var = NULL, title = NULL) {
  aes_args <- aes(x = .data[[xpc]], y = .data[[ypc]], colour = .data[[colour_var]])
  if (!is.null(shape_var)) aes_args <- modifyList(aes_args, aes(shape = .data[[shape_var]]))
  p <- ggplot(df, aes_args) +
    geom_point(size = 2, alpha = 0.8) +
    labs(title = title %||% paste(xpc, "vs", ypc, "by", colour_var),
         x = xpc, y = ypc, colour = colour_var, shape = shape_var) +
    theme_bw(base_size = 12)
  if (colour_var == "race") p <- p + scale_colour_manual(values = colours)
  p
}
`%||%` <- function(a, b) if (is.null(a)) b else a

## ---- 5. PC1-PC2 and PC3-PC4 coloured by race -------------------------------
p12 <- plot_pc_pair(pc, "PC1", "PC2", "race",
                    title = "PC1 vs PC2 by self-reported race")
ggsave(file.path(fig_dir, "pca_PC1_PC2_byRace.png"), p12,
       width = 12, height = 10, dpi = 150,  scale = 1)

p34 <- plot_pc_pair(pc, "PC3", "PC4", "race", 
                    title = "PC3 vs PC4 by self-reported race")
ggsave(file.path(fig_dir, "pca_PC3_PC4_byRace.png"), p34,
       width = 12, height = 10, dpi = 150,  scale = 1)

p56 <- plot_pc_pair(pc, "PC5", "PC6", "race", 
                    title = "PC5 vs PC6 by self-reported race")
ggsave(file.path(fig_dir, "pca_PC5_PC6_byRace.png"), p56,
       width = 12, height = 10, dpi = 150,  scale = 1)
p78 <- plot_pc_pair(pc, "PC7", "PC8", "race", 
                    title = "PC7 vs PC8 by self-reported race")
ggsave(file.path(fig_dir, "pca_PC7_PC8_byRace.png"), p78,
       width = 12, height = 10, dpi = 150,  scale = 1)


## ---- 6. PC1-PC2 with Hispanic ethnicity as shape ---------------------------

# shape for Hispanic ethnicity: filled circle vs filled triangle
p12h <- ggplot(pc, aes(x = PC1, y = PC2, colour = race, shape = hisp)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_colour_manual(values = race_cols, name = "Race") +
  scale_shape_manual(values = hisp_shapes, name = "Ethnicity") +
  labs(title = "PC1 vs PC2 by self-reported Race and Hispanic ethnicity",
       x = "PC1", y = "PC2") +
  theme_bw(base_size = 14) +
  guides(colour = guide_legend(override.aes = list(shape = 16)))  # legend dots
ggsave(file.path(fig_dir, "pca_PC1_PC2_byRace_hisp.png"), p12h,
       width = 12, height = 10, dpi = 150,  scale = 1)

# ---- 6. PC1-PC2 with Cohort and datacollectionCenter as color ---------------------------
cohorts <- unique(pc$cohort)
cohort_cols <- setNames(as.vector(pals::kelly(length(cohorts) + 2)[-c(1,2)]), cohorts)
p12cohort <- ggplot(pc, aes(x = PC1, y = PC2, colour = cohort)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_colour_manual(values = cohort_cols, name = "Cohorts") +
  labs(title = "PC1 vs PC2 by Cohort",
       x = "PC1", y = "PC2") +
  theme_bw(base_size = 14) +
  guides(colour = guide_legend(override.aes = list(shape = 16)))  # legend dots
ggsave(file.path(fig_dir, "pca_PC1_PC2_bycohort.png"), p12cohort,
       width = 12, height = 10, dpi = 150,  scale = 1)


data_col_center <- unique(pc$dataContributionGroup)
data_col_center_cols <- setNames(as.vector(pals::kelly(length(data_col_center) + 2)[-c(1,2)]), data_col_center)
p12data <- ggplot(pc, aes(x = PC1, y = PC2, colour = dataContributionGroup)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_colour_manual(values = data_col_center_cols, name = "Data Contribution Centers") +
  labs(title = "PC1 vs PC2 by Data Contribution Center",
       x = "PC1", y = "PC2") +
  theme_bw(base_size = 14) +
  guides(colour = guide_legend(override.aes = list(shape = 16)))  # legend dots
ggsave(file.path(fig_dir, "pca_PC1_PC2_bydatacenter.png"), p12data,
       width = 12, height = 10, dpi = 150,  scale = 1)


## ---- 7. Pairs panel: PC1-PC4 (the informative PCs from the scree) ----------
# long format so all informative PC pairs can be eyeballed in one figure
pc_pairs <- bind_rows(
  pc %>% transmute(race, xv = PC1, yv = PC2, pair = "PC1 vs PC2"),
  pc %>% transmute(race, xv = PC2, yv = PC3, pair = "PC2 vs PC3"),
  pc %>% transmute(race, xv = PC3, yv = PC4, pair = "PC3 vs PC4")
)
p_panel <- ggplot(pc_pairs, aes(xv, yv, colour = race)) +
  geom_point(size = 1.4, alpha = 0.75) +
  scale_colour_manual(values = race_cols) +
  facet_wrap(~ pair, scales = "free", nrow = 1) +
  labs(title = "Informative PCs (1-4) by self-reported race",
       x = NULL, y = NULL, colour = "Race") +
  theme_bw(base_size = 11)
ggsave(file.path(fig_dir, "pca_PC1-6_panel_byRace.png"), p_panel,
       width = 15, height = 5, dpi = 150)

cat("PCA-by-ancestry plots written to", fig_dir, "\n")
cat("Check: do the ancestry groups separate along PC1-PC4? If yes, the PCs\n")
cat("capture ancestry (not batch). PC5+ being noise (per scree) means ~4\n")
cat("informative PCs. Also confirm no isolated outliers far from all clusters.\n")