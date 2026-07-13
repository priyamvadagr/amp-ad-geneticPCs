################################################################################
# detect_PCA_outliers.R
# Detect PCA outlier samples using the KNN-based Probabilistic LOF statistic
# from Prive et al. 2020 (Bioinformatics), implemented in bigutilsr.
# Threshold is chosen by VISUAL INSPECTION of
# the statistic's histogram + PC scatter coloured by the statistic (as the
# paper recommends), NOT an automatic cutoff.
#
################################################################################

library(bigutilsr)
library(ggplot2)
library(dplyr)

fig_dir <- "/home/ec2-user/AMP-AD_genetic_PCs/Results/PCA/Outliers"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

## ---- 1. PC scores (use the informative PCs from your scree, ~4) -------------
pc <- read.csv("~/data/pca_results/PCA/pca_eigenvectors.csv",
               stringsAsFactors = FALSE)
pc$sample.id <- as.character(pc$sample.id)

# use the informative PCs only (scree elbow ~4). Adjust n_pc if you keep more.
n_pc <- 4
U <- as.matrix(pc[, paste0("PC", seq_len(n_pc))])

## ---- 2. Probabilistic LOF statistic (KNN-based) ----------------------------
# prob_dist returns per-point KNN distance statistics; the LOF-style ratio is
# in $lof. Uses K nearest neighbours (default ~ min(30, n/2)); set explicitly.
K <- 30
llof <- bigutilsr::prob_dist(U, ncores = 4)     # KNN prob-dist object
stat <- llof$dist.self / llof$dist.nn           # pd_j / mean(pd of neighbours)
# (bigutilsr also exposes LOF directly via LOF(); prob_dist gives the components)

# The paper's statistic already incorporates the sqrt; bigutilsr's LOF() wraps it.
# Simpler + matches the paper: use bigutilsr::LOF() on the PCs.
lof <- bigutilsr::LOF(U, seq_k = K)  # log-LOF statistic per sample
pc$lof <- lof
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

# An automatic starting suggestion via Tukey's rule (skewness+MT-adjusted),
# also from bigutilsr — use as a GUIDE, then confirm visually.
thr <- bigutilsr::tukey_mc_up(pc$lof)
cat("Suggested LOF threshold (tukey_mc_up):", round(thr, 3), "\n")
thr <- 2 # picked 2 based on visual inspection
pc$outlier <- pc$lof > thr
cat("Outliers above suggested threshold:", sum(pc$outlier), "\n")
print(pc %>% filter(outlier) %>% select(sample.id, lof))

## ---- 3. VISUAL threshold selection (histogram) -----------------------------
# The paper: choose the threshold by eye from the histogram + coloured PC plots.

p_hist <- ggplot(pc, aes(lof)) +
  geom_histogram(bins = 60, fill = "grey70", colour = "grey40") +
  geom_vline(xintercept = thr, color = 'red', linetype = 'dashed') +
  labs(title = "Distribution of the (log) Local Outlier Factor statistic",
       subtitle = "Choose the outlier threshold by eye: look for the right-tail break",
       x = "log-LOF statistic", y = "count") +
  theme_bw(base_size = 12)
ggsave(file.path(fig_dir, "pca_LOF_histogram.png"), p_hist,
       width = 8, height = 5, dpi = 150)



## ---- 4. PC scatter coloured by the statistic (paper's Fig. 2 style) --------
races <- unique(pc$race)
race_cols <- c("White" = "#CC79A7",  # orange
               "Asian" = "#0039A6",  # dark blue
               "Black or African American" = "#009E73",  # green
               "Not specified" = "#85144b",  # pink/magenta
               "American Indian or Alaska Native" = "#D55E00")  # vermillion/red-orange  


race_cols["Missing/Unknown"] <- "#374057"   # mute the unknowns
hisp_shapes <- c("Hispanic" = 17, "Non-Hispanic" = 16, "Unknown" = 4)

p_scores <- ggplot(pc, aes(PC5, PC6, colour = lof)) +
  geom_point(size = 1.8) +
  scale_colour_viridis_c(option = "plasma") +
  labs(title = "PC1 vs PC2 coloured by LOF statistic",
       subtitle = "High-LOF points (bright) sitting apart from clusters are candidate outliers",
       colour = "log-LOF") +
  theme_bw(base_size = 12)
ggsave(file.path(fig_dir, "pca_PC1_PC2_byLOF.png"), p_scores,
       width = 8, height = 6, dpi = 150)

p_scores_pc34 <- ggplot(pc, aes(PC3, PC4, colour = lof)) +
  geom_point(size = 1.8) +
  scale_colour_viridis_c(option = "plasma") +
  labs(title = "PC3 vs PC4 coloured by LOF statistic",
       subtitle = "High-LOF points (bright) sitting apart from clusters are candidate outliers",
       colour = "log-LOF") +
  theme_bw(base_size = 12)
ggsave(file.path(fig_dir, "pca_PC3_PC4_byLOF.png"), p_scores_pc34,
       width = 8, height = 6, dpi = 150)


p_score_cutoff_34 <- ggplot(pc, aes(PC3, PC4, colour = lof > thr)) +
  geom_point(size = 1.8) +
  scale_colour_viridis_d(option = "plasma",
                         labels = c(paste0("LOF \u2264 ", thr), paste0("LOF > ", thr)),
                         name = "Outlier status") +
  labs(title = "PC3 vs PC4 coloured by LOF threshold",
       x = "PC3", y = "PC4") +
  theme_bw(base_size = 12)
ggsave(file.path(fig_dir, "pca_PC3_PC4_abovecutOff.png"), p_score_cutoff_34,
       width = 8, height = 6, dpi = 150)

p_score_cutoff_12 <- ggplot(pc, aes(PC1, PC2, colour = lof > thr)) +
  geom_point(size = 1.8) +
  scale_colour_viridis_d(option = "plasma",
                         labels = c(paste0("LOF \u2264 ", thr), paste0("LOF > ", thr)),
                         name = "Outlier status") +
  labs(title = "PC1 vs PC2 coloured by LOF threshold",
       x = "PC1", y = "PC2") +
  theme_bw(base_size = 12)
ggsave(file.path(fig_dir, "pca_PC1_PC2_abovecutOff.png"), p_score_cutoff_12,
       width = 8, height = 6, dpi = 150)

p_score_cutoff_56 <- ggplot(pc, aes(PC5, PC6, colour = lof > thr)) +
  geom_point(size = 1.8) +
  scale_colour_viridis_d(option = "plasma",
                         labels = c(paste0("LOF \u2264 ", thr), paste0("LOF > ", thr)),
                         name = "Outlier status") +
  labs(title = "PC5 vs PC6 coloured by LOF threshold",
       x = "PC5", y = "PC6") +
  theme_bw(base_size = 12)
ggsave(file.path(fig_dir, "pca_PC5_PC6_abovecutOff.png"), p_score_cutoff_56,
       width = 8, height = 6, dpi = 150)


p_flag_pc1_pc2 <- ggplot(pc, aes(x = PC1, y = PC2, colour = race, shape = hisp)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_point(data = subset(pc, outlier), colour = "black", size = 4, shape = 1) +
  scale_colour_manual(values = race_cols, name = "Self-reported Race") +
  scale_shape_manual(values = hisp_shapes, name = "Hispanic Ethnicity") +
  labs(title = "PC2 vs PC3 by self-reported Race and Hispanic ethnicity",
       x = "PC1", y = "PC2") +
  theme_bw(base_size = 14) +
  guides(colour = guide_legend(override.aes = list(shape = 16)))  # legend dots
  ggsave(file.path(fig_dir, "pca_PC1_PC2_outliers.png"), p_flag_pc1_pc2,
       width = 12, height = 10, dpi = 150)

p_flag_pc3_pc4 <- ggplot(pc, aes(x = PC3, y = PC4, colour = race, shape = hisp)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_point(data = subset(pc, outlier), colour = "black", size = 4, shape = 1) +
  scale_colour_manual(values = race_cols, name = "Self-reported Race") +
  scale_shape_manual(values = hisp_shapes, name = "Hispanic Ethnicity") +
  labs(title = "PC3 vs PC4 by self-reported Race and Hispanic ethnicity",
       x = "PC3", y = "PC4") +
  theme_bw(base_size = 14) +
  guides(colour = guide_legend(override.aes = list(shape = 16)))  # legend dots
  ggsave(file.path(fig_dir, "pca_PC3_PC4_outliers.png"), p_flag_pc3_pc4,
       width = 12, height = 10, dpi = 150)

p_flag_pc5_pc6 <- ggplot(pc, aes(x = PC5, y = PC6, colour = race, shape = hisp)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_point(data = subset(pc, outlier), colour = "black", size = 4, shape = 1) +
  scale_colour_manual(values = race_cols, name = "Self-reported Race") +
  scale_shape_manual(values = hisp_shapes, name = "Hispanic Ethnicity") +
  labs(title = "PC5 vs PC6 by self-reported Race and Hispanic ethnicity",
       x = "PC5", y = "PC6") +
  theme_bw(base_size = 14) +
  guides(colour = guide_legend(override.aes = list(shape = 16)))  # legend dots
  ggsave(file.path(fig_dir, "pca_PC5_PC6_outliers.png"), p_flag_pc5_pc6,
       width = 12, height = 10, dpi = 150)

write.csv(pc %>% select(sample.id, lof, outlier),
          "~/data/pca_results/PCA/pca_lof_outliers.csv", row.names = FALSE)

cat("\nOutputs in", fig_dir, "\n")
cat("Decide the final threshold from the histogram + coloured PC plots, then\n")
cat("re-set pc$outlier <- pc$lof > <your_threshold>. Confirm flagged points\n")
cat("sit APART from ancestry clusters (not within a genuine ancestry group).\n")
