################################################################################
# plot_het_outlier_qc.R
# Plot per-sample WGS QC metrics with the 10 heterozygosity-flagged samples
# highlighted against the full-cohort distribution.
#
# Goal: visually confirm whether the het outliers are anomalous on any
# technical metric (-> real QC problem) or sit within the normal spread
# (-> heterozygosity flag is likely ancestry, not artifact).
################################################################################

library(dplyr)
library(ggplot2)
library(tidyr)

out_dir <- "~/AMP-AD_genetic_PCs/Results/Filtering/Heterozygosity"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ---- 1. Identify the outlier samples ---------------------------------------
het_df <- read.table('~/data/pca_results/QC/Heterozygosity_flag.txt')
het_outlier_ids <- het_df[het_df$flag_het, "sample.id"]

amp_ad_wgs_df <- read.table('~/data/processed/metadata/DivCO/DivCo_assay_WGS_metadata.csv',
                            sep = ',', header = T)

## ---- 2. Build a plotting frame: metrics + outlier flag + het value ---------
# Numeric QC metrics to visualize
metrics <- c("Mean.Coverage", "Percent.Sequence.Contamination",
             "Percent.PF.Aligned", "Percent.Total.Dup.",
             "Strand.Balance", "AT.Dropout", "GC.Dropout",
             "Median.Insert.Size", "Concentration..ng.uL.")

# GQN is stored as character -> coerce to numeric for plotting
amp_ad_wgs_df$GQN <- suppressWarnings(as.numeric(amp_ad_wgs_df$GQN))
metrics <- c(metrics, "GQN")

# keep only samples that are in our genotyped/het table, tag the outliers,
# and attach the inbreeding coefficient for the F-vs-metric scatterplots
qc_plot <- amp_ad_wgs_df %>%
  filter(specimenID %in% het_df$sample.id) %>%
  mutate(is_outlier = specimenID %in% het_outlier_ids) %>%
  left_join(het_df, by = c("specimenID" = "sample.id"))

cat("Samples in plotting frame:", nrow(qc_plot),
    "| outliers:", sum(qc_plot$is_outlier), "\n")


#het_qc <- qc_plot[qc_plot$Percent.Sequence.Contamination >= 20 & qc_plot$flag_het == TRUE, ]

## ---- 3. Faceted distributions: each metric, outliers overlaid --------------
# Long format for faceting
long_df <- qc_plot %>%
  select(specimenID, is_outlier, all_of(metrics)) %>%
  pivot_longer(cols = all_of(metrics),
               names_to = "metric", values_to = "value")

# Violin/box of full distribution with outlier points highlighted
p_dist <- ggplot(long_df, aes(x = metric, y = value)) +
  geom_violin(fill = "grey90", colour = "grey60", scale = "width") +
  geom_jitter(data = subset(long_df, !is_outlier),
              width = 0.15, size = 0.6, alpha = 0.25, colour = "grey50") +
  geom_jitter(data = subset(long_df, is_outlier),
              width = 0.12, size = 2.2, colour = "#D6604D") +
  facet_wrap(~ metric, scales = "free", ncol = 3) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(title = "WGS QC metrics: het outliers (red) vs. full cohort",
       x = NULL, y = "value")

ggsave(file.path(out_dir, "het_outliers_metric_distributions.png"),
       p_dist, width = 11, height = 9, dpi = 150)

## ---- 4. Key scatter: inbreeding F vs. contamination & coverage -------------
# These two are the most diagnostic: contamination would show as negative F
# (excess het); coverage rules out low-depth artifacts.
p_contam <- ggplot(qc_plot,
                   aes(x = Percent.Sequence.Contamination, y = inbreeding)) +
  geom_point(aes(colour = is_outlier, size = is_outlier), alpha = 0.7) +
  scale_colour_manual(values = c(`FALSE` = "grey60", `TRUE` = "#D6604D")) +
  scale_size_manual(values = c(`FALSE` = 1, `TRUE` = 2.5)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  theme_bw(base_size = 12) +
  labs(title = "Inbreeding F vs. sequence contamination",
       subtitle = "Contamination would drive NEGATIVE F (excess heterozygosity)",
       x = "% Sequence contamination", y = "Inbreeding coefficient (F)",
       colour = "Het outlier", size = "Het outlier")

ggsave(file.path(out_dir, "F_vs_contamination.png"),
       p_contam, width = 7, height = 5.5, dpi = 150)

p_cov <- ggplot(qc_plot, aes(x = Mean.Coverage, y = inbreeding)) +
  geom_point(aes(colour = is_outlier, size = is_outlier), alpha = 0.7) +
  scale_colour_manual(values = c(`FALSE` = "grey60", `TRUE` = "#D6604D")) +
  scale_size_manual(values = c(`FALSE` = 1, `TRUE` = 2.5)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  theme_bw(base_size = 12) +
  labs(title = "Inbreeding F vs. mean coverage",
       x = "Mean coverage (X)", y = "Inbreeding coefficient (F)",
       colour = "Het outlier", size = "Het outlier")

ggsave(file.path(out_dir, "F_vs_coverage.png"),
       p_cov, width = 7, height = 5.5, dpi = 150)

## ---- 5. Tabular comparison: outliers vs. cohort medians --------------------
comparison <- long_df %>%
  group_by(metric) %>%
  summarise(
    cohort_median = median(value, na.rm = TRUE),
    cohort_min    = min(value, na.rm = TRUE),
    cohort_max    = max(value, na.rm = TRUE),
    outlier_min   = min(value[is_outlier], na.rm = TRUE),
    outlier_max   = max(value[is_outlier], na.rm = TRUE),
    .groups = "drop"
  )
print(comparison)
write.csv(comparison, file.path(out_dir, "het_outlier_metric_comparison.csv"),
          row.names = FALSE)

cat("\nPlots written to", out_dir, "\n")
cat("Interpretation: if red points sit WITHIN the grey distributions on every\n")
cat("technical metric, the het flags are likely ancestry, not artifact.\n")