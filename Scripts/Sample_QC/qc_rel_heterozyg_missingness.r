################################################################################
# qc_relatedness.R
# Full sample QC for the DivCo + AMP-AD 1.0 genotype set, prior to PCA.
#
# Steps (designed to be run INTERACTIVELY, block by block):
#   1. Sample call rate (missingness per individual)
#   2. Heterozygosity outliers
#   3. KING-robust relatedness -> flag duplicates (~0.5) and relatives
#
# Relatives/duplicates are only FLAGGED here; you decide who to drop.
################################################################################

## ---- 0. Setup --------------------------------------------------------------
library(SNPRelate)
library(dplyr)
library(ggplot2)

gds_file <- "~/data/processed/DivCo_AMP1.0_GRCh38_1019Samples_filtered.gds"
out_dir  <- "~/data/pca_results/QC"
plots_dir <- "~/AMP-AD_genetic_PCs/Results/Filtering"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
set.seed(1)   # LD pruning is stochastic; seed for reproducibility
genofile <- snpgdsOpen(gds_file)

# sanity check: dimensions
snpgdsSummary(gds_file)
all_samples <- read.gdsn(index.gdsn(genofile, "sample.id"))
length(all_samples)

# read in ids of filtered snps 
pruned <- readRDS(paste0(out_dir, '/', 'pruned_snp_ids.rds'))

## ---- 1. Sample call rate (per-individual missingness) ----------------------
# Flag individuals with high missingness (poorly genotyped samples).
samp_miss <- snpgdsSampMissRate(genofile, snp.id = pruned)
samp_miss_df <- data.frame(sample.id = all_samples, miss.rate = samp_miss)

miss_threshold <- 0.05   # flag samples missing >5% of pruned SNPs
samp_miss_df$flag_missing <- samp_miss_df$miss.rate > miss_threshold

summary(samp_miss_df$miss.rate)
cat("Samples flagged for high missingness:",
    sum(samp_miss_df$flag_missing), "\n")


# inspect the distribution
png(file = paste0(plots_dir, '/', 'sample_missingness.png'))
print(
samp_miss_df |>
  ggplot(aes(x = miss.rate)) +
  geom_histogram(bins = 50, fill = "#4E79A7",
                 color = "white", linewidth = 0.3) +
  geom_vline(xintercept = 0.01, linetype = "dashed", color = "#E15759", linewidth = 0.8) +
  annotate("text", x = 0.01, y = Inf, label = "1% threshold", 
           hjust = -0.1, vjust = 2, color = "#E15759", size = 4) +
  scale_x_continuous(labels = scales::percent_format()) +
  labs(title = "Per-sample missing rate",
       x = "Missing rate",
       y = "Number of samples") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank()))
dev.off()

# save the missingness information 
write.table(samp_miss_df, file = '~/data/pca_results/QC/Sample_missingness.txt')


## ---- 2. Heterozygosity outliers --------------------------------------------
# Excess or deficit heterozygosity flags contamination (high het) or
# inbreeding/poor quality (low het). Flag samples > 'het_sd' SDs from the mean.
het <- snpgdsIndInb(genofile, snp.id = pruned, method = "mom.visscher")
het_df <- data.frame(sample.id = all_samples, inbreeding = het$inbreeding)

het_sd  <- 4
het_mean <- mean(het_df$inbreeding, na.rm = TRUE)
het_sdev <- sd(het_df$inbreeding, na.rm = TRUE)
het_df$flag_het <- abs(het_df$inbreeding - het_mean) > het_sd * het_sdev

cat("Heterozygosity (inbreeding coeff) summary:\n"); summary(het_df$inbreeding)
cat("Samples flagged as het outliers (>", het_sd, "SD):",
    sum(het_df$flag_het), "\n")


pdf(file = paste0(plots_dir, '/', 'sample_heterozygosity.pdf'))
print(
het_df |>
  ggplot(aes(x = inbreeding)) +
  geom_histogram(bins = 50, fill = "#4E79A7",
                 color = "white", linewidth = 0.3) +
  labs(title = "Per-sample Heterozygosity",
       x = "Inbreeding coefficient (F)",
       y = "Number of samples") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank()))
dev.off()

# save heterozygosity information 
write.table(het_df, file = '~/data/pca_results/QC/Heterozygosity_flag.txt')

## ---- 4. KING-robust relatedness --------------------------------------------
# KING-robust is robust to population structure -> appropriate for this
# ancestrally diverse cohort. Critical given DivCo + AMP-AD 1.0 merges
# OVERLAPPING donors: duplicates appear at kinship ~0.5.
king <- snpgdsIBDKING(genofile, snp.id = pruned, autosome.only = TRUE,
                      num.thread = 2)

# pairwise kinship table
kin <- snpgdsIBDSelection(king)   # cols: ID1, ID2, IBS0, kinship

# Relationship thresholds (standard):
#   > 0.354  : duplicate / monozygotic (also same individual under 2 IDs)
#   0.177-0.354 : 1st-degree (parent-offspring, full sibs)
#   0.0884-0.177: 2nd-degree
#   0.0442-0.0884: 3rd-degree
#   < 0.0442 : unrelated
kin$relationship <- cut(
  kin$kinship,
  breaks = c(-Inf, 0.0442, 0.0884, 0.177, 0.354, Inf),
  labels = c("unrelated","3rd-degree","2nd-degree","1st-degree","dup/MZ")
)

# Flagged pairs (3rd-degree or closer) -- for MANUAL review
related_pairs <- kin %>% filter(kinship > 0.0442) %>% arrange(desc(kinship))
dup_pairs     <- kin %>% filter(kinship > 0.354) %>% arrange(desc(kinship))

cat("Total pairs evaluated:", nrow(kin), "\n")
cat("Related pairs (3rd-degree+):", nrow(related_pairs), "\n")
cat("Duplicate/MZ pairs (kinship>0.354):", nrow(dup_pairs), "\n")

# inspect the distribution and the flagged pairs
hist(kin$kinship, breaks = 100, main = "Pairwise KING kinship",
     xlab = "kinship coefficient")
print(head(related_pairs, 30))
print(dup_pairs)   # these are the likely overlapping-donor duplicates

# write the flagged pairs out for manual review
write.csv(related_pairs, file.path(out_dir, "Sample_relatedness_flag.csv"),
          row.names = FALSE)
write.csv(dup_pairs, file.path(out_dir, "duplicate_pairs_flagged.csv"),
          row.names = FALSE)

## ---- 5. Assemble QC summary + cleaned sample list --------------------------
# Combine the per-sample flags. NOTE: relatives/duplicates are NOT auto-removed
# here -- you review related_pairs_flagged.csv and decide which member to drop.
qc <- samp_miss_df %>%
  left_join(het_df, by = "sample.id")

# samples appearing in any 2nd degree related/duplicate pair (for your reference)

related_pairs <- related_pairs[related_pairs$relationship != '3rd-degree', ]
flagged_related_ids <- unique(c(related_pairs$ID1, related_pairs$ID2))
qc$in_related_pair <- qc$sample.id %in% flagged_related_ids

cat("\n--- QC SUMMARY ---\n")
cat("Total samples:               ", nrow(qc), "\n")
cat("Fail call-rate (>5% miss):   ", sum(qc$flag_missing), "\n")
cat("Het outliers (>4 SD):        ", sum(qc$flag_het), "\n")
cat("Fail technical QC (either):  ", sum(qc$fail_technical), "\n")
cat("In a related/duplicate pair: ", sum(qc$in_related_pair), "\n")

