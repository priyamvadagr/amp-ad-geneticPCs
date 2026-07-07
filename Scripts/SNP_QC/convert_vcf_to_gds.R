#!/usr/bin/env Rscript
# convert_vcf_to_gds.R
# Convert filtered VCF to GDS format for SNPRelate PCA

library(SNPRelate)

# --- Paths ---
vcf_file <- "~/data/processed/DivCo_AMP1.0_GRCh38_1019Samples_filtered.vcf.gz"
gds_file <- "~/data/processed/DivCo_AMP1.0_GRCh38_1019Samples_filtered.gds"

# --- Convert ---
snpgdsVCF2GDS(
  vcf.fn      = vcf_file,
  out.fn      = gds_file,
  method      = "biallelic.only",   # keep only biallelic SNPs
  snpfirstdim = FALSE,              # samples in rows; standard for SNPRelate
  verbose     = TRUE
)

# --- Verify ---
snpgdsSummary(gds_file)