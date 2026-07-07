################################################################################
# lea_tracy_widom.R
# Determine the number of significant PCs using LEA's Tracy-Widom test.
#
# LEA::tracy.widom() needs a pcaProject object from LEA::pca(), which runs on a
# .lfmm genotype matrix (individuals x SNPs; entries 0/1/2; 9 = missing).
# So we: (1) extract the QC-passed, LD-pruned genotype matrix from the GDS,
#        (2) write it as .lfmm, (3) run LEA::pca(), (4) run tracy.widom().
#
# Install (Bioconductor):
#   if (!require("BiocManager")) install.packages("BiocManager")
#   BiocManager::install("LEA")
################################################################################

if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install("LEA")

library(SNPRelate)
library(LEA)
library(dplyr)

out_dir <- "~/data/pca_results/PCA"
lea_dir <- file.path(path.expand(out_dir), "LEA")
dir.create(lea_dir, showWarnings = FALSE, recursive = TRUE)

## ---- 1. Pull QC-passed samples + pruned SNPs -------------------------------
pruned <- readRDS("~/data/pca_results/QC/pruned_snp_ids.rds")
flags  <- read.csv("~/data/pca_results/QC/cohort_keep_drop_flags.csv",
                   stringsAsFactors = FALSE)
keep_samples <- flags %>% filter(!drop) %>% pull(sample.id) %>% as.character()

gds_file <- "~/data/processed/DivCo_AMP1.0_GRCh38_1019Samples_filtered.gds"
genofile <- snpgdsOpen(gds_file)

gds_samples  <- read.gdsn(index.gdsn(genofile, "sample.id"))
keep_samples <- intersect(keep_samples, gds_samples)
cat("QC-passed samples:", length(keep_samples), "; pruned SNPs:", length(pruned), "\n")

## ---- 2. Extract genotype matrix (samples x SNPs), recode for LEA -----------
# snpgdsGetGeno returns counts of the reference allele (0/1/2), NA for missing.
geno <- snpgdsGetGeno(genofile, sample.id = keep_samples, snp.id = pruned,
                      with.id = TRUE)
G <- geno$genotype                  # rows = samples, cols = SNPs
# compute variance per SNP ignoring missing (NA) values
snp_var <- apply(G, 2, function(x) var(x, na.rm = TRUE))
nonconstant <- which(snp_var > 0 & !is.na(snp_var))
cat("Dropping", ncol(G) - length(nonconstant),
    "constant/zero-variance SNPs;", length(nonconstant), "remain\n")
G <- G[, nonconstant]
pruned_kept <- geno$snp.id[nonconstant]   # keep track of which SNPs remain

# now recode missing and write
G[is.na(G)] <- 9

## ---- 3. Write .lfmm and run LEA's PCA --------------------------------------
lfmm_file <- file.path(lea_dir, "qc_pruned_geno.lfmm")
LEA::write.lfmm(G, lfmm_file)

# LEA's PCA (scale = TRUE standardises each SNP, as in standard genetic PCA)
pc <- LEA::pca(lfmm_file, scale = TRUE)

## ---- 4. Tracy-Widom test ---------------------------------------------------
tw <- LEA::tracy.widom(pc)

tw_df <- data.frame(
  PC          = seq_along(tw$pvalues),
  eigenvalue  = tw$eigenvalues[seq_along(tw$pvalues)],
  tw_stat     = tw$twstats,
  p_value     = tw$pvalues,
  pct_variance = tw$percentage[seq_along(tw$pvalues)] * 100
)

first_ns <- which(tw_df$p_value >= 0.05)[1]   # first non-significant PC
n_sig    <- first_ns - 1
n_sig
cat("\n=== LEA Tracy-Widom test ===\n")
print(head(tw_df, 15))
cat("\nSignificant PCs (p < 0.05):", n_sig, "\n")

write.csv(tw_df, file.path(out_dir, "lea_tracy_widom_results.csv"), row.names = FALSE)

## ---- 5. Cross-check: % variance vs your SNPRelate scree --------------------
cat("\nTop-10 % variance (LEA) — compare to your SNPRelate scree:\n")
print(round(tw_df$pct_variance[1:10], 3))
cat("\nThese should closely match the SNPRelate PCA (same genotypes).\n")
cat("Use the TW count alongside the scree elbow (~4) and PC scatter plots.\n")