################################################################################
# run_pca.R
# Population-structure PCA on the QC-passed, LD-pruned dataset.
#
# Inputs:
#   GDS file (filtered, autosomal, common biallelic SNPs)
#   pruned_snp_ids.rds   : LD-pruned SNP IDs (374,792 SNPs)
#   sample list     : samples passing all QC (contamination, relatedness, sex)
#
# Output: PCA eigenvectors + eigenvalues, scree data, and the % variance.
################################################################################

library(SNPRelate)
library(dplyr)

out_dir  <- "~/data/pca_results/PCA"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

gds_file <- "~/data/processed/DivCo_AMP1.0_GRCh38_1019Samples_filtered.gds"

## ---- 1. Load the pruned SNP set and the QC keep-list -----------------------
pruned <- readRDS("~/data/pca_results/QC/pruned_snp_ids.rds")
cat("Pruned SNPs:", length(pruned), "\n")

# keep-list: all samples NOT flagged for drop in the master list
flags <- read.csv("~/data/pca_results/QC/cohort_keep_drop_flags.csv",
                  stringsAsFactors = FALSE)
keep_samples <- flags %>% filter(!drop) %>% pull(sample.id) %>% as.character()
cat("Samples passing QC (PCA set):", length(keep_samples),
    "of", nrow(flags), "\n")

## ---- 2. Open GDS and run PCA -----------------------------------------------
genofile <- snpgdsOpen(gds_file)

# sanity: make sure keep_samples are actually in the GDS
gds_samples <- read.gdsn(index.gdsn(genofile, "sample.id"))
keep_samples <- intersect(keep_samples, gds_samples)
cat("Keep-list samples found in GDS:", length(keep_samples), "\n")

pca <- snpgdsPCA(
  genofile,
  sample.id   = keep_samples,   # QC-passed samples only
  snp.id      = pruned,         # LD-pruned, long-range-LD-excluded SNPs
  num.thread  = 4,
  eigen.cnt   = 32,             # compute top 32 PCs (enough for structure + scree)
  autosome.only = TRUE,
  verbose     = TRUE
)

## ---- 3. Variance explained + assemble the PC table -------------------------
# pct variance per PC
var_pct <- pca$varprop * 100
cat("\nVariance explained by top 10 PCs (%):\n")
print(round(var_pct[1:10], 3))

# eigenvectors -> data frame keyed by sample.id
pc_df <- as.data.frame(pca$eigenvect)
n_pc  <- ncol(pc_df)
colnames(pc_df) <- paste0("PC", seq_len(n_pc))
pc_df <- cbind(sample.id = pca$sample.id, pc_df)

## ---- 4. Save outputs -------------------------------------------------------
write.csv(pc_df, file.path(out_dir, "pca_eigenvectors.csv"), row.names = FALSE)

eig_df <- data.frame(PC = seq_along(pca$eigenval),
                     eigenvalue = pca$eigenval,
                     var_pct = pca$varprop * 100)
write.csv(eig_df, file.path(out_dir, "pca_eigenvalues.csv"), row.names = FALSE)

saveRDS(pca, file.path(out_dir, "pca_result.rds"))   # full object for later use

snpgdsClose(genofile)

cat("\nPCA complete. Outputs written to", out_dir, "\n")
cat("  pca_eigenvectors.csv : per-sample PC scores\n")
cat("  pca_eigenvalues.csv  : eigenvalues + variance explained (for scree)\n")
cat("  pca_result.rds       : full snpgdsPCA object\n")                            