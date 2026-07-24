################################################################################
# project_dropped_samples.R
# Project ALL QC-dropped samples onto the clean 974-sample PC space using the
# SNP loadings (projection weights). The clean PCA defines the axes; dropped
# samples are placed into that fixed space without altering it.
#
# NOTE: snpgdsPCASampLoading uses simple loading-based projection, which is
# shrinkage-biased (Prive 2020) -- projected points are pulled toward the
# origin, more so for higher PCs. Fine for qualitative "where do they land"
# diagnostics; the bias is conservative (makes dropped samples look MORE
# normal, so off-cluster landings are real).
#
# Inputs : pc_snp_loadings.rds (from estimate_snp_weights.R)
#          clean PCA result, filtered GDS, cohort_keep_drop_flags.csv
# Outputs: projected_dropped_samples.csv  (dropped samples' PC1-PC4 + reason)
#          overlay plots (clean space + projected dropped samples)
################################################################################

library(SNPRelate)
library(ggplot2)

pca_dir  <- "~/data/synapse_data/PCA"
qc_dir   <- "~/data/pca_results/QC"
fig_dir  <- "/home/ec2-user/AMP-AD_genetic_PCs/Results/PCA"
gds_file <- "~/data/processed/DivCo_AMP1.0_GRCh38_1019Samples_filtered.gds"
n_pc     <- 32
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

## ---- 1. Load clean PCA, loadings, and keep/drop flags ----------------------
pca      <- readRDS(file.path(pca_dir, "pca_result.rds"))
snp_load <- readRDS(file.path(pca_dir, "pc_snp_loadings.rds"))

flags <- read.csv(file.path(qc_dir, "cohort_keep_drop_flags.csv"),
                  stringsAsFactors = FALSE)
flags$sample.id <- as.character(flags$sample.id)
dropped_ids <- flags$sample.id[flags$drop %in% TRUE]
cat("Dropped samples to project:", length(dropped_ids), "\n")

## ---- 2. Project the dropped samples onto the clean space -------------------
# snpgdsPCASampLoading projects samples (by sample.id) onto the PC space
# defined by snp_load. The dropped samples are in the SAME GDS, so their
# genotypes are already aligned to the loading SNPs / allele coding.
genofile <- snpgdsOpen(gds_file)

samp_proj <- snpgdsPCASampLoading(snp_load, genofile,
                                  sample.id = dropped_ids,
                                  num.thread = 4)
# samp_proj$eigenvect : [dropped samples x PCs]
proj <- data.frame(sample.id = samp_proj$sample.id,
                   samp_proj$eigenvect[, seq_len(n_pc)],
                   stringsAsFactors = FALSE)
colnames(proj)[-1] <- paste0("PC", seq_len(n_pc))
proj$set <- "dropped (projected)"

snpgdsClose(genofile)

## ---- 3. Attach drop reason (for colouring) ---------------------------------
proj <- merge(proj,
              flags[, c("sample.id",
                        intersect("reason", colnames(flags)))],
              by = "sample.id", all.x = TRUE)

## ---- 4. Clean samples (the reference space) --------------------------------
clean <- data.frame(sample.id = pca$sample.id,
                    pca$eigenvect[, seq_len(n_pc)],
                    stringsAsFactors = FALSE)
colnames(clean)[-1] <- paste0("PC", seq_len(n_pc))
clean$set <- "clean (PCA)"
if ("reason" %in% colnames(proj)) clean$reason <- NA

## ---- 5. Save projected coordinates -----------------------------------------
write.csv(proj, file.path(pca_dir, "projected_dropped_samples.csv"),
          row.names = FALSE)

## ---- 6. Overlay plots ------------------------------------------------------
plot_df <- rbind(clean[, c("sample.id","PC1","PC2","PC3","PC4","set")],
                 proj[,  c("sample.id","PC1","PC2","PC3","PC4","set")])

mk <- function(xx, yy) {
  ggplot() +
    geom_point(data = subset(plot_df, set == "clean (PCA)"),
               aes_string(xx, yy), colour = "grey75", size = 1.4, alpha = 0.6) +
    geom_point(data = subset(plot_df, set == "dropped (projected)"),
               aes_string(xx, yy), colour = "red", size = 2.2, shape = 17) +
    labs(title = paste0(xx, " vs ", yy,
                        ": dropped samples projected onto clean PC space"),
         subtitle = "grey = clean 974 (defines axes); red = QC-dropped, projected",
         x = xx, y = yy) +
    theme_bw(base_size = 12)
}
ggsave(file.path(fig_dir, "projected_dropped_PC1_PC2.png"),
       mk("PC1","PC2"), width = 9, height = 7, dpi = 150)
ggsave(file.path(fig_dir, "projected_dropped_PC3_PC4.png"),
       mk("PC3","PC4"), width = 9, height = 7, dpi = 150)

cat("\nProjected", nrow(proj), "dropped samples.\n")
cat("Coords : ", file.path(pca_dir, "projected_dropped_samples.csv"), "\n")
cat("Plots  : ", fig_dir, "/projected_dropped_PC*.png\n", sep = "")
cat("Shrinkage bias pulls projected points toward origin (conservative).\n")