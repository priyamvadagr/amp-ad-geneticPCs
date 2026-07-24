################################################################################
# estimate_snp_weights.R
# Compute SNP loadings ("projection weights") from the clean 974-sample PCA,
# exported with variant identifiers + allele coding so new/held-out samples
# can be projected onto this exact PC space later.
#
# Loadings are computed on the SAME LD-pruned SNP set the PCA used.
# NOTE: simple loading-based projection (snpgdsPCASampLoading) is shrinkage-
# biased (Prive 2020); fine for qualitative placement, use OADP for unbiased.
#
# Inputs : clean PCA result (pca_result.rds), filtered GDS
# Outputs: pc_snp_loadings.rds            (SNPLoading object, for projecting)
#          pc_projection_weights.csv      (loadings + variant IDs/alleles)
################################################################################

library(SNPRelate)

pca_dir  <- "~/data/synapse_data/PCA"
gds_file <- "~/data/processed/DivCo_AMP1.0_GRCh38_1019Samples_filtered.gds"
n_pc     <- 32   # all PCs

## ---- 1. Load the clean PCA result ------------------------------------------
pca <- readRDS(file.path(pca_dir, "pca_result.rds"))
cat("PCA samples:", length(pca$sample.id),
    "| PCA SNPs:", length(pca$snp.id), "\n")

genofile <- snpgdsOpen(gds_file)

## ---- 2. Compute SNP loadings from the PCA ----------------------------------
# Uses the SNPs recorded in the PCA object (the LD-pruned set). Deterministic.
snp_load <- snpgdsPCASNPLoading(pca, genofile, num.thread = 4)
# snp_load$snploading : [eigenvector x SNP] matrix
# snp_load$snp.id     : SNP ids, aligned to columns of snploading
cat("Loadings computed for", length(snp_load$snp.id), "SNPs,",
    nrow(snp_load$snploading), "PCs\n")

saveRDS(snp_load, file.path(pca_dir, "pc_snp_loadings.rds"))

## ---- 3. Attach variant IDs + allele coding (ESSENTIAL for projection) ------
# Pull chr/pos/alleles/rsID from the GDS, keyed by snp.id, so downstream
# projection can align variants and allele orientation correctly.
all_snp_id  <- read.gdsn(index.gdsn(genofile, "snp.id"))
snp_chr     <- read.gdsn(index.gdsn(genofile, "snp.chromosome"))
snp_pos     <- read.gdsn(index.gdsn(genofile, "snp.position"))
snp_allele  <- read.gdsn(index.gdsn(genofile, "snp.allele"))   # "REF/ALT"
# rsID if present in the GDS
rs <- tryCatch(read.gdsn(index.gdsn(genofile, "snp.rs.id")),
               error = function(e) rep(NA_character_, length(all_snp_id)))

snp_info <- data.frame(snp.id = all_snp_id, chr = snp_chr, pos = snp_pos,
                       allele = snp_allele, rs.id = rs,
                       stringsAsFactors = FALSE)


## ---- 4. Assemble the projection-weights table ------------------------------
load_mat <- t(snp_load$snploading[seq_len(n_pc), , drop = FALSE])  # SNPs x PCs
colnames(load_mat) <- paste0("PC", seq_len(n_pc), "_loading")

weights <- data.frame(snp.id = snp_load$snp.id, load_mat,
                      stringsAsFactors = FALSE)
weights <- merge(snp_info, weights, by = "snp.id")   # add chr/pos/allele/rsID
# keep the PCA's SNP order
weights <- weights[match(snp_load$snp.id, weights$snp.id), ]

write.csv(weights, file.path(pca_dir, "pc_projection_weights.csv"),
          row.names = FALSE)

snpgdsClose(genofile)

cat("\nWrote:\n  ", file.path(pca_dir, "pc_snp_loadings.rds"),
    "  (SNPLoading object -> feed to snpgdsPCASampLoading to project)\n  ",
    file.path(pca_dir, "pc_projection_weights.csv"),
    "  (loadings + variant IDs/alleles)\n")
cat("\nColumns: snp.id, chr, pos, allele (REF/ALT), rs.id, PC1_loading..PC",
    n_pc, "_loading\n", sep = "")
cat("To project new samples: align them to this SNP set + allele coding,\n")
cat("then snpgdsPCASampLoading(snp_load, newGDS). Mind shrinkage bias.\n")