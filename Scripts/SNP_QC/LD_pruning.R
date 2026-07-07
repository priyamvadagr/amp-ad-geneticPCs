## LD pruning + long-range-LD exclusion
# Uses the GRCh38-native high-LD-regions list from plinkQC
# (derived from Anderson et al. 2010, lifted hg18 -> hg38)


## ---- 0. Setup --------------------------------------------------------------
library(SNPRelate)
library(dplyr)

gds_file <- "~/data/processed/DivCo_AMP1.0_GRCh38_1019Samples_filtered.gds"
out_dir  <- "~/data/pca_results/QC"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
set.seed(1)   # LD pruning is stochastic; seed for reproducibility
genofile <- snpgdsOpen(gds_file)

# sanity check: dimensions
snpgdsSummary(gds_file)
all_samples <- read.gdsn(index.gdsn(genofile, "sample.id"))
length(all_samples)
 
# load and process LD region file 
lrld_file <- "~/AMP-AD_genetic_PCs/Data/high-LD-regions-hg38-GRCh38.txt"
lrld_raw <- read.table(lrld_file, header = FALSE, stringsAsFactors = FALSE,
                       fill = TRUE)

# keep first three columns as chr/start/end (+ optional 4th as id)
lrld <- data.frame(
  chr   = as.character(lrld_raw[[1]]),
  start = as.numeric(lrld_raw[[2]]),
  end   = as.numeric(lrld_raw[[3]]),
  id    = as.character(lrld_raw[[4]]),
  stringsAsFactors = FALSE
)

cat("Long-range-LD regions loaded:", nrow(lrld), "\n")
print(head(lrld))

# ---- Map SNPs -> chr/pos and drop those inside LRLD regions ------------------
snp_chr <- read.gdsn(index.gdsn(genofile, "snp.chromosome"))
snp_pos <- read.gdsn(index.gdsn(genofile, "snp.position"))
snp_id  <- read.gdsn(index.gdsn(genofile, "snp.id"))

# normalize chromosome naming to match the lrld$chr style ("chr6")
snp_chr_chr <- ifelse(grepl("^chr", snp_chr), snp_chr, paste0("chr", snp_chr))

in_lrld <- rep(FALSE, length(snp_id))
for (i in seq_len(nrow(lrld))) {
  in_lrld <- in_lrld | (snp_chr_chr == lrld$chr[i] &
                        snp_pos    >= lrld$start[i] &
                        snp_pos    <= lrld$end[i])
}
keep_snp <- snp_id[!in_lrld]
cat("Long-range-LD regions used:", nrow(lrld), "\n")
cat("SNPs excluded in long-range-LD regions:", sum(in_lrld), "\n")
cat("SNPs remaining before pruning:", length(keep_snp), "\n")

# ---- LD pruning on the remaining SNPs ---------------------------------------
snpset <- snpgdsLDpruning(
  genofile,
  snp.id        = keep_snp,
  method        = "corr",
  ld.threshold  = 0.2,
  maf           = 0.01,
  missing.rate  = 0.01,
  slide.max.bp  = 500000,
  autosome.only = TRUE
)
pruned <- unlist(snpset, use.names = FALSE)
cat("SNPs after LD pruning:", length(pruned), "\n")
saveRDS(pruned, file.path(out_dir, "pruned_snp_ids.rds"))