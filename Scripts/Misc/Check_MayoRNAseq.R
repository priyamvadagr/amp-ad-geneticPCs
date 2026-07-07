# Check MayoRNAseq data from syn21442783
library(dplyr)
MayoRNAseq_wgs_df <- read.table('/home/ec2-user/data/metadata/AMP-AD_1.0_WG/AMP-AD_Mayo_WGS_QualityControlSampleMetrics.csv', sep = ',', header = T)
WGS_processed <- read.table('/home/ec2-user/data/processed/metadata/Mayo/Mayo_assay_WGS_metadata.csv', sep = ',', header = T)
sample_ids <- WGS_processed$specimenID
sample_ids_new_dat <- paste0('s_', sample_ids)
length(intersect(sample_ids_new_dat, MayoRNAseq_wgs_df$sampleID))
MayoRNAseq_wgs_df <- MayoRNAseq_wgs_df[MayoRNAseq_wgs_df$sampleID %in% sample_ids_new_dat,]
summary(MayoRNAseq_wgs_df$PercentDuplicateReads)
