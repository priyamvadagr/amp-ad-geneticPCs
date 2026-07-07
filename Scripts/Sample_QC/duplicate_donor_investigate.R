#Investigate duplicate donor
#Duplicate donor IDs are R3057101 and R4323048
#Both are from DivCO 
library(dplyr)
DivCO_metadata_dir <- '~/data/processed/metadata/DivCO'
amp_ad_wgs_df <- read.table(paste0(DivCO_metadata_dir, '/', 'DivCo_assay_WGS_metadata.csv'), sep = ',', header = T)
amp_ad_biospec_df <- read.table(paste0(DivCO_metadata_dir, '/', 'DivCo_biospecimen_metadata.csv'), sep = ',', header = T)
amp_ad_ind_df <- read.table(paste0(DivCO_metadata_dir, '/', 'DivCo_individual_metadata.csv'), sep = ',', header = T)

dup_pair_ids <- c('R3057101', 'R4323048')

#Subset to duplicate IDs
dup_wgs_df <- amp_ad_wgs_df[amp_ad_wgs_df$specimenID %in% dup_pair_ids, ]
dup_biospec_df <- amp_ad_biospec_df[amp_ad_biospec_df$specimenID %in% dup_pair_ids, ]
dup_ind_df <- amp_ad_ind_df[amp_ad_ind_df$individualID %in% amp_ad_biospec_df[amp_ad_biospec_df$specimenID %in% dup_pair_ids, 'individualID'],]

write.table(dup_ind_df, file = '~/data/pca_results/QC/duplicate_individual_metadata.txt', sep = '\t', quote = F)
write.table(dup_wgs_df, file = '~/data/pca_results/QC/duplicate_wgs_metadata.txt', sep = '\t', quote = F)

amp_ad_ind_df[amp_ad_ind_df$individualID == 'R2262117',]
