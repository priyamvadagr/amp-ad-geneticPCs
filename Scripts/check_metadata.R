library(dplyr)
sample_list <- read.table('~/data/processed/sample_ids.txt')
sample_list <- sample_list$V1
unprocessed_sample_list <- read.table('~/data/processed/sample_ids_unprocessed_chr9.txt')
unprocessed_sample_list <- unprocessed_sample_list$V1

AMP_AD_DivCo_df <- read.table('~/data/metadata/AMP-AD_DivCo/AMP-AD_DiverseCohorts_individual_metadata.csv', 
                                sep = ',', header = T)
length(sample_list[which(sample_list %in% AMP_AD_DivCo_df$individualID)])
AMP_AD_DivCo_wgs_df <- read.table('~/data/metadata/AMP-AD_DivCo/AMP-AD_DiverseCohorts_assay_WGS_metadata.csv', 
                                sep = ',', header = T)
divco_overlap <- sample_list[which(sample_list %in% AMP_AD_DivCo_wgs_df$specimenID)] #only 743 samples present of the 762 from the DivCo
length(unique(AMP_AD_DivCo_wgs_df$specimenID))
AMP_AD_DivCo_biospec_df <- read.table('~/data/metadata/AMP-AD_DivCo/AMP-AD_DiverseCohorts_biospecimen_metadata.csv', 
                                sep = ',', header = T)
length(unique(AMP_AD_DivCo_biospec_df$individualID))
length(sample_list[which(sample_list %in% AMP_AD_DivCo_biospec_df$specimenID)]) #this is also 743 

#AMP-AD 1.0 
rosmap_df <- read.table('~/data/metadata/AMP-AD_1.0/ROSMAP_assay_wholeGenomeSeq_metadata.csv', 
                                sep = ',', header = T)
rosmap_overlap_n <- length(sample_list[which(sample_list %in% rosmap_df$specimenID)]) #12 overlap 
rosmap_overlap <- sample_list[which(sample_list %in% rosmap_df$specimenID)] #12 overlap 

msbb_wgs_df <- read.table('~/data/metadata/AMP-AD_1.0/MSBB_assay_wholeGenomeSeq_metadata.csv', 
                                sep = ',', header = T)
msbb_overlap <- sample_list[which(sample_list %in% as.character(msbb_wgs_df$specimenID))]
msbb_df <- read.table('~/data/metadata/AMP-AD_1.0/MSBB_individual_metadata.csv', 
                                sep = ',', header = T)
length(sample_list[which(sample_list %in% as.character(msbb_df$individualID))]) #no overlap with sample 
mayo_wgs_df <- read.table('~/data/metadata/AMP-AD_1.0/MayoRNAseq_assay_wholeGenomeSeq_metadata.csv',
                                sep = ',', header = T)
mayo_overlap <- sample_list[which(sample_list %in% as.character(mayo_wgs_df$specimenID))] #94 samples in wgs data 

mayo_df <- read.table('~/data/metadata/AMP-AD_1.0/MayoRNAseq_individual_metadata.csv',
                                sep = ',', header = T)
length(sample_list[which(sample_list %in% as.character(mayo_df$individualID))]) #specimenID and individualID same for mayodf

combined_overlap <- c(mayo_overlap, msbb_overlap, rosmap_overlap, divco_overlap)

length(setdiff(sample_list, combined_overlap))
divco_ampad1_idlink <- read.table('~/data/metadata/DivCo2_AMPAD1_IDlink2.csv', sep = ',', header = T)
unique(divco_ampad1_idlink$study_WGS_AMPAD1)
mayo_divco_ampad1_idlink <- divco_ampad1_idlink[divco_ampad1_idlink$study_WGS_AMPAD1 == "MAYO_AMPAD1",]
rosmap_divco_ampad1_idlink <- divco_ampad1_idlink[divco_ampad1_idlink$study_WGS_AMPAD1 == "ROSMAP_AMPAD1",]
#Filter to samples in list 
mayo_divco_ampad1_idlink <- mayo_divco_ampad1_idlink[mayo_divco_ampad1_idlink$individualID %in% sample_list,]
length(set_intersection(combined_overlap, mayo_divco_ampad1_idlink$specimenID)) #all Mayo IDs already present
rosmap_divco_ampad1_idlink <- rosmap_divco_ampad1_idlink[rosmap_divco_ampad1_idlink$individualID %in% sample_list,]
length(set_intersection(combined_overlap, rosmap_divco_ampad1_idlink$specimenID)) #only 12 rosmap ids overlap
rosmap_divco_ampad1_common <- rosmap_divco_ampad1_idlink[rosmap_divco_ampad1_idlink$specimenID %in% sample_list,]
setdiff(rosmap_divco_ampad1_common$specimenID, rosmap_divco_ampad1_common$individualID) #specimenID == idividualID for these 
#Need to map the individual IDs from DivCo (which is the specimen ID in the wgs DivCo) to AMP AD 1.0 specimenID 
rosmap_divco_ampad1_no_overlap <- rosmap_divco_ampad1_idlink[!(rosmap_divco_ampad1_idlink$specimenID %in% sample_list), ]
#Update the rosmap overlap 
combined_overlap <- c(combined_overlap, rosmap_divco_ampad1_no_overlap$individualID)
length(setdiff(sample_list, combined_overlap))
