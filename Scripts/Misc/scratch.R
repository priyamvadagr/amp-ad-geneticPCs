x <- rosmap_bio_df[rosmap_bio_df$individualID %in% rosmap_divco_ampad1_common$specimenID_WGS_AMPAD1,]

#Samples with rosmap individual ids in WGS 
nn <- sample_list[sample_list %in% rosmap_bio_df$individualID]
#Of these, how many come from the link list?
x <- nn[nn %in% rosmap_overlap] #12 with identical biospecimen and individualIDs
y <- nn[nn %in% rosmap_divco_ampad1_no_overlap$individualID] #170 with individualIDs not equal to specimenID from linked list 
z <- nn[nn %in% AMP_AD_DivCo_wgs_df$specimenID] #175 
ll <- AMP_AD_DivCo_biospec_df[AMP_AD_DivCo_biospec_df$specimenID %in% z, ]
intersect(x, z)

x <- rosmap_bio_df[rosmap_bio_df$individualID %in% c(rosmap_divco_ampad1_no_overlap$individualID, rosmap_divco_ampad1_common$specimenID_WGS_AMPAD1), ]

x <- rosmap_bio_df[rosmap_bio_df$specimenID %in% c(rosmap_divco_ampad1_idlink$specimenID_WGS_AMPAD1), ] 
intersect(y, z)
sum(length(x), length(y), length(z)) == length(nn)

#Samples with Mayo data in the wgs 
mm <- sample_list[sample_list %in% mayo_bio_df$specimenID]
oo <- sample_list[sample_list %in% mayo_bio_df$individualID]
oo == mayo_divco_ampad1_idlink$specimenID


nn[nn %in% rosmap_overlap]

xx <- rosmap_bio_df[rosmap_bio_df$specimenID %in% nn[nn %in% rosmap_overlap],]
zz <- rosmap_wgs_df[rosmap_wgs_df$specimenID %in% rosmap_divco_ampad1_no_overlap$specimenID_WGS_AMPAD1, ]
length(sample_list[sample_list %in% rosmap_divco_ampad1_no_overlap$individualID])
y1 <- rosmap_bio_df[rosmap_bio_df$specimenID %in% rosmap_divco_ampad1_no_overlap$specimenID_WGS_AMPAD1,]
y2 <- rosmap_bio_df[rosmap_bio_df$individualID %in% rosmap_divco_ampad1_no_overlap$individualID,]



amp_divco <- AMP_AD_DivCo_wgs_df[!(AMP_AD_DivCo_wgs_df$specimenID %in% sample_list), ]
intersect(amp_divco$specimenID, rosmap_bio_df$individualID)
intersect(amp_divco$specimenID, mayo_wgs_df$individualID)


v <- rosmap_bio_df[rosmap_bio_df$specimenID %in% rosmap_individual_IDs,]
unique(v$assay)
dim(v[v$assay == "wholeGenomeSeq",]) #so 12 with identical individual IDs and specimen IDs 

N <- rosmap_bio_df[rosmap_bio_df$specimenID %in% z, ]
dim(N[N$assay == "wholeGenomeSeq",])

dim(rosmap_bio_df[rosmap_bio_df$specimenID %in% rosmap_divco_ampad1_idlink$specimenID_WGS_AMPAD1 & rosmap_bio_df$assay == 'wholeGenomeSeq', ])



amp_ad_wgs_df <- read.table('~/data/processed//metadata/DivCO/DivCo_assay_WGS_metadata.csv',
sep = ',', header = T)
amp_ad_wgs_df[amp_ad_wgs_df$specimenID %in% het_df[het_df$flag_het, 'sample.id'],]



DivCO_old <- read.table(
  "~/data/metadata/AMP-AD_DivCo/AMP-AD_DiverseCohorts_individual_metadata.csv",
  sep = ",", header = TRUE
)

DivCO_new <- read.table(
  "~/data/metadata/harmonized_metadata/AMP-AD_DiverseCohorts_individual_metadata_harmonized.csv",
  sep = ",", header = TRUE,
  quote = "",          # disable quote parsing
  fill = TRUE,         # fill missing fields
  comment.char = ""    # ignore # as comment character
)

dim(DivCO_new)
dim(DivCO_old)
setdiff(unique(DivCO_old[DivCO_old$individualID %in% AMP_AD_DivCo_biospec_df$individualID, 'individualID']),
 unique(DivCO_new[DivCO_new$individualID %in% AMP_AD_DivCo_biospec_df$individualID, 'individualID'])) 


MayoRNA_seq_new <- read.table(
  "~/data/metadata/harmonized_metadata/MayoRNAseq_individual_metadata_harmonized.csv",
  sep = ",", header = TRUE,
  quote = "",          # disable quote parsing
  fill = TRUE,         # fill missing fields
  comment.char = ""    # ignore # as comment character
)

MayoRNA_seq_old <- read.table(
  "~/data/metadata/AMP-AD_1.0/MayoRNAseq_individual_metadata.csv",
  sep = ",", header = TRUE
)

setdiff(unique(MayoRNA_seq_new[MayoRNA_seq_new$individualID %in% mayo_overlap, 'individualID']), 
unique(MayoRNA_seq_old[MayoRNA_seq_old$individualID %in% mayo_overlap, 'individualID']))


Rosmap_new <- read.table(
  "~/data/metadata/harmonized_metadata/ROSMAP_clinical_harmonized.csv",
  sep = ",", header = TRUE,
  quote = "",          # disable quote parsing
  fill = TRUE,         # fill missing fields
  comment.char = ""    # ignore # as comment character
)
Rosmap_old <- rosmap_individual_df <- read.table(
  "~/data/metadata/AMP-AD_1.0/ROSMAP_clinical.csv",
  sep = ",", header = TRUE
)
setdiff(unique(Rosmap_new[Rosmap_new$individualID %in% rosmap_divco_ampad1_idlink, 'individualID']),
unique(Rosmap_old[Rosmap_old$individualID %in% rosmap_divco_ampad1_idlink, 'individualID']))


setdiff(unique(Rosmap_new[Rosmap_new$individualID %in% rosmap_divco_ampad1_idlink$individualID, 'individualID']),
unique(Rosmap_old[Rosmap_old$individualID %in% rosmap_divco_ampad1_idlink$individualID, 'individualID']))


amp_ad_ind_df <- amp_ad_df[amp_ad_df <- ]