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
