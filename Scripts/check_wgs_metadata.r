
colnames(qc_assay_map) #final specimenID is the name in wgs_qc
length(setdiff(qc_assay_map$final_specimenid, wgs_qc$Sample))
length(setdiff(qc_assay_map$filename.sample, wgs_qc$Sample)) #135 samples from the map file not present in wgs_qc file, these were the ones processed separately 
length(setdiff(wgs_qc$Sample, qc_assay_map$filename.sample)) #34 samples from qc file not present in mapping file 
length(setdiff(qc_assay_map$final_specimenid, wgs_assay_metadata$specimenID)) #293 samples from QC map in wgs_metadata 
y <- intersect((qc_assay_map$final_specimenid, wgs_assay_metadata$specimenID)) #293 samples from QC map in wgs_metadata 


# add 'WGS' to the specimenID in mapping file
x <- paste0(setdiff(qc_assay_map$final_specimenid, wgs_assay_metadata$specimenID), '_WGS')

setdiff(x, wgs_assay_metadata$specimenID) # all files from the 

length(setdiff(setdiff(qc_assay_map$final_specimenid, wgs_assay_metadata$specimenID), wgs_qc$Sample))

setdiff(qc_assay_map[qc_assay_map$final_specimenid %in% setdiff(qc_assay_map$final_specimenid, wgs_assay_metadata$specimenID), 'filename.sample'], wgs_qc$Sample)
# add a third column 


length(setdiff(wgs_assay_metadata$specimenID, qc_assay_map$final_specimenid)) #all samples from wgs qc are in the wgs assay metadata 



qc_assay_map[qc_assay_map$filename.sample %in% setdiff(qc_assay_map$filename.sample, wgs_qc$Sample),'final_specimenid']


dropped_samples <- read.table('/home/ec2-user/data/metadata/AMP-AD_DivCO-WGS/dropped_samples.tsv', sep = '\t', header = T)

length(intersect(dropped_samples$Sample.ID, setdiff(wgs_qc$Sample, qc_assay_map$filename.sample)))
setdiff(setdiff(wgs_qc$Sample, qc_assay_map$filename.sample), dropped_samples$Sample.ID)


#14 samples from wgs qc with no source
z <- wgs_qc[wgs_qc$Sample %in% setdiff(setdiff(wgs_qc$Sample, qc_assay_map$filename.sample), dropped_samples$Sample.ID), 'Sample']



wgs_assay_metadata[wgs_assay_metadata$Estimated.Library.Size ==  1535073070, ]


sample_ids <- read.table('/home/ec2-user/data/processed/sample_ids.txt')
sample_ids <- sample_ids$V1

setdiff(setdiff(setdiff(wgs_qc$Sample, qc_assay_map$filename.sample), dropped_samples$Sample.ID), sample_ids)

Mayo_biospec_df <- read.table('/home/ec2-user/data/metadata/AMP-AD_1.0/MayoRNAseq_biospecimen_metadata.csv', sep = ',', header = T)

intersect(z, as.character(Mayo_biospec_df$individualID))

head(Mayo_biospec_df$individualID)

rosmap_biospec_df <- read.table('/home/ec2-user/data/metadata/AMP-AD_1.0/ROSMAP_biospecimen_metadata.csv', sep = ',', header = T)

intersect(z, rosmap_biospec_df$individualID)
intersect(z, rosmap_biospec_df$specimenID)

intersect(z, qc_assay_map$final_specimenid)

intersect(sample_ids, z)

hist(wgs_assay_df$Mean.Coverage)
summary(wgs_assay_df$Mean.Coverage)
