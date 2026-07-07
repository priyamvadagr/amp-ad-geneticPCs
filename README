# AMP-AD / DivCo Genetic PCA Workflow

Population-structure PCA workflow for the DivCo + AMP-AD 1.0 joint-called WGS
dataset: variant filtering, sample-level QC (contamination, relatedness,
heterozygosity, sex check), LD pruning, PCA, and outlier detection.


## Data

### Genetic Data

syn67158725 - DivCO + AMP AD 1.0 joint called samples

### Metadata

```
AMP-AD 1.0    MSBB, Mayo, and ROSMAP metadata
  MSBB
    syn22360825 - MSBB_assay_wholeGenomeSeq_metadata.csv
    syn6101474  - MSBB_individual_metadata.csv
    syn21893059 - MSBB_biospecimen_metadata.csv
  Mayo
    syn23481994 - MayoRNAseq_assay_wholeGenomeSeq_metadata.csv
    syn73713766 - MayoRNAseq_individual_metadata_harmonized.csv
    syn20827192 - MayoRNAseq_biospecimen_metadata.csv
    syn21442783 - AMP-AD_Mayo_WGS_QualityControlSampleMetrics.csv
                  (% contamination metric; file itself lives in
                  AMP-AD_1.0_WG/, see below)
  ROSMAP
    syn21314542 - ROSMAP_assay_wholeGenomeSeq_metadata.csv
    syn73713768 - ROSMAP_clinical_harmonized.csv
    syn21323366 - ROSMAP_biospecimen_metadata.csv

AMP-AD Diverse Cohorts metadata
  syn51757644 - AMP-AD_DiverseCohorts_assay_WGS_metadata.csv
  syn51757645 - AMP-AD_DiverseCohorts_biospecimen_metadata.csv
  syn73713769 - AMP-AD_DiverseCohorts_individual_metadata_harmonized.csv
```

## Repository structure

```
Data/                     Small reference/annotation files used by the scripts
                          (not genotype or metadata; see Data/README)

Scripts/
  SNP_QC/                 VCF -> GDS conversion, LD pruning
  Sample_QC/              Sample missingness, heterozygosity, relatedness
                          (KING), sex check, and drop-list decisions
  PCA/                    Run PCA, scree/Tracy-Widom PC selection, PCA
                          outlier detection (LOF), plotting by ancestry
  Misc/                   One-off metadata reconciliation / investigation
                          scripts (DivCo <-> AMP-AD 1.0 <-> Synapse ID
                          mapping, sex-record creation, scratch analysis)
  check_metadata.R,
  check_wgs_metadata.r    Ad hoc checks reconciling sample lists across
                          metadata files

Results/
  Filtering/              QC diagnostic plots: heterozygosity, kinship,
                          sequence contamination, sex check
  Metadata_reconcile/     Rendered output of the metadata reconciliation
                          notebooks
  PCA/                    Scree plots, PCs colored by ancestry/cohort/
                          datacenter, PCA outlier plots

logs/                     Run logs from long-running steps (GDS conversion,
                          LD pruning)
renv.lock                 R package environment (renv)
.Rprofile                 renv activation
```

Raw genotype/metadata files and other large or controlled-access data are not
committed (see `.gitignore`) and are expected under `~/data/`.

## Workflow

### 1. Data download

```
synapse get -r syn2580853 --downloadLocation ~/data/raw/
```

### 2. Data processing

Concatenate all autosome files:

```
bcftools concat --threads 8 ~/data/raw/*.chr{1..22}.vcf.gz \
  -Oz -o ~/data/processed/DivCo_AMP1.0_GRCh38_1019Samples_JointCalls_Recalibrated_Annotated_05-25-2025.all_auto.vcf.gz
```

Filter variants using the following criteria:
- VQSR-PASS
- MAF > 0.01
- Variant missingness > 1%

```
bcftools view \
  -f PASS \
  -m2 -M2 -v snps \
  -q 0.01:minor \
  -e 'F_MISSING > 0.01' \
  --threads 8 \
  ~/data/processed/*.all_auto.vcf.gz \
  -Oz -o ~/data/processed/DivCo_AMP1.0_GRCh38_1019Samples_filtered.vcf.gz

tabix -p vcf ~/data/processed/DivCo_AMP1.0_GRCh38_1019Samples_filtered.vcf.gz
```

Convert the filtered VCF to GDS for SNPRelate (`Scripts/SNP_QC/convert_vcf_to_gds.R`).

### 3. Sex check

Convert to a bed file:

```
plink \
  --vcf ~/data/processed/DivCo_AMP1.0_GRCh38_1019Samples_chrX_filtered.vcf.gz \
  --split-x hg38 \
  --double-id \
  --make-bed \
  --out ~/data/processed/chrX_filt_check
```

Compare each sample's recorded sex to the sex inferred from X-chromosome
genotypes, flagging mismatches:

```
plink \
  --bfile ~/data/processed/chrX_filt_check \
  --update-sex ~/data/processed/recorded_sex.txt \
  --check-sex \
  --out ~/data/processed/sexcheck_filt_vs_recorded
```

`recorded_sex.txt` is built by `Scripts/Misc/create_sex_record_plink.R`; results
are plotted by `Scripts/Sample_QC/plot_Xsex_F_byAncestry.R`.

### 4. Sample QC (`Scripts/Sample_QC/`)

`qc_rel_heterozyg_missingness.r` computes, per sample: call rate
(missingness), heterozygosity outliers, and KING-robust relatedness
(flagging duplicates/relatives). It is run interactively; downstream scripts
consume its outputs:

- `plot_heterozygosity_by_race.r`, `Plot_het_outliers.R` — heterozygosity
  diagnostics by race/technical metrics
- `plot_relatedness.r` — kinship (IBS0 vs. kinship) diagnostics
- `decide_related_drops.R` — resolves which member of a related pair/triad
  to drop
- `duplicate_donor_investigate.R` — investigates specific duplicate donor IDs
- `Create_master_drop_list.R` — consolidates all exclusions (contamination,
  relatedness, sex mismatch) into one master drop list with reasons

### 5. SNP QC / LD pruning (`Scripts/SNP_QC/`)

`LD_pruning.R` LD-prunes the QC-passed genotypes and excludes long-range LD
regions (`Data/high-LD-regions-hg38-GRCh38.txt`).

### 6. PCA (`Scripts/PCA/`)

- `runPCA.R` — population-structure PCA on the QC-passed, LD-pruned dataset
- `Scree_plot.R` — scree plot + Tracy-Widom test to pick the number of
  informative PCs
- `TW_pca.R` — Tracy-Widom test via LEA (alternative implementation)
- `detect_PCA_outliers.R` — ancestry-aware PCA outlier detection (KNN-based
  probabilistic LOF, Prive et al. 2020)
- `Plot_PCs.R` — PCs colored by self-reported race/ethnicity, cohort, and
  data center to confirm they capture ancestry rather than batch

### Metadata reconciliation (`Scripts/Misc/`)

Ad hoc notebooks/scripts mapping sample IDs between DivCo, AMP-AD 1.0, and
Synapse metadata (`check_metadata.rmd`, `reconcile_NYG_QC_with_synapse.rmd`,
`Check_MayoRNAseq.R`). Rendered notebook output lives in
`Results/Metadata_reconcile/`.
