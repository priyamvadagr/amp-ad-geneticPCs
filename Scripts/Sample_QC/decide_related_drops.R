################################################################################
# decide_related_drops_v3.R
#
# Workflow:
#   1. Genuine related pairs (2nd-degree+, excl. duplicate)
#   2. DROP contaminated (>2%) samples first
#   3. Plot Tier-1 QC for the remaining related-pair samples (no pair colouring)
#   4. Build the relatedness graph; separate TRIADS (clusters >=3) from
#      DISJOINT PAIRS (isolated 2-sample components)
#   5. Resolve disjoint pairs by specimen count (keep the one with more)
#   6. Triads handled separately (graph removal / manual)
#
# Requires: kin, het_df, samp_miss_df + DivCo metadata files
################################################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(igraph)
library(ggraph)

out_dir <- "~/data/pca_results/QC"
fig_dir <- "/home/ec2-user/AMP-AD_genetic_PCs/Results/Filtering/Kinship"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

kin          <- read.table('~/data/pca_results/QC/Sample_relatedness_flag.csv',
                           sep = ',', header = TRUE)
het_df       <- read.table('~/data/pca_results/QC/Heterozygosity_flag.txt')
samp_miss_df <- read.table('~/data/pca_results/QC/Sample_missingness.txt')
samples_to_drop_sex <- read.table('~/data/pca_results/QC/suggested_drop_sex_mismatch.txt')

all_sample_ids <- het_df$sample.id

## ---- 1. Genuine related pairs (2nd-degree+, excl. duplicate) ---------------
to_resolve <- kin %>%
  filter(kinship > 0.0884, relationship != 'dup/MZ') %>%
  arrange(desc(kinship)) %>%
  select(ID1, ID2, IBS0, kinship, relationship = matches("relationship"))
cat("Pairs to resolve (2nd-degree+, excl. duplicate):", nrow(to_resolve), "\n")

#check if any of these samples have been falgged for execcessive heterogeneity or sex mismatch 
het_flag_sample_ids <- het_df[het_df$flag_het == TRUE, 'sample.id']
samples_to_drop_sex  <- samples_to_drop_sex$sample.id
print(intersect(unique(c(to_resolve$ID1, to_resolve$ID2)), c(het_flag_sample_ids, samples_to_drop_sex)))

## ---- 2. Per-sample QC + DROP contaminated (>2%)  ----------------------
amp_ad_wgs_df <- read.table('~/data/processed/metadata/DivCO/DivCo_assay_WGS_metadata.csv',
                            sep = ',', header = TRUE)
amp_ad_wgs_df$GQN <- suppressWarnings(as.numeric(amp_ad_wgs_df$GQN))

seq_metrics <- c("Mean.Coverage", "Percent.Sequence.Contamination",
                 "Percent.PF.Aligned", "Percent.Total.Dup.", "Strand.Balance",
                 "AT.Dropout", "GC.Dropout", "Median.Insert.Size", "GQN")

sample_qc <- samp_miss_df %>%
  select(sample.id, miss.rate) %>%
  left_join(het_df %>% select(sample.id, inbreeding, flag_het), by = "sample.id") %>%
  left_join(amp_ad_wgs_df %>% select(specimenID, any_of(seq_metrics)),
            by = c("sample.id" = "specimenID")) %>%
  mutate(seq_contaminated = !is.na(Percent.Sequence.Contamination) &
                            Percent.Sequence.Contamination > 2)

# samples in related pairs
related_ids <- unique(c(to_resolve$ID1, to_resolve$ID2))

# contaminated samples among the related-pair members -> dropped first
contam_related <- sample_qc %>%
  filter(sample.id %in% related_ids, seq_contaminated) %>% pull(sample.id)
cat("Contaminated (>2%) samples among related pairs (dropped first):",
    length(contam_related), "\n"); print(contam_related)

# remaining related-pair members after contamination drop
related_remaining <- setdiff(related_ids, contam_related)

# pairs that still have BOTH members remaining (still need resolving);
# pairs where a member was dropped for contamination are already resolved
to_resolve <- to_resolve %>%
  mutate(both_remaining = ID1 %in% related_remaining & ID2 %in% related_remaining)
cat("Pairs already resolved by contamination drop:",
    sum(!to_resolve$both_remaining), "\n")
pairs_open <- to_resolve %>% filter(both_remaining)

## ---- 3. PLOT QC for remaining related-pair samples (no pair colour) -
plot_metrics  <- c("miss.rate", "Percent.Sequence.Contamination", "inbreeding")
metric_labels <- c(miss.rate = "Sample missingness",
                   Percent.Sequence.Contamination = "% Sequence contamination",
                   inbreeding = "Inbreeding coefficient (F)")

cohort_long <- sample_qc %>%
  select(sample.id, any_of(plot_metrics)) %>%
  pivot_longer(-sample.id, names_to = "metric", values_to = "value") %>%
  filter(!is.na(value))

# remaining related-pair samples (after contamination drop), de-duplicated
related_long <- sample_qc %>%
  filter(sample.id %in% related_remaining) %>%
  select(sample.id, any_of(plot_metrics)) %>%
  pivot_longer(-sample.id, names_to = "metric", values_to = "value") %>%
  filter(!is.na(value))

p_dist <- ggplot(cohort_long, aes(x = metric, y = value)) +
  geom_violin(fill = "grey92", colour = "grey65", scale = "width") +
  geom_jitter(width = 0.15, size = 0.5, alpha = 0.18, colour = "grey65") +
  # related-pair samples in a single highlight colour (no per-pair colouring)
  geom_jitter(data = related_long, colour = "#B2182B", size = 2.4,
              width = 0.12, alpha = 0.85) +
  facet_wrap(~ metric, scales = "free", ncol = 3,
             labeller = labeller(metric = metric_labels)) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  labs(title = "Remaining related-pair samples vs. cohort",
       subtitle = "After dropping contaminated (>2%) samples from related pairs. Red = related-pair members",
       x = NULL, y = "value")

ggsave(file.path(fig_dir, "related_pairs_QC_violin_tier1.png"),
       p_dist, width = 12, height = 5, dpi = 150)

## ---- 4. Build graph; separate TRIADS (>=3) from DISJOINT PAIRS -------------
g <- graph_from_data_frame(
  pairs_open %>% select(ID1, ID2, kinship, relationship), directed = FALSE)

comp <- components(g)
V(g)$cluster      <- comp$membership
V(g)$cluster_size <- comp$csize[comp$membership]

# disjoint pairs = components of exactly 2 samples
disjoint_clusters <- which(comp$csize == 2)
triad_clusters    <- which(comp$csize >= 3)

disjoint_ids <- V(g)$name[comp$membership %in% disjoint_clusters]
triad_ids    <- V(g)$name[comp$membership %in% triad_clusters]

cat("\nDisjoint pairs (2-sample components):", length(disjoint_clusters), "\n")
cat("Triad+ clusters (>=3 samples):", length(triad_clusters), "\n")

## ---- 5. Resolve DISJOINT PAIRS by specimen count ---------------------------
amp_ad_biospec_df <- read.table(
  '~/data/metadata/AMP-AD_DivCo/AMP-AD_DiverseCohorts_biospecimen_metadata.csv',
  sep = ',', header = TRUE)

specimen_counts <- amp_ad_biospec_df %>%
  group_by(individualID) %>%
  summarise(
    n_specimens = n(),
    assay_types = paste(sort(unique(assay)), collapse = ", "),
    .groups = "drop"
  )
# map the disjoint-pair samples -> individualID -> specimen count, with cluster id

disjoint_df <- data.frame(sample.id = disjoint_ids,
                          cluster = comp$membership[match(disjoint_ids, V(g)$name)]) %>%
  left_join(amp_ad_biospec_df %>% select(specimenID, individualID) %>% distinct(),
            by = c("sample.id" = "specimenID")) %>%
    left_join(sample_qc %>% select(sample.id, miss.rate, inbreeding), by = "sample.id") %>%
  left_join(amp_ad_ind_df %>% select(sex, specimenID), by = c("sample.id" = "specimenID")) %>%          
  left_join(specimen_counts, by = "individualID")

# per disjoint pair (cluster), keep the member with MORE specimens; flag ties
disjoint_decision <- disjoint_df %>%
  group_by(cluster) %>%
  arrange(cluster,
          desc(n_specimens),                       # more specimens first
          desc(tolower(sex) == "male"),            # if tied, male first
          miss.rate,                          # then lower missingness
          abs(inbreeding))  %>%                    # then F closer to 0)         
  mutate(keep = row_number() == 1,
         decided_by_specimens = n_distinct(n_specimens) > 1) %>%
  ungroup() %>%
  arrange(cluster, desc(n_specimens))

cat("\nDisjoint-pair resolution by specimen count",
    "(keep = TRUE -> retain; tie -> pick male -> lexicographic sorting):\n")
print(disjoint_decision %>%
        select(cluster, sample.id, individualID, n_specimens, keep, tie), n = 100)
write.csv(disjoint_decision,
          file.path(out_dir, "disjoint_pairs_specimen_decision.csv"), row.names = FALSE)

# suggested drops from disjoint pairs (the non-kept member; ties resolved manually)
disjoint_drops <- disjoint_decision %>%
  filter(!keep & !tie) %>% pull(sample.id)
disjoint_ties  <- disjoint_decision %>% filter(tie) %>% distinct(cluster) %>% pull(cluster)
cat("\nDisjoint-pair drops (specimen-based):\n"); print(disjoint_drops)
cat("Disjoint pairs tied on specimens:",
    paste(disjoint_ties, collapse = ", "), "\n")

## ---- 6. TRIADS: visualise and handle separately ----------------------------
amp_ad_ind_df <- read.table(
  '~/data/processed/metadata/DivCO/DivCo_individual_metadata.csv',
  sep = ',', header = TRUE)

  g_triads <- induced_subgraph(g, V(g)[comp$membership %in% triad_clusters])

  nid <- data.frame(specimenID = V(g_triads)$name, stringsAsFactors = FALSE)
  nid$individualID <- amp_ad_biospec_df$individualID[
    match(nid$specimenID, amp_ad_biospec_df$specimenID)]
  nid$sex  <- amp_ad_ind_df$sex[ match(nid$individualID, amp_ad_ind_df$individualID)]
  nid$race <- amp_ad_ind_df$race[match(nid$individualID, amp_ad_ind_df$individualID)]
  nid$n_specimens <- specimen_counts$n_specimens[
    match(nid$individualID, specimen_counts$individualID)]
  V(g_triads)$sex         <- nid$sex
  V(g_triads)$race        <- nid$race
  V(g_triads)$n_specimens <- nid$n_specimens

  triad_plot <- ggraph(g_triads, layout = "fr") +
    geom_edge_link(aes(label = relationship), colour = "grey60",
                   angle_calc = "along", label_dodge = unit(2.5, "mm"),
                   label_size = 5) +
    geom_node_point(aes(colour = sex), size = 14) +
    geom_node_text(aes(label = paste0(name, "\n(n=", n_specimens, ")")),
                   size = 4, vjust = -1.4) +
    geom_node_text(aes(label = race), size = 4, vjust = 2.6,
                   fontface = "italic", colour = "grey30") +
    scale_colour_manual(values = c("male" = "#4981BF", "female" = "#B2182B"),
                        na.value = "grey70", name = "Sex") +
    scale_x_continuous(expand = expansion(mult = 0.2)) +
    scale_y_continuous(expand = expansion(mult = 0.1)) +
    facet_nodes(~cluster, scales = "free", ncol = 3) +
    coord_cartesian(clip = "off") +
    labs(title = "Kinship triads (clusters >= 3 samples)",
         subtitle = "Node colour = sex; n = specimen count; italic = race. Resolve by removing the high-degree / lower-specimen node.") +
    theme_graph(base_size = 12) +  theme(strip.text = element_text(margin = margin(b = 2, t = 2)),
        panel.spacing = unit(1, "lines"))


  ggsave(file.path(fig_dir, "related_pairs_triads.png"),
         triad_plot, width = 10, height = 10, dpi = 150)
  cat("\nTriad plot written; resolve triads by removing the central/high-degree\n")
  cat("node (or the lower-specimen member) so no edge remains.\n")

#resolve triads via specimen count
triad_df <- data.frame(sample.id = triad_ids,
                          cluster = comp$membership[match(triad_ids, V(g)$name)]) %>%
  left_join(amp_ad_biospec_df %>% select(specimenID, individualID) %>% distinct(),
            by = c("sample.id" = "specimenID")) %>%
  left_join(amp_ad_ind_df %>% select(sex, specimenID), by = c("sample.id" = "specimenID")) %>% 
  left_join(sample_qc %>% select(sample.id, miss.rate, inbreeding), by = "sample.id") %>%         
  left_join(specimen_counts, by = "individualID")
triad_decision <- triad_df %>%      # triad_df = triad samples with sex, n_specimens, QC
  group_by(cluster) %>%
  arrange(cluster,
          desc(n_specimens),                 # most specimens
          desc(tolower(sex) == "male"),       # then male
          miss.rate,                          # then lower missingness
          abs(inbreeding),                    # then F closer to 0
          sample.id) %>%                       # then lexicographic
    mutate(keep = row_number() == 1,
         decided_by_specimens = n_distinct(n_specimens) > 1) %>%        # keep ONLY the top one
  ungroup()


## ---- 7. Create master list of samples droppes --------

# The dup/MZ pair (e.g. R3057101 / R4323048) was excluded up front. Pull both
# members of any dup/MZ pair from the original kinship table.
dup_pairs <- kin %>% filter(relationship == 'dup/MZ')
dup_ids   <- unique(c(dup_pairs$ID1, dup_pairs$ID2))
cat("Duplicate/swap samples removed at start:", length(dup_ids), "\n")
print(dup_ids)

## Relatedness drops from DISJOINT PAIRS 
# the member NOT kept in each disjoint pair
disjoint_drop_ids <- disjoint_decision %>% filter(!keep) %>% pull(sample.id)
cat("Relatedness drops from disjoint pairs:", length(disjoint_drop_ids), "\n")
print(disjoint_drop_ids)

## Relatedness drops from Triads 
# the member NOT kept in each triad
triad_drop_ids <- triad_decision %>% filter(!keep) %>% pull(sample.id)
cat("Relatedness drops from triads:", length(triad_drop_ids), "\n")
print(triad_drop_ids)

cat("Relatedness pairs dropped due to specimen count:",
 length(disjoint_decision$decided_by_specimens[disjoint_decision$decided_by_specimens == TRUE])/2 +  
 length(triad_decision$decided_by_specimens[triad_decision$decided_by_specimens == TRUE])/3, "\n")

all_related_ids <- unique(c(related_ids, dup_ids))   # every sample in any flagged pair
 
reason_for <- function(id) {
  if (id %in% dup_ids)            return("duplicate/swap")
  if (id %in% contam_related)     return("contamination+realtedness")
  if (id %in% disjoint_drop_ids)  return("relatedness")
  if (id %in% triad_drop_ids)     return("relatedness")
  return(NA_character_)   # kept
}
 
kinship_drop_flags <- data.frame(sample.id = all_sample_ids,
                                 stringsAsFactors = FALSE) %>%
  rowwise() %>%
  mutate(drop_reason = reason_for(sample.id),
         kinship_related_drop = !is.na(drop_reason)) %>%
  ungroup()
 

write.table(kinship_drop_flags, '~/data/pca_results/QC/Sample_relatedness_drop_flag.csv', sep = ',')


dim(triad_decision)

contam_related

dat[dat$FID == 'NYBB_348WGS',]
