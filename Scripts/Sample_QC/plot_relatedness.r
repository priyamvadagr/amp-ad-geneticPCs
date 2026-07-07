## ---- Relatedness diagnostic: IBS0 vs kinship, with relationship zones ------
library(ggplot2)
library(dplyr)

kin <- read.table('~/data/pca_results/QC/Sample_relatedness_flag.csv',
 sep = ',', header = T)

# Standard KING relationship thresholds
thr <- c(unrelated = 0.0442, third = 0.0884, second = 0.177, first = 0.354)
 
# Assign each pair to a tier (reuse the 'relationship' factor if already present)
kin <- kin %>%
  mutate(tier = cut(
    kinship,
    breaks = c(-Inf, 0.0442, 0.0884, 0.177, 0.354, Inf),
    labels = c("unrelated (<0.0442)",
               "3rd-degree (0.0442-0.0884)",
               "2nd-degree (0.0884-0.177)",
               "1st-degree (0.177-0.354)",
               "duplicate/MZ (>0.354)")
  ))
 
tier_cols <- c(
  "unrelated (<0.0442)"        = "grey75",
  "3rd-degree (0.0442-0.0884)" = "#92C5DE",
  "2nd-degree (0.0884-0.177)"  = "#F4A582",
  "1st-degree (0.177-0.354)"   = "#D6604D",
  "duplicate/MZ (>0.354)"      = "#B2182B"
)
 
## ---- Plot 1: full distribution (log y-axis so rare tiers are visible) ------
# A log y-axis is essential: the 3rd-degree bin has ~220k pairs while the
# higher tiers have <25, so on a linear axis the real relatives are invisible.
p_full <- ggplot(kin, aes(x = kinship, fill = tier)) +
  geom_histogram(binwidth = 0.005, colour = NA) +
  geom_vline(xintercept = thr, linetype = "dashed", colour = "grey30") +
  annotate("text", x = thr, y = Inf,
           label = c("3rd","2nd","1st","dup"),
           vjust = 1.5, hjust = -0.1, size = 3, colour = "grey30") +
  scale_fill_manual(values = tier_cols, name = "Relationship") +
  scale_y_log10() +
  labs(title = "Pairwise KING kinship by relationship tier",
       subtitle = "log y-axis; dashed lines = tier thresholds",
       x = "Kinship coefficient", y = "Pair count (log scale)") +
  theme_bw(base_size = 12) +
  theme(legend.position = "right")
 
ggsave(file.path(plots_dir, "kinship_tiers_full.png"),
       p_full, width = 9, height = 5, dpi = 150)
 
## ---- Plot 2: zoom on the relatedness range (kinship > 0.0442) --------------
# Drops the giant 'unrelated' mass and the structure-inflated bulk near the
# threshold so the genuine relatives (2nd/1st/dup) are clearly visible.
p_zoom <- kin %>%
  filter(kinship > 0.0442) %>%
  ggplot(aes(x = kinship, fill = tier)) +
  geom_histogram(binwidth = 0.005, colour = NA) +
  geom_vline(xintercept = thr[c("third","second","first")],
             linetype = "dashed", colour = "grey30") +
  scale_fill_manual(values = tier_cols, name = "Relationship") +
  scale_y_log10() +
  labs(title = "Kinship for flagged pairs (> 0.0442), zoomed",
       subtitle = "The 3rd-degree pile is structure inflation; 2nd-degree+ are genuine relatives",
       x = "Kinship coefficient", y = "Pair count (log scale)") +
  theme_bw(base_size = 12)
 
ggsave(file.path(plots_dir, "kinship_tiers_zoom.png"),
       p_zoom, width = 9, height = 5, dpi = 150)


p_new <- ggplot(kin, aes(x = IBS0, y = kinship, colour = tier)) +
  geom_point(alpha = 0.4, size = 1) +
  geom_hline(yintercept = c(0.0884, 0.177, 0.354), linetype = "dashed")
  ggsave(file.path(plots_dir, "kinship_tiers_zoom_v2.png"),
       p_new, width = 9, height = 5, dpi = 150)
 
## ---- Counts per tier (for the figure caption / methods) --------------------
tier_counts <- kin %>% count(tier)
print(tier_counts)
write.csv(tier_counts, file.path(out_dir, "kinship_tier_counts.csv"),
          row.names = FALSE)
 
cat("\nPlots written to", out_dir, "\n")

## ---- Plot 3: classify each pair by BOTH kinship and IBS0 ------------------------
# Standard KING decision rules:
#   duplicate/MZ : kinship > 0.354
#   1st-degree   : 0.177 < kinship <= 0.354
#       -> parent-offspring (PO): IBS0 ~ 0 (very small)
#       -> full-sib (FS)        : IBS0 clearly > 0
#   2nd-degree   : 0.0884 < kinship <= 0.177
#   3rd-degree   : 0.0442 < kinship <= 0.0884
#   unrelated    : kinship <= 0.0442
#
# IBS0 cutoff separating PO from FS within the 1st-degree band.
# KING's own heuristic: PO have IBS0 below ~0.1 * (something); in practice a
# small absolute threshold works. Tune by eye from your data (see plot).
ibs0_po_cut <- 0.0025   # pairs below this within 1st-degree => parent-offspring

kin <- kin %>%
  mutate(
    relationship = case_when(
      kinship > 0.354                       ~ "duplicate/MZ",
      kinship > 0.177 & IBS0 <  ibs0_po_cut ~ "1st: parent-offspring",
      kinship > 0.177 & IBS0 >= ibs0_po_cut ~ "1st: full-sib",
      kinship > 0.0884                      ~ "2nd-degree",
      kinship > 0.0442                      ~ "3rd-degree",
      TRUE                                  ~ "unrelated"
    )
  )

rel_levels <- c("unrelated","3rd-degree","2nd-degree",
                "1st: full-sib","1st: parent-offspring","duplicate/MZ")
kin$relationship <- factor(kin$relationship, levels = rel_levels)

rel_cols <- c(
  "unrelated"             = "grey80",
  "3rd-degree"            = "#4DAF4A",  # green
  "2nd-degree"            = "#377EB8",  # blue
  "1st: full-sib"         = "#FF7F00",  # orange
  "1st: parent-offspring" = "#984EA3",  # purple
  "duplicate/MZ"          = "#E41A1C"   # red
)
 

## ---- 2. The diagnostic scatter with annotated zones ------------------------
# Plot only the related pairs (kinship > 0.0442) so the structure is visible;
# the ~220k unrelated pairs would just be a smear at the bottom-right.
plot_df <- kin %>% filter(kinship > 0.0442)

p <- ggplot(plot_df, aes(x = IBS0, y = kinship, colour = relationship)) +
  # kinship tier guides
  geom_hline(yintercept = c(0.0884, 0.177, 0.354),
             linetype = "dashed", colour = "grey50") +
  # PO/FS divider within the 1st-degree band
  geom_vline(xintercept = ibs0_po_cut,
             linetype = "dotted", colour = "grey40") +
  geom_point(size = 2.4, alpha = 0.8) +
  scale_colour_manual(values = rel_cols, name = "Inferred relationship",
                      drop = FALSE) +
  # zone labels
  annotate("text", x = max(plot_df$IBS0)*0.95, y = 0.066,
           label = "3rd-degree", hjust = 1, size = 3, colour = "grey40") +
  annotate("text", x = max(plot_df$IBS0)*0.95, y = 0.13,
           label = "2nd-degree", hjust = 1, size = 3, colour = "grey40") +
  annotate("text", x = ibs0_po_cut/2, y = 0.40,
           label = "dup/MZ", size = 3, colour = "grey40") +
  annotate("text", x = ibs0_po_cut*0.4, y = 0.26,
           label = "parent-\noffspring", size = 3, colour = "grey40") +
  annotate("text", x = ibs0_po_cut*2.2, y = 0.26,
           label = "full-sib", size = 3, colour = "grey40") +
  labs(title = "Relatedness diagnostic: IBS0 vs kinship",
       subtitle = "Parent-offspring (IBS0~0) separate from full-sibs (IBS0>0) at the same kinship",
       x = "IBS0 (proportion of markers sharing zero alleles)",
       y = "KING kinship coefficient") +
  theme_bw(base_size = 12)

ggsave(file.path(plots_dir, "relatedness_IBS0_vs_kinship.png"),
       p, width = 9, height = 6, dpi = 150)

## ---- 3. Report the classified relative pairs -------------------------------
classified <- kin %>%
  filter(kinship > 0.0884) %>%               # 2nd-degree and closer (the real ones)
  arrange(desc(kinship)) %>%
  select(ID1, ID2, IBS0, kinship, relationship)

print(classified, n = 50)
write.csv(classified, file.path(out_dir, "classified_relative_pairs.csv"),
          row.names = FALSE)

cat("\nRelationship counts (kinship > 0.0884):\n")
print(table(classified$relationship))
cat("\nPlot written to", file.path(out_dir, "relatedness_IBS0_vs_kinship.png"), "\n")