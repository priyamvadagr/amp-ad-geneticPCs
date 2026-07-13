################################################################################
# Scree_plot_and_TW_test.r
# Determine the number of informative PCs via:
#   1. Scree plot (eigenvalues / % variance by PC) with elbow inspection
#   2. Tracy-Widom test (Patterson, Price & Reich 2006) for significant PCs
#
# Input: pca_result.rds (or pca_eigenvalues.csv) from run_pca.R
################################################################################

library(ggplot2)
library(dplyr)
library(RMTstat)

out_dir <- "~/data/pca_results/PCA"
fig_dir <- "/home/ec2-user/AMP-AD_genetic_PCs/Results/PCA"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

pca <- readRDS(file.path(out_dir, "pca_result.rds"))

## ---- 1. Eigenvalue / variance table ----------------------------------------
eig <- data.frame(
  PC         = seq_along(pca$eigenval),
  eigenvalue = pca$eigenval,
  var_pct    = pca$varprop * 100
) %>% filter(!is.na(eigenvalue) & eigenvalue > 0)

n_show <- min(32, nrow(eig))   # plot the top PCs


var_tbl <- data.frame(
  PC          = seq_along(pca$eigenval),
  eigenvalue  = pca$eigenval,
  var_pct     = pca$varprop * 100,
  cum_var_pct = cumsum(pca$varprop * 100)
)
var_tbl <- var_tbl[!is.na(var_tbl$eigenvalue), ]

write.csv(var_tbl, file.path(out_dir, "pca_variance_explained.csv"), row.names = FALSE)
head(var_tbl, 10)

## ---- 2. SCREE PLOTS --------------------------------------------------------
# (a) % variance explained
p_var <- ggplot(eig[1:n_show, ], aes(PC, var_pct)) +
  geom_line(colour = "grey50") +
  geom_point(size = 2, colour = "#377EB8") +
  geom_vline(xintercept = 6, color = 'red', linetype = 'dashed') +
  scale_x_continuous(breaks = seq(0, n_show, 2)) +
  labs(title = "Scree plot: variance explained per PC",
       x = "Principal component", y = "% variance explained") +
  theme_bw(base_size = 12)
ggsave(file.path(fig_dir, "scree_variance.png"), p_var,
       width = 8, height = 5, dpi = 150)

# (b) eigenvalues on log scale (clarifies the noise-floor elbow)
p_eig <- ggplot(eig[1:n_show, ], aes(PC, eigenvalue)) +
  geom_line(colour = "grey50") +
  geom_point(size = 2, colour = "#E41A1C") +
  geom_vline(xintercept = 6, color = 'red', linetype = 'dashed') +
  scale_y_log10() +
  scale_x_continuous(breaks = seq(0, n_show, 2)) +
  labs(title = "Scree plot: eigenvalues (log scale)",
       subtitle = "Elbow = transition to the noise floor; PCs before it are informative",
       x = "Principal component", y = "Eigenvalue (log10)") +
  theme_bw(base_size = 12)
ggsave(file.path(fig_dir, "scree_eigenvalue_log.png"), p_eig,
       width = 8, height = 5, dpi = 150)

## ---- # Quick elbow heuristic (for cross-check) ----------------------------
# Proportion-of-variance drop: flag where successive PC variance ratio levels off
eig <- eig %>% mutate(var_ratio = var_pct / lead(var_pct))
cat("\nTop PCs and % variance (eyeball the elbow):\n")
print(head(eig[, c("PC","eigenvalue","var_pct", "var_ratio")], 20))
cat("\nScree plots written to", fig_dir, "\n")
cat("Decide # informative PCs from: Tracy-Widom count + scree elbow + visual\n")
cat("inspection of PC scatter plots (whether each PC shows real structure).\n")

## ---- 3. Subset eigenvector file to significant PCs ----------------------------

pc_df <- read.table(file.path(out_dir, "pca_eigenvectors.csv"), header = T, sep = ',')
sig_pc_df <- pc_df[, c('sample.id', 'PC1', 'PC2', 'PC3', 'PC4')]

write.table(sig_pc_df, file.path(out_dir, "sig_pca_eigenvectors.csv"), row.names = F)




