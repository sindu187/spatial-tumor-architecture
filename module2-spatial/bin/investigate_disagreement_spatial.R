#!/usr/bin/env Rscript
# Check per-cluster purity against pathologist Classification labels --
# same approach as module1-scrna's investigation, applied to spatial spots.
# Usage: investigate_disagreement_spatial.R <best_clustered_spatial.rds>

args <- commandArgs(trailingOnly = TRUE)
seurat_obj <- readRDS(args[1])

suppressMessages({ library(Seurat); library(dplyr) })

meta <- seurat_obj@meta.data
meta$Classification[meta$Classification %in% c("Artefact", "Classification")] <- NA

cluster_purity <- meta %>%
  filter(!is.na(Classification)) %>%
  group_by(seurat_clusters, Classification) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(seurat_clusters) %>%
  mutate(cluster_total = sum(n), fraction = n / cluster_total) %>%
  slice_max(fraction, n = 1) %>%
  arrange(fraction)

write.csv(cluster_purity, "cluster_purity_spatial.csv", row.names = FALSE)
cat("--- Per-cluster purity against pathologist Classification ---\n")
print(as.data.frame(cluster_purity))

# for each cluster, ALL pathology labels present (not just
# the dominant one) -- shows whether "impure" clusters split cleanly
# between 2-3 real categories, or are genuinely scattered/ambiguous.
full_breakdown <- meta %>%
  filter(!is.na(Classification)) %>%
  group_by(seurat_clusters, Classification) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(seurat_clusters) %>%
  mutate(cluster_total = sum(n), fraction = n / cluster_total) %>%
  arrange(seurat_clusters, desc(fraction))

write.csv(full_breakdown, "full_breakdown_spatial.csv", row.names = FALSE)
cat("\n--- Full per-cluster breakdown (all pathology labels present) ---\n")
print(as.data.frame(full_breakdown))

saveRDS(seurat_obj, "investigated_spatial.rds")
