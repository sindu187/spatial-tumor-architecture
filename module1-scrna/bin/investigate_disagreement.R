#!/usr/bin/env Rscript
# Check per-cluster purity against published labels, and whether any
# "impure" clusters resolve cleanly at the finer celltype_minor level
# (i.e. real sub-structure, not misclassification).
# Usage: investigate_disagreement.R <best_clustered.rds>

args <- commandArgs(trailingOnly = TRUE)
seurat_obj <- readRDS(args[1])

suppressMessages({ library(Seurat); library(dplyr) })

cluster_purity <- seurat_obj@meta.data %>%
  group_by(seurat_clusters, celltype_major) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(seurat_clusters) %>%
  mutate(cluster_total = sum(n), fraction = n / cluster_total) %>%
  slice_max(fraction, n = 1) %>%
  arrange(fraction)

write.csv(cluster_purity, "cluster_purity.csv", row.names = FALSE)
cat("--- Per-cluster purity ---\n")
print(as.data.frame(cluster_purity))

minor_check <- seurat_obj@meta.data %>%
  group_by(seurat_clusters, celltype_minor) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(seurat_clusters) %>%
  mutate(cluster_total = sum(n), fraction = n / cluster_total) %>%
  arrange(seurat_clusters, desc(fraction))

write.csv(minor_check, "minor_type_breakdown.csv", row.names = FALSE)

saveRDS(seurat_obj, "investigated.rds")
cat("\nSaved cluster_purity.csv, minor_type_breakdown.csv, investigated.rds\n")
