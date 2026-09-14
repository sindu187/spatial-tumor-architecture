#!/usr/bin/env Rscript
# Cluster the base object at ONE resolution and score it against the
# published celltype_major labels. Designed to run as one of several
# PARALLEL process instances (one per resolution) in Nextflow.
# Usage: cluster_at_resolution.R <base_object.rds> <resolution>

args <- commandArgs(trailingOnly = TRUE)
seurat_obj <- readRDS(args[1])
resolution <- as.numeric(args[2])

suppressMessages({ library(Seurat); library(mclust) })

seurat_obj <- FindClusters(seurat_obj, resolution = resolution, verbose = FALSE)
n_clusters <- length(unique(Idents(seurat_obj)))
ari <- adjustedRandIndex(Idents(seurat_obj), seurat_obj$celltype_major)

# Tag the output filenames with the resolution so parallel runs don't collide
res_tag <- gsub("\\.", "p", sprintf("%.2f", resolution))
saveRDS(seurat_obj, paste0("clustered_res", res_tag, ".rds"))

result <- data.frame(resolution = resolution, n_clusters = n_clusters, ari = ari)
write.csv(result, paste0("ari_res", res_tag, ".csv"), row.names = FALSE)

cat(sprintf("resolution=%.2f  clusters=%d  ARI=%.3f\n", resolution, n_clusters, ari))
