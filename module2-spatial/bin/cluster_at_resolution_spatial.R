#!/usr/bin/env Rscript
# Cluster spots at ONE resolution, score against pathologist Classification

# Usage: cluster_at_resolution_spatial.R <base_spatial.rds> <resolution>

args <- commandArgs(trailingOnly = TRUE)
seurat_obj <- readRDS(args[1])
resolution <- as.numeric(args[2])

suppressMessages({ library(Seurat); library(mclust) })

seurat_obj <- FindClusters(seurat_obj, resolution = resolution, verbose = FALSE)
n_clusters <- length(unique(Idents(seurat_obj)))

# Exclude spots with no usable pathology label (NA, "Artefact", or a stray

valid <- !is.na(seurat_obj$Classification) &
  !(seurat_obj$Classification %in% c("Artefact", "Classification"))

ari <- adjustedRandIndex(
  Idents(seurat_obj)[valid],
  seurat_obj$Classification[valid]
)

res_tag <- gsub("\\.", "p", sprintf("%.2f", resolution))
saveRDS(seurat_obj, paste0("clustered_res", res_tag, ".rds"))

result <- data.frame(resolution = resolution, n_clusters = n_clusters, ari = ari, n_valid_spots = sum(valid))
write.csv(result, paste0("ari_res", res_tag, ".csv"), row.names = FALSE)

cat(sprintf("resolution=%.2f  clusters=%d  ARI=%.3f  (n=%d valid spots)\n",
            resolution, n_clusters, ari, sum(valid)))
