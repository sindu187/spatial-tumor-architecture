#!/usr/bin/env Rscript
# Filter spots by QC thresholds.
# Usage: qc_filter_spatial.R <raw_spatial.rds> <min_features> <max_mito>

args <- commandArgs(trailingOnly = TRUE)
seurat_obj <- readRDS(args[1])
min_features <- as.numeric(args[2])
max_mito <- as.numeric(args[3])

suppressMessages(library(Seurat))

cat("Spots before filtering:", ncol(seurat_obj), "\n")
seurat_obj <- subset(seurat_obj, subset = nFeature_Spatial > min_features & percent.mito < max_mito)
cat("Spots after filtering:", ncol(seurat_obj), "\n")

saveRDS(seurat_obj, "filtered_spatial.rds")
