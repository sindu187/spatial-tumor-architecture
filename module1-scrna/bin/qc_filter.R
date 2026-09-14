#!/usr/bin/env Rscript
# Filter cells by QC thresholds.
# Usage: qc_filter.R <raw_unfiltered.rds> <min_features> <max_features> <max_mito>

args <- commandArgs(trailingOnly = TRUE)
seurat_obj <- readRDS(args[1])
min_features <- as.numeric(args[2])
max_features <- as.numeric(args[3])
max_mito <- as.numeric(args[4])

suppressMessages(library(Seurat))

cat("Cells before filtering:", ncol(seurat_obj), "\n")
seurat_obj <- subset(
  seurat_obj,
  subset = nFeature_RNA > min_features & nFeature_RNA < max_features & percent.mito < max_mito
)
cat("Cells after filtering:", ncol(seurat_obj), "\n")

saveRDS(seurat_obj, "filtered.rds")
