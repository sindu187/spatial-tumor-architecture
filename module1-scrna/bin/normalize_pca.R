#!/usr/bin/env Rscript
# Normalize, find variable features, scale, PCA, build neighbor graph.
# This produces the "base" object that every resolution in the sweep
# branches off from -- computed once, reused by all parallel branches.
# Usage: normalize_pca.R <filtered.rds> <npcs>

args <- commandArgs(trailingOnly = TRUE)
seurat_obj <- readRDS(args[1])
npcs <- as.numeric(args[2])

suppressMessages(library(Seurat))

seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = 2000, verbose = FALSE)
seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
seurat_obj <- RunPCA(seurat_obj, npcs = npcs, verbose = FALSE)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:npcs, verbose = FALSE)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:npcs, verbose = FALSE)

saveRDS(seurat_obj, "base_object.rds")
cat("Base object ready (PCA + neighbor graph computed).\n")
