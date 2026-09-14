#!/usr/bin/env Rscript
# SCTransform (the recommended normalization for Visium spatial data --
# handles spot-to-spot technical variance better than standard log-
# normalization), then PCA and neighbor graph.
# Usage: sctransform_pca.R <filtered_spatial.rds> <npcs>

args <- commandArgs(trailingOnly = TRUE)
seurat_obj <- readRDS(args[1])
npcs <- as.numeric(args[2])

suppressMessages(library(Seurat))

seurat_obj <- SCTransform(seurat_obj, assay = "Spatial", verbose = FALSE)
seurat_obj <- RunPCA(seurat_obj, assay = "SCT", npcs = npcs, verbose = FALSE)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:npcs, verbose = FALSE)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:npcs, verbose = FALSE)

saveRDS(seurat_obj, "base_spatial.rds")
cat("Base spatial object ready (SCTransform + PCA + neighbor graph).\n")
