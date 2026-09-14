#!/usr/bin/env Rscript
# Build a spatial neighbor graph from the real Visium spot coordinates.
# Visium's hex grid gives each interior spot exactly 6 nearest neighbors --
# we use k-nearest-neighbors (k=6) on the actual pixel coordinates rather
# than a black-box built-in function, so the neighbor definition is
# explicit and inspectable.
# Usage: build_spatial_neighbors.R <spatial_rds> <k_neighbors>

args <- commandArgs(trailingOnly = TRUE)
seurat_obj <- readRDS(args[1])
k <- as.numeric(args[2])

suppressMessages({ library(Seurat); library(spdep) })

coords <- GetTissueCoordinates(seurat_obj)
coords <- coords[, c("x", "y")]  # keep pixel coordinates; row names are
                                   # already correct barcodes from
                                   # GetTissueCoordinates itself -- do not
                                   # override them.

# Restrict to spots that actually have a CNV score. Some spots in this
# object never went through inferCNV in Module 4 (they lacked a usable
# pathology label at that stage), leaving cnv_score = NA -- the neighbor
# graph must be built on exactly the same spot set we'll later test, or
# the two will be misaligned downstream.
has_score <- !is.na(seurat_obj$cnv_score)
valid_spots <- intersect(rownames(coords), colnames(seurat_obj)[has_score])
cat("Spots with valid coordinates AND a CNV score:", length(valid_spots),
    "(of", nrow(coords), "with coordinates,", sum(has_score), "with a CNV score)\n")
coords <- coords[valid_spots, ]

cat("Building k=", k, " nearest-neighbor spatial graph for", nrow(coords), "spots\n")

knn_result <- knearneigh(as.matrix(coords), k = k)
nb <- knn2nb(knn_result)
listw <- nb2listw(nb, style = "W")  # row-standardized weights

saveRDS(list(nb = nb, listw = listw, coords = coords), "spatial_neighbors.rds")
cat("Spatial neighbor graph built and saved.\n")
