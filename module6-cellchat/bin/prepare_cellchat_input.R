#!/usr/bin/env Rscript
# Prepare CellChat's required inputs from Module 3's deconvolution output:
# normalized expression matrix, cell-type group labels (Module 3's
# anchor-based predicted.id per spot), spatial coordinates, and the
# Space Ranger scale factors CellChat needs for its spatial (contact-
# range-aware) communication mode.
# Usage: prepare_cellchat_input.R <deconvolved_rds> <scalefactors_json>

args <- commandArgs(trailingOnly = TRUE)
seurat_obj <- readRDS(args[1])
scalefactors <- jsonlite::fromJSON(args[2])

suppressMessages(library(Seurat))

DefaultAssay(seurat_obj) <- "Spatial"
seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
data_matrix <- GetAssayData(seurat_obj, assay = "Spatial", layer = "data")

meta <- data.frame(
  group = seurat_obj$predicted.id,
  row.names = colnames(seurat_obj)
)

coords <- GetTissueCoordinates(seurat_obj)[, c("x", "y")]
coords <- as.matrix(coords)

# Restrict everything to spots with a valid predicted cell-type group
valid_spots <- intersect(rownames(coords), rownames(meta)[!is.na(meta$group)])
data_matrix <- data_matrix[, valid_spots]
meta <- meta[valid_spots, , drop = FALSE]
coords <- coords[valid_spots, ]

cat("Prepared CellChat input:", nrow(data_matrix), "genes x", ncol(data_matrix), "spots\n")
cat("\n--- Spots per predicted cell-type group ---\n")
print(table(meta$group))

saveRDS(
  list(data_matrix = data_matrix, meta = meta, coords = coords, scalefactors = scalefactors),
  "cellchat_input.rds"
)
