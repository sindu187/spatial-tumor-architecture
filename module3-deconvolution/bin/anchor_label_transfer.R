#!/usr/bin/env Rscript
# Method A: anchor-based label transfer. pairs of
# mutually-similar cells between the scRNA reference and spatial query --
# then uses them to transfer cell-type identity, producing a prediction
# score per cell type for every spot (roughly, an estimated proportion).
#
# Both reference and query are put through standard log-normalization
# (rather than the SCTransform assay used for spatial clustering in
# Module 2) to keep the two datasets on a directly comparable footing
# for the anchor-finding step -- a standard simplification for
# cross-modality integration.
#
# Usage: anchor_label_transfer.R <scrna_rds> <spatial_rds>

args <- commandArgs(trailingOnly = TRUE)
scrna_obj <- readRDS(args[1])
spatial_obj <- readRDS(args[2])

suppressMessages(library(Seurat))

DefaultAssay(scrna_obj) <- "RNA"  # already log-normalized in Module 1

DefaultAssay(spatial_obj) <- "Spatial"
spatial_obj <- NormalizeData(spatial_obj, verbose = FALSE)

cat("Finding transfer anchors between scRNA reference and spatial query...\n")
anchors <- FindTransferAnchors(
  reference = scrna_obj,
  query = spatial_obj,
  normalization.method = "LogNormalize",
  verbose = FALSE
)

cat("Transferring cell-type labels...\n")
predictions <- TransferData(
  anchorset = anchors,
  refdata = scrna_obj$our_annotation,
  verbose = FALSE
)

write.csv(predictions, "anchor_predictions.csv", row.names = TRUE)

spatial_obj <- AddMetaData(spatial_obj, metadata = predictions)
saveRDS(spatial_obj, "spatial_with_anchor.rds")

cat("\n--- Predicted dominant cell type per spot (counts) ---\n")
print(table(spatial_obj$predicted.id))
cat("\nDone.\n")
