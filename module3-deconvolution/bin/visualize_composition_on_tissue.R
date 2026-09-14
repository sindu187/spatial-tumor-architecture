#!/usr/bin/env Rscript
# Visualize the anchor-based predicted dominant cell type per spot,
# overlaid on the real tissue image, alongside the pathology ground
# truth we already validated in Module 2 
#
# Usage: visualize_composition_on_tissue.R <spatial_with_anchor_rds> <patient_id>

args <- commandArgs(trailingOnly = TRUE)
spatial_obj <- readRDS(args[1])
patient_id <- args[2]

suppressMessages({ library(Seurat); library(ggplot2); library(patchwork) })

p1 <- SpatialDimPlot(spatial_obj, group.by = "predicted.id") +
  ggtitle("Predicted dominant cell type (Method A: anchor transfer)")

p2 <- SpatialDimPlot(spatial_obj, group.by = "Classification") +
  ggtitle("Pathologist annotation (ground truth, from Module 2)")

combined <- p1 + p2 + plot_annotation(title = paste0(patient_id, ": predicted cell type vs. pathology"))
ggsave("predicted_celltype_vs_pathology.png", plot = combined, width = 16, height = 8, dpi = 150)

# Also show a couple of individual cell-type prediction score heatmaps,
# which give a continuous (not just dominant-label) view of composition
score_cols <- grep("^prediction\\.score\\.", names(spatial_obj@meta.data), value = TRUE)
score_cols <- setdiff(score_cols, "prediction.score.max")
featured <- head(score_cols, 4)  # first 4 for a readable panel

p3 <- SpatialFeaturePlot(spatial_obj, features = featured, ncol = 2)
ggsave("celltype_score_heatmaps.png", plot = p3, width = 12, height = 10, dpi = 150)

saveRDS(spatial_obj, "final_deconvolved.rds")
cat("Module 3 complete for", patient_id, "\n")
