#!/usr/bin/env Rscript
# Plot our independent clusters and the pathologist's own annotations,
# both overlaid on the real tissue image, side by side.
# Usage: visualize_on_tissue.R <best_clustered_spatial.rds> <patient_id>

args <- commandArgs(trailingOnly = TRUE)
seurat_obj <- readRDS(args[1])
patient_id <- args[2]

suppressMessages({ library(Seurat); library(ggplot2); library(patchwork) })

p1 <- SpatialDimPlot(seurat_obj, group.by = "seurat_clusters", label = TRUE, label.size = 3) +
  ggtitle("Our independent clusters")

p2 <- SpatialDimPlot(seurat_obj, group.by = "Classification") +
  ggtitle("Pathologist annotation (ground truth)")

combined <- p1 + p2 + plot_annotation(title = paste0(patient_id, ": expression clusters vs. pathology"))
ggsave("clusters_vs_pathology.png", plot = combined, width = 16, height = 8, dpi = 150)

saveRDS(seurat_obj, "final_spatial_annotated.rds")
cat("Module 2 complete for", patient_id, "\n")
