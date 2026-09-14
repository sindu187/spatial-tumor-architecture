#!/usr/bin/env Rscript
# Plot the inferred CNV score spatially, alongside pathology annotation.
# Usage: visualize_cnv_on_tissue.R <spatial_rds> <cnv_scores_csv> <patient_id>

args <- commandArgs(trailingOnly = TRUE)
seurat_obj <- readRDS(args[1])
scores_df <- read.csv(args[2])
patient_id <- args[3]

suppressMessages({ library(Seurat); library(ggplot2); library(patchwork) })

rownames(scores_df) <- scores_df$spot
seurat_obj$cnv_score <- NA
common <- intersect(colnames(seurat_obj), rownames(scores_df))
seurat_obj$cnv_score[common] <- scores_df[common, "cnv_score"]

p1 <- SpatialFeaturePlot(seurat_obj, features = "cnv_score") +
  ggtitle("Inferred CNV score (deviation from neutral)")

p2 <- SpatialDimPlot(seurat_obj, group.by = "Classification") +
  ggtitle("Pathologist annotation (ground truth)")

combined <- p1 + p2 + plot_annotation(title = paste0(patient_id, ": inferred CNV vs. pathology"))
ggsave("cnv_vs_pathology.png", plot = combined, width = 16, height = 8, dpi = 150)

saveRDS(seurat_obj, "final_cnv_annotated.rds")
cat("Module 4 complete for", patient_id, "\n")
