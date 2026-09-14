#!/usr/bin/env Rscript
# Visualize the spatial autocorrelation result: (1) a Moran scatterplot
# -- the classic way to visualize spatial autocorrelation, plotting each
# spot's own value against the average of its neighbors' values -- and
# (2) the CNV score map itself, annotated with the test statistic.
# Usage: visualize_spatial_clustering.R <moran_data_for_plotting.rds> <moran_test_results.csv> <spatial_rds> <patient_id>

args <- commandArgs(trailingOnly = TRUE)
moran_data <- readRDS(args[1])
results <- read.csv(args[2])
seurat_obj <- readRDS(args[3])
patient_id <- args[4]

suppressMessages({ library(Seurat); library(spdep); library(ggplot2); library(patchwork) })

cnv_row <- results[results$variable == "cnv_score", ]

# Moran scatterplot: each spot's value (x) vs. its spatial lag -- the
# average value among its spatial neighbors (y). A clear positive slope
# is the classic visual signature of spatial autocorrelation.
png("moran_scatterplot.png", width = 700, height = 600)
moran.plot(
  moran_data$cnv_score, moran_data$listw,
  xlab = "CNV score (spot)", ylab = "CNV score (spatial lag: neighbor average)",
  main = sprintf("Moran scatterplot: CNV score\nMoran's I = %.3f, p = %.4f",
                  cnv_row$morans_i, cnv_row$p_value)
)
dev.off()

# Spatial map for context, alongside pathology
p1 <- SpatialFeaturePlot(seurat_obj, features = "cnv_score") +
  ggtitle(sprintf("CNV score (Moran's I = %.3f, p = %.4f)", cnv_row$morans_i, cnv_row$p_value))

p2 <- SpatialDimPlot(seurat_obj, group.by = "Classification") +
  ggtitle("Pathologist annotation")

combined <- p1 + p2 + plot_annotation(title = paste0(patient_id, ": spatial clustering of CNV signal"))
ggsave("cnv_spatial_clustering.png", plot = combined, width = 16, height = 8, dpi = 150)

cat("Module 5 complete for", patient_id, "\n")
