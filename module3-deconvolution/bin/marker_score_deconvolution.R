#!/usr/bin/env Rscript
# Method B: marker-gene module scoring. For each cell type, takes its top
# marker genes (from Module 1's independently validated FindAllMarkers
# results) as a signature gene set, then scores how strongly each spot's
# expression matches that signature (Seurat's AddModuleScore). This is a
# genuinely different statistical approach from Method A's anchor-based
# integration -- no reference-cell matching, just direct signature
# scoring -- making it a meaningful independent comparison.
#
# Usage: marker_score_deconvolution.R <spatial_rds> <all_cluster_markers.csv> <cluster_labels.csv> <n_marker_genes>

args <- commandArgs(trailingOnly = TRUE)
spatial_obj <- readRDS(args[1])
markers <- read.csv(args[2])
cluster_labels <- read.csv(args[3])
n_genes <- as.numeric(args[4])

suppressMessages({ library(Seurat); library(dplyr) })

# Attach the final cell-type label to each cluster's markers
markers <- merge(markers, cluster_labels, by.x = "cluster", by.y = "seurat_clusters")

# Top N marker genes per cell-type label
gene_sets <- markers %>%
  group_by(assigned_label) %>%
  slice_max(avg_log2FC, n = n_genes) %>%
  summarise(genes = list(gene))

gene_set_list <- setNames(gene_sets$genes, gene_sets$assigned_label)

cat("Gene sets built for", length(gene_set_list), "cell types (", n_genes, "genes each)\n")

DefaultAssay(spatial_obj) <- "Spatial"
spatial_obj <- NormalizeData(spatial_obj, verbose = FALSE)

spatial_obj <- AddModuleScore(
  spatial_obj,
  features = gene_set_list,
  name = "score_",
  seed = 42
)

# AddModuleScore names columns score_1, score_2, ... in input order --
# rename them back to the actual cell-type labels for clarity
score_cols <- paste0("score_", seq_along(gene_set_list))
names(spatial_obj@meta.data)[names(spatial_obj@meta.data) %in% score_cols] <- names(gene_set_list)

marker_scores <- spatial_obj@meta.data[, names(gene_set_list)]
write.csv(marker_scores, "marker_scores.csv", row.names = TRUE)

saveRDS(spatial_obj, "spatial_with_marker_scores.rds")
cat("Done. Marker-based scores computed for", nrow(marker_scores), "spots.\n")
