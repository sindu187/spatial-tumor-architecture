#!/usr/bin/env Rscript
# Find marker genes for every cluster (independent biological validation
# of cluster identity.
# Usage: find_markers.R <investigated.rds>

args <- commandArgs(trailingOnly = TRUE)
seurat_obj <- readRDS(args[1])

suppressMessages({ library(Seurat); library(dplyr) })
Idents(seurat_obj) <- "seurat_clusters"

all_markers <- FindAllMarkers(
  seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5, verbose = FALSE
)
write.csv(all_markers, "all_cluster_markers.csv", row.names = FALSE)

top_markers <- all_markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 5) %>%
  summarise(top_genes = paste(gene, collapse = ", "))

cat("--- Top 5 markers per cluster ---\n")
print(as.data.frame(top_markers))

saveRDS(seurat_obj, "with_markers.rds")
