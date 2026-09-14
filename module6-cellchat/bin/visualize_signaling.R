#!/usr/bin/env Rscript
# Visualize CellChat's inferred communication network: an aggregated
# circle plot (overall signaling strength between cell-type groups) and
# a bubble plot for the specific pathways of interest.
# Usage: visualize_signaling.R <cellchat_obj.rds> <patient_id>

args <- commandArgs(trailingOnly = TRUE)
cellchat <- readRDS(args[1])
patient_id <- args[2]

suppressMessages(library(CellChat))

png("aggregated_communication_network.png", width = 900, height = 900)
groupSize <- as.numeric(table(cellchat@idents))
netVisual_circle(
  cellchat@net$weight, vertex.weight = groupSize, weight.scale = TRUE,
  label.edge = FALSE, title.name = paste0(patient_id, ": aggregated cell-cell communication")
)
dev.off()

pathways_of_interest <- c("PD-L1", "CSF", "CXCL", "TGFb")
available_pathways <- intersect(pathways_of_interest, cellchat@netP$pathways)

if (length(available_pathways) > 0) {
  p <- netVisual_bubble(cellchat, signaling = available_pathways, remove.isolate = FALSE)
  ggplot2::ggsave("pathways_of_interest_bubble.png", plot = p, width = 10, height = 8, dpi = 150)
  cat("Bubble plot saved for pathways:", paste(available_pathways, collapse = ", "), "\n")
} else {
  cat("None of the pathways of interest (PD-L1, CSF, CXCL, TGFb) were found as significant ",
      "in this dataset -- see all_significant_interactions.csv for what WAS found.\n")
}

cat("Module 6 complete for", patient_id, "\n")
