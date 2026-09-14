#!/usr/bin/env Rscript
# Run CellChat in spatial mode: infer ligand-receptor communication
# between predicted cell-type groups, constrained by real physical
# distance between spots (not just co-expression irrespective of
# location) using the human CellChatDB.
# Usage: run_cellchat.R <cellchat_input.rds> <min_cells>

args <- commandArgs(trailingOnly = TRUE)
input <- readRDS(args[1])
min_cells <- as.numeric(args[2])

suppressMessages(library(CellChat))

# CellChat's spatial mode needs the Visium spot diameter (in the same
# pixel units as the coordinates) to define what "physically adjacent"
# means -- taken from Space Ranger's own scale factors file, not guessed.
spot_diameter_fullres <- input$scalefactors$spot_diameter_fullres

cellchat <- createCellChat(
  object = input$data_matrix,
  meta = input$meta,
  group.by = "group",
  datatype = "spatial",
  coordinates = input$coords,
  spatial.factors = data.frame(ratio = 1, tol = spot_diameter_fullres / 2)
)

cellchat@DB <- CellChatDB.human

cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat, do.fast = FALSE)
cellchat <- identifyOverExpressedInteractions(cellchat)

cellchat <- computeCommunProb(
  cellchat,
  type = "truncatedMean", trim = 0.1,
  distance.use = TRUE,
  contact.knn.k = 6  # matches Module 5's k=6 (Visium hex grid: each spot has 6 real neighbors)
)

cellchat <- filterCommunication(cellchat, min.cells = min_cells)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

saveRDS(cellchat, "cellchat_obj.rds")
cat("CellChat run complete.\n")
