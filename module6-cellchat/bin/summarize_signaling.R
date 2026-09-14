#!/usr/bin/env Rscript
# Summarize CellChat's significant ligand-receptor results, with
# particular attention to axes relevant to this TNBC tumor's biology:
# PD-L1/PD-1 (immune checkpoint), CSF1/CSF1R (macrophage recruitment),
# CXCL12/CXCR4 (tumor-stroma), TGFB1/TGFBR1 (CAF-tumor crosstalk).
# Usage: summarize_signaling.R <cellchat_obj.rds>

args <- commandArgs(trailingOnly = TRUE)
cellchat <- readRDS(args[1])

suppressMessages(library(CellChat))

all_results <- subsetCommunication(cellchat)
write.csv(all_results, "all_significant_interactions.csv", row.names = FALSE)
cat("Total significant ligand-receptor interactions found:", nrow(all_results), "\n")

pathways_of_interest <- c("PD-L1", "CSF", "CXCL", "TGFb")
highlighted <- all_results[grepl(paste(pathways_of_interest, collapse = "|"),
                                    all_results$pathway_name, ignore.case = TRUE), ]

write.csv(highlighted, "highlighted_interactions.csv", row.names = FALSE)
cat("\n--- Interactions matching pathways of interest (checkpoint, macrophage, stroma) ---\n")
print(highlighted)

pathway_ranking <- aggregate(prob ~ pathway_name, data = all_results, FUN = sum)
pathway_ranking <- pathway_ranking[order(-pathway_ranking$prob), ]
write.csv(pathway_ranking, "pathway_ranking.csv", row.names = FALSE)

cat("\n--- Top 10 pathways by total communication probability ---\n")
print(head(pathway_ranking, 10))
