#!/usr/bin/env Rscript
# Investigate why the anchor method assigns zero score to CD4+ T cells
# everywhere. Hypothesis: CD4+ signal is being "absorbed" by competing
# T-cell-related categories (CD8+, Cycling T-cells) during anchor-based
# integration, rather than being genuinely absent.

anchor <- read.csv("results/01_anchor_transfer/anchor_predictions.csv", row.names = 1)
marker <- read.csv("results/02_marker_scoring/marker_scores.csv", row.names = 1)

common_spots <- intersect(rownames(anchor), rownames(marker))
anchor <- anchor[common_spots, ]
marker <- marker[common_spots, ]

# Spots where the INDEPENDENT marker method thinks CD4+ signal is
# strongest -- these are the spots where, if the hypothesis is right,
# we'd expect the anchor method to be suppressing real CD4+ signal.
top_cd4_spots <- rownames(marker)[order(-marker$T.cells.CD4.)][1:20]

cat("--- Anchor method's predicted.id for the TOP 20 marker-CD4-scoring spots ---\n")
print(table(anchor[top_cd4_spots, "predicted.id"]))

cat("\n--- Anchor prediction scores for T-cell-related categories, ")
cat("averaged over these top-CD4 (by marker score) spots ---\n")
tcell_cols <- c("prediction.score.T.cells.CD4.", "prediction.score.T.cells.CD8.", "prediction.score.Cycling.T.cells")
print(colMeans(anchor[top_cd4_spots, tcell_cols]))

cat("\n--- Correlation: marker-based CD4 score vs. anchor-based CD8/Cycling scores ---\n")
cat("marker CD4 vs anchor CD8:     ", cor(marker$T.cells.CD4., anchor$prediction.score.T.cells.CD8., method = "spearman"), "\n")
cat("marker CD4 vs anchor Cycling: ", cor(marker$T.cells.CD4., anchor$prediction.score.Cycling.T.cells, method = "spearman"), "\n")
