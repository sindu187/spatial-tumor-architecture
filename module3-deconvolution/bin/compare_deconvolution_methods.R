#!/usr/bin/env Rscript
# Compare the two independent deconvolution methods: for each cell type
# present in both, correlate Method A's prediction score against Method
# B's marker score across all spots. High correlation = both methods
# independently agree this cell type is more/less present in the same
# spots -- a real validation, since the two methods share no machinery.
#
# Usage: compare_deconvolution_methods.R <anchor_predictions.csv> <marker_scores.csv>

args <- commandArgs(trailingOnly = TRUE)
anchor <- read.csv(args[1], row.names = 1)
marker <- read.csv(args[2], row.names = 1)

# Match cell type names: anchor predictions are columns like
# "prediction.score.T.cells.CD4." (Seurat sanitizes names); marker score
# columns are the raw assigned_label strings. Match by best substring fit.
anchor_score_cols <- grep("^prediction\\.score\\.", names(anchor), value = TRUE)
anchor_score_cols <- setdiff(anchor_score_cols, "prediction.score.max")

cat("Anchor method cell types found:", length(anchor_score_cols), "\n")
cat("Marker method cell types found:", ncol(marker), "\n")

common_spots <- intersect(rownames(anchor), rownames(marker))
cat("Common spots for comparison:", length(common_spots), "\n\n")

results <- list()
for (marker_type in names(marker)) {
  # loose match: strip non-alphanumeric chars from both sides and compare
  clean_marker <- gsub("[^A-Za-z0-9]", "", marker_type)
  matches <- anchor_score_cols[gsub("[^A-Za-z0-9]", "", gsub("prediction.score.", "", anchor_score_cols)) == clean_marker]

  if (length(matches) == 1) {
    a_vals <- anchor[common_spots, matches]
    m_vals <- marker[common_spots, marker_type]
    correlation <- cor(a_vals, m_vals, method = "spearman", use = "complete.obs")
    results[[marker_type]] <- correlation
  }
}

comparison_df <- data.frame(
  cell_type = names(results),
  spearman_correlation = unlist(results)
)
comparison_df <- comparison_df[order(-comparison_df$spearman_correlation), ]

write.csv(comparison_df, "method_agreement.csv", row.names = FALSE)
cat("--- Method A vs. Method B agreement (Spearman correlation per cell type) ---\n")
print(comparison_df)
cat("\nMedian correlation across matched cell types:", round(median(comparison_df$spearman_correlation, na.rm = TRUE), 3), "\n")
