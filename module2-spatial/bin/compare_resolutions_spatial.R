#!/usr/bin/env Rscript
# Fan-in: gather all parallel per-resolution results, pick the best by ARI.

args <- commandArgs(trailingOnly = TRUE)
sep_idx <- which(args == "--")
ari_files <- args[1:(sep_idx - 1)]
rds_files <- args[(sep_idx + 1):length(args)]

results <- do.call(rbind, lapply(ari_files, read.csv))
results <- results[order(-results$ari), ]

cat("--- Spatial resolution comparison (all parallel branches) ---\n")
print(results)
write.csv(results, "resolution_comparison_spatial.csv", row.names = FALSE)

best_resolution <- results$resolution[1]
res_tag <- gsub("\\.", "p", sprintf("%.2f", best_resolution))
best_rds <- rds_files[grepl(paste0("res", res_tag, "\\.rds$"), rds_files)]

cat("\nBest resolution:", best_resolution, "(ARI =", results$ari[1], ")\n")
file.copy(best_rds, "best_clustered_spatial.rds", overwrite = TRUE)
