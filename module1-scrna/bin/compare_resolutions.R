#!/usr/bin/env Rscript
# Fan-in step: gather ALL the parallel per-resolution results, build a
# comparison table, and select the object corresponding to the highest ARI.
# Usage: compare_resolutions.R <ari_csv_1> <ari_csv_2> ... -- <rds_1> <rds_2> ...
# (files passed via Nextflow's .collect() -- see main.nf)

args <- commandArgs(trailingOnly = TRUE)
sep_idx <- which(args == "--")
ari_files <- args[1:(sep_idx - 1)]
rds_files <- args[(sep_idx + 1):length(args)]

results <- do.call(rbind, lapply(ari_files, read.csv))
results <- results[order(-results$ari), ]

cat("--- Resolution comparison (all parallel branches) ---\n")
print(results)

write.csv(results, "resolution_comparison.csv", row.names = FALSE)

best_resolution <- results$resolution[1]
res_tag <- gsub("\\.", "p", sprintf("%.2f", best_resolution))
best_rds <- rds_files[grepl(paste0("res", res_tag, "\\.rds$"), rds_files)]

cat("\nBest resolution:", best_resolution, "(highest ARI =", results$ari[1], ")\n")

file.copy(best_rds, "best_clustered.rds", overwrite = TRUE)
cat("Copied", best_rds, "-> best_clustered.rds\n")
