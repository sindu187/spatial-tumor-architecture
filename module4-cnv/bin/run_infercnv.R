#!/usr/bin/env Rscript
# Run inferCNV, using pathology-annotated normal-like spots as the
# reference baseline. HMM-based CNV state calling is disabled -- we only
# need the underlying denoised, smoothed relative expression signal for
# our own quantitative validation, not formal per-region state calls,
# and HMM adds substantial runtime.
# Usage: run_infercnv.R <counts_matrix> <annotation_file> <gene_order_file> <ref_groups_csv> <cutoff> <num_threads>

args <- commandArgs(trailingOnly = TRUE)
counts_file <- args[1]; annotation_file <- args[2]; gene_order_file <- args[3]
ref_groups <- strsplit(args[4], ",")[[1]]
cutoff <- as.numeric(args[5])
num_threads <- as.numeric(args[6])

suppressMessages(library(infercnv))

cat("Reference (normal-like) groups:", paste(ref_groups, collapse = " | "), "\n")

infercnv_obj <- CreateInfercnvObject(
  raw_counts_matrix = counts_file,
  annotations_file  = annotation_file,
  gene_order_file   = gene_order_file,
  ref_group_names   = ref_groups
)

infercnv_obj <- infercnv::run(
  infercnv_obj,
  cutoff = cutoff,
  out_dir = "infercnv_output",
  cluster_by_groups = TRUE,
  denoise = TRUE,
  HMM = FALSE,
  num_threads = num_threads
)

saveRDS(infercnv_obj, "final_infercnv_obj.rds")
cat("\ninferCNV run complete.\n")
