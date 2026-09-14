#!/usr/bin/env Rscript
# Subset one patient's cells from the combined Wu et al. matrix.
# Usage: subset_patient.R <mtx> <genes> <barcodes> <metadata> <patient_id>

args <- commandArgs(trailingOnly = TRUE)
mtx_path <- args[1]; genes_path <- args[2]; barcodes_path <- args[3]
metadata_path <- args[4]; patient_id <- args[5]

suppressMessages({ library(Matrix); library(Seurat) })

cat("Loading full matrix...\n")
mat <- readMM(mtx_path)
genes <- read.csv(genes_path, header = FALSE)$V1
barcodes <- read.csv(barcodes_path, header = FALSE)$V1
metadata <- read.csv(metadata_path, row.names = 1)
rownames(mat) <- genes
colnames(mat) <- barcodes

cid_cells <- rownames(metadata)[metadata$orig.ident == patient_id]
cat(sprintf("%s: %d cells found\n", patient_id, length(cid_cells)))
if (length(cid_cells) == 0) stop("No cells found for patient_id = ", patient_id)

mat_subset <- mat[, colnames(mat) %in% cid_cells]
meta_subset <- metadata[colnames(mat_subset), ]

seurat_obj <- CreateSeuratObject(counts = mat_subset, meta.data = meta_subset, project = patient_id)
seurat_obj[["percent.mt.recomputed"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

saveRDS(seurat_obj, "raw_unfiltered.rds")
cat("Saved raw_unfiltered.rds (", ncol(seurat_obj), "cells )\n")
