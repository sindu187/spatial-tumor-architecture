#!/usr/bin/env Rscript
# Load Visium spatial data (Space Ranger output) and attach pathologist
# annotations from the study's own metadata as validation.
# Usage: load_spatial.R <spaceranger_dir> <metadata_csv> <patient_id>

args <- commandArgs(trailingOnly = TRUE)
spaceranger_dir <- args[1]; metadata_path <- args[2]; patient_id <- args[3]

suppressMessages(library(Seurat))

cat("Loading Visium data from:", spaceranger_dir, "\n")


# (filtered_feature_bc_matrix.h5).load counts and the tissue image separately instead.
counts <- Read10X(
  data.dir = file.path(spaceranger_dir, "filtered_feature_bc_matrix"),
  gene.column = 1  # this release has only gene symbols (1 column), not the
                    # standard 10x 3-column format (Ensembl ID, symbol, type)
)
seurat_obj <- CreateSeuratObject(counts = counts, assay = "Spatial", project = patient_id)

image <- Read10X_Image(image.dir = file.path(spaceranger_dir, "spatial"))
image <- image[Cells(seurat_obj)]
DefaultAssay(image) <- "Spatial"
seurat_obj[[patient_id]] <- image

seurat_obj[["percent.mito"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

cat("Loaded", ncol(seurat_obj), "spots x", nrow(seurat_obj), "genes\n")

# --- Attach pathology ---
pathology <- read.csv(metadata_path, row.names = 1)
common_barcodes <- intersect(colnames(seurat_obj), rownames(pathology))
cat("Spots with pathology annotation available:", length(common_barcodes), "/", ncol(seurat_obj), "\n")

seurat_obj$Classification <- NA
seurat_obj$Classification[common_barcodes] <- pathology[common_barcodes, "Classification"]

cat("\n--- Pathology classification breakdown ---\n")
print(table(seurat_obj$Classification, useNA = "ifany"))

saveRDS(seurat_obj, "raw_spatial.rds")
cat("\nSaved raw_spatial.rds\n")
