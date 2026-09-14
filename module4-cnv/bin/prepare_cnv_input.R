#!/usr/bin/env Rscript
# Build inferCNV's three required inputs: raw counts matrix (genes x
# spots, ordered to match the gene order file), the gene order file
# restricted to genes actually present, and an annotation file mapping
# each spot to its pathology category (excluding spots with no usable
# label).
# Usage: prepare_cnv_input.R <spatial_rds> <gene_order_file>

args <- commandArgs(trailingOnly = TRUE)
seurat_obj <- readRDS(args[1])
gene_order <- read.table(args[2], sep = "\t", stringsAsFactors = FALSE)
colnames(gene_order) <- c("gene", "chr", "start", "end")

suppressMessages(library(Seurat))

meta <- seurat_obj@meta.data
meta$Classification[meta$Classification %in% c("Artefact", "Classification")] <- NA
valid_spots <- rownames(meta)[!is.na(meta$Classification)]

counts <- GetAssayData(seurat_obj, assay = "Spatial", layer = "counts")
counts <- counts[, colnames(counts) %in% valid_spots]

common_genes <- intersect(rownames(counts), gene_order$gene)
counts <- counts[common_genes, ]
gene_order <- gene_order[gene_order$gene %in% common_genes, ]

standard_chrs <- paste0("chr", c(1:22, "X", "Y"))
gene_order$chr <- factor(gene_order$chr, levels = standard_chrs)
gene_order <- gene_order[order(gene_order$chr, gene_order$start), ]
counts <- counts[gene_order$gene, ]

annotation <- data.frame(spot = colnames(counts), group = meta[colnames(counts), "Classification"])

write.table(as.matrix(counts), "counts_matrix.txt", sep = "\t", quote = FALSE, col.names = NA)
write.table(gene_order, "gene_order_filtered.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(annotation, "annotation_file.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

cat("Prepared CNV input:", nrow(counts), "genes x", ncol(counts), "spots\n")
cat(nrow(annotation), "spots with valid pathology annotation\n")
cat("\n--- Spots per pathology category ---\n")
print(table(annotation$group))
