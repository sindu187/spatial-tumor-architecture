#!/usr/bin/env Rscript
# Build the gene order file inferCNV needs (gene, chr, start, end),
# using official Bioconductor annotation packages rather than an
# unofficial downloaded file -- fully reproducible from source packages.
# Usage: build_gene_order_file.R <spatial_rds>

args <- commandArgs(trailingOnly = TRUE)
seurat_obj <- readRDS(args[1])

suppressMessages({
  library(Seurat)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  library(GenomicFeatures)
})

genes_in_data <- rownames(seurat_obj)
cat("Genes in spatial data:", length(genes_in_data), "\n")

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
gene_ranges <- genes(txdb)  # GRanges keyed by Entrez ID

symbol_to_entrez <- AnnotationDbi::select(
  org.Hs.eg.db, keys = genes_in_data, keytype = "SYMBOL", columns = "ENTREZID"
)
symbol_to_entrez <- symbol_to_entrez[!is.na(symbol_to_entrez$ENTREZID) & !duplicated(symbol_to_entrez$SYMBOL), ]

matched <- gene_ranges[names(gene_ranges) %in% symbol_to_entrez$ENTREZID]
df <- data.frame(
  entrez = names(matched),
  chr    = as.character(GenomicRanges::seqnames(matched)),
  start  = GenomicRanges::start(matched),
  end    = GenomicRanges::end(matched)
)
df <- merge(df, symbol_to_entrez, by.x = "entrez", by.y = "ENTREZID")
df <- df[!duplicated(df$SYMBOL), ]

standard_chrs <- paste0("chr", c(1:22, "X", "Y"))
df <- df[df$chr %in% standard_chrs, ]
df$chr <- factor(df$chr, levels = standard_chrs)
df <- df[order(df$chr, df$start), ]

out <- df[, c("SYMBOL", "chr", "start", "end")]
write.table(out, "gene_order.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

cat("Gene order file built for", nrow(out), "of", length(genes_in_data), "genes\n")
cat("(remainder lacked a mappable genomic position or were on non-standard contigs)\n")
