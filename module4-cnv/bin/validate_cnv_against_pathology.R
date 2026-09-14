#!/usr/bin/env Rscript
# Quantitatively test whether inferred CNV signal actually discriminates
# malignant from normal spots 
# Usage: validate_cnv_against_pathology.R <final_infercnv_obj.rds> <annotation_file> <ref_groups_csv> <malignant_groups_csv>

args <- commandArgs(trailingOnly = TRUE)
infercnv_obj <- readRDS(args[1])
annotation_file <- args[2]
ref_groups <- strsplit(args[3], ",")[[1]]
malignant_groups <- strsplit(args[4], ",")[[1]]

suppressMessages(library(pROC))

# inferCNV's final expr.data is the denoised, smoothed relative
# expression matrix, centered around 1.0 = neutral copy number. A
# per-spot CNV score summarizes how far, on average, that spot's genome-
# wide signal deviates from the neutral baseline.
expr <- infercnv_obj@expr.data
cnv_score <- colMeans((expr - 1)^2)

ann <- read.table(annotation_file, sep = "\t", stringsAsFactors = FALSE)
colnames(ann) <- c("spot", "group")
rownames(ann) <- ann$spot

df <- data.frame(spot = names(cnv_score), cnv_score = cnv_score, group = ann[names(cnv_score), "group"])
df$category <- ifelse(df$group %in% malignant_groups, "Malignant",
                ifelse(df$group %in% ref_groups, "Normal-like", "Other"))

sub <- df[df$category %in% c("Malignant", "Normal-like"), ]
cat("Malignant spots:", sum(sub$category == "Malignant"), "\n")
cat("Normal-like spots:", sum(sub$category == "Normal-like"), "\n\n")

wilcox_result <- wilcox.test(cnv_score ~ category, data = sub)
cat("--- Wilcoxon rank-sum test: CNV score, Malignant vs. Normal-like ---\n")
print(wilcox_result)

roc_obj <- roc(response = sub$category, predictor = sub$cnv_score,
                levels = c("Normal-like", "Malignant"), direction = "<", quiet = TRUE)
cat("\nAUC (does CNV score discriminate malignant from normal spots?):", round(as.numeric(auc(roc_obj)), 3), "\n")

write.csv(df, "cnv_scores_per_spot.csv", row.names = FALSE)

summary_stats <- data.frame(
  wilcox_p_value = wilcox_result$p.value,
  auc = as.numeric(auc(roc_obj)),
  n_malignant = sum(sub$category == "Malignant"),
  n_normal = sum(sub$category == "Normal-like")
)
write.csv(summary_stats, "cnv_validation_summary.csv", row.names = FALSE)

png("cnv_score_boxplot.png", width = 800, height = 600)
boxplot(cnv_score ~ category, data = sub,
        ylab = "CNV score (mean squared deviation from neutral)",
        main = "Inferred CNV score by pathology category")
dev.off()

cat("\nSaved cnv_scores_per_spot.csv, cnv_validation_summary.csv, cnv_score_boxplot.png\n")
