#!/usr/bin/env Rscript
# Test whether lumping DCIS and Invasive cancer into one Malignant
# category diluted the discrimination signal -- check each separately
# against the normal-like baseline.

suppressMessages(library(pROC))

df <- read.csv("results/04_validate/cnv_scores_per_spot.csv")

ref_groups <- c("Normal + stroma + lymphocytes", "Stroma", "Stroma + adipose tissue")

for (malignant_type in c("DCIS", "Invasive cancer + lymphocytes")) {
  sub <- df[df$group %in% c(malignant_type, ref_groups), ]
  sub$category <- ifelse(sub$group == malignant_type, "Malignant", "Normal-like")

  wilcox_result <- wilcox.test(cnv_score ~ category, data = sub)
  roc_obj <- roc(response = sub$category, predictor = sub$cnv_score,
                  levels = c("Normal-like", "Malignant"), direction = "<", quiet = TRUE)

  cat("===", malignant_type, "vs. Normal-like ===\n")
  cat("n malignant:", sum(sub$category == "Malignant"),
      " n normal:", sum(sub$category == "Normal-like"), "\n")
  cat("Wilcoxon p-value:", format(wilcox_result$p.value, scientific = TRUE), "\n")
  cat("AUC:", round(as.numeric(auc(roc_obj)), 3), "\n\n")
}
