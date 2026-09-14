#!/usr/bin/env Rscript
# Test whether CNV score shows significant spatial autocorrelation --
# do CNV-high spots cluster together in physical space, or are they
# scattered randomly, Uses a permutation-based (Monte Carlo) Moran's I
# test, which makes no distributional assumptions: the observed
# statistic is compared against its distribution under many random
# relabelings of the same values across the same spatial positions.
# Usage: test_spatial_autocorrelation.R <spatial_neighbors.rds> <spatial_rds> <n_permutations>

args <- commandArgs(trailingOnly = TRUE)
neighbors <- readRDS(args[1])
seurat_obj <- readRDS(args[2])
n_sim <- as.numeric(args[3])

suppressMessages({ library(Seurat); library(spdep) })

cnv_score <- seurat_obj$cnv_score
names(cnv_score) <- colnames(seurat_obj)
cnv_score <- cnv_score[rownames(neighbors$coords)]  # match order to the neighbor graph

cat("Running permutation-based Moran's I test (", n_sim, "permutations) on CNV score...\n")
moran_result <- moran.mc(cnv_score, neighbors$listw, nsim = n_sim, zero.policy = TRUE)

cat("\n--- Moran's I test: CNV score spatial autocorrelation ---\n")
print(moran_result)

# For comparison/context: same test on the binary pathology malignant
# indicator, since we already know from Module 4 this should track the
# real tumor regions -- a sanity check that our neighbor graph and test
# setup are behaving sensibly, using a variable with an expected answer.
malignant_groups <- c("DCIS", "Invasive cancer + lymphocytes")
is_malignant <- as.numeric(seurat_obj$Classification %in% malignant_groups)
names(is_malignant) <- colnames(seurat_obj)
is_malignant <- is_malignant[rownames(neighbors$coords)]

moran_malignant <- moran.mc(is_malignant, neighbors$listw, nsim = n_sim, zero.policy = TRUE)
cat("\n--- Moran's I test: pathology malignant/normal label (sanity check) ---\n")
print(moran_malignant)

results <- data.frame(
  variable = c("cnv_score", "pathology_malignant_label"),
  morans_i = c(moran_result$statistic, moran_malignant$statistic),
  p_value  = c(moran_result$p.value, moran_malignant$p.value),
  n_permutations = n_sim
)
write.csv(results, "moran_test_results.csv", row.names = FALSE)

saveRDS(list(cnv_score = cnv_score, coords = neighbors$coords, listw = neighbors$listw),
        "moran_data_for_plotting.rds")

cat("\nSaved moran_test_results.csv\n")
