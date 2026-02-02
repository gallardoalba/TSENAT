# Minimal runner to validate Tsallis differential analysis

library(SummarizedExperiment)

## load package (R CMD check runs tests with package built/installed)
if (interactive()) {
  # during local interactive runs, source R files for quick testing
  r_files <- list.files("R", pattern = "\\.R$", full.names = TRUE)
  for (f in r_files) source(f)
} else {
  library(abyssEdge)
}

# create synthetic transcript-level read counts
set.seed(42)
mat <- matrix(rpois(60, lambda = 10), ncol = 6)
colnames(mat) <- paste0("S", 1:6)
# gene assignment must match nrow(mat) == length(genes)
genes <- c(rep("G1", 3), rep("G2", 3), rep("G3", 2), rep("G4", 2))

# calculate Tsallis diversity (q = 2)
se <- calculate_diversity(mat, genes, q = 2)

# convert to data.frame expected by calculate_difference
library(SummarizedExperiment)
div_df <- as.data.frame(SummarizedExperiment::assay(se))
div_df <- cbind(genes = rowData(se)$genes, div_df)

# sample categories for difference
samples <- c(rep("Healthy", 3), rep("Pathogenic", 3))

# run difference analysis
res <- calculate_difference(div_df, samples, control = "Healthy", method = "mean", test = "wilcoxon")

print(head(res))
cat("Rows:", nrow(res), "\n")
