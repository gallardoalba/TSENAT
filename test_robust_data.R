library(SummarizedExperiment)
set.seed(42)
# Create robust test data: each cell has 10-100 counts (avoids sparsity filtering)
counts_matrix <- matrix(
  sample(10:100, 2000, replace = TRUE), nrow = 50, ncol = 40,
  dimnames = list(paste0("ENST", 1:50), paste0("Sample", 1:40))
)
se <- SummarizedExperiment(
  assays = list(counts = counts_matrix),
  colData = data.frame(condition = rep(c("A", "B"), 20))
)
analysis <- TSENATAnalysis(se)
analysis <- calculate_diversity_s4(analysis, q = 1.0)
cat("Success! Diversity results:\n")
cat("Number of genes:", nrow(analysis@diversity_results[[1]]), "\n")
cat("Result structure:\n")
str(analysis@diversity_results[[1]], max.level = 2)
