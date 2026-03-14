# ============================================================================
# TSENAT Appendix B: Comparing Flexible Parametric (GAM) and Rank-Based (K-W) Methods
# ============================================================================
# Purpose: Validates q-value × group interaction detection using two complementary
#          statistical approaches: GAM and Kruskal-Wallis
# Source: /home/nouser/galaxy/tools_source/TSENAT/vignettes/TSENAT_appendix_B.Rmd
# ============================================================================

# Setup: Load packages
suppressPackageStartupMessages({
    library(devtools)
    devtools::load_all(".")
    library(ggplot2)
    library(SummarizedExperiment)
    library(dplyr)
    library(gridExtra)
})

set.seed(42)

# Setup: Load and prepare data
# Load preprocessed dataset
data(readcounts)
readcounts <- as.matrix(salmon_dataset)
mode(readcounts) <- "numeric"

# Load TPM data
tpm_data <- as.matrix(salmon_tpm)
mode(tpm_data) <- "numeric"

# Subset for faster analysis (set to NULL for all genes)
n_genes_subset <- 100
if (!is.null(n_genes_subset) && nrow(readcounts) > n_genes_subset) {
  readcounts <- readcounts[1:n_genes_subset, ]
  tpm_data <- tpm_data[1:n_genes_subset, ]
  cat("NOTE: Analyzing", n_genes_subset, "genes for speed. Set n_genes_subset=NULL for all genes.\n\n")
}

# Load metadata
metadata_df <- read.table(
    system.file("extdata", "metadata.tsv", package = "TSENAT"),
    header = TRUE, sep = "\t"
)

gff3_dataset <- system.file("extdata", "annotation.gff3.gz", package = "TSENAT")

# Build SummarizedExperiment with counts assay
se <- build_se(readcounts, gff3_dataset, metadata = metadata_df)

# Add TPM as an assay, ensuring dimension alignment
# Select TPM data for genes in SE and samples in SE
tpm_subset <- tpm_data[rownames(se), colnames(se)]
assay(se, "tpm", withDimnames = FALSE) <- tpm_subset

se <- tryCatch(
  filter_se(se, stringency = "loose"),
  error = function(e) se
)

# Bayesian Pseudocount Regularization
# Estimate WLFC pseudocounts integrated workflow
result <- estimate_wlfc_pseudocounts(se, verbose = FALSE)
scalar_pseudocount <- result$scalar_pseudocount

# Calculate Tsallis entropy with Bayesian pseudocount regularization
qvec <- seq(0.1, 2, by = 0.05)
ts_se <- calculate_diversity(
  se, 
  q = qvec, 
  norm = TRUE, 
  pseudocount = scalar_pseudocount,
  metadata = metadata_df,
  shrinkage = "empirical_bayes"
)

# ============================================================================
# PART 1: Flexible Parametric Approach (GAM)
# ============================================================================
# GAM analysis is computationally intensive:
# - Per-gene fitting: Each gene gets a smooth GAM model fitted across all q-values.
# - Model complexity: Generalized Additive Models use iterative smoothing optimization.
# - Scale: With ~500 genes × 39 q-values = 19,500 model evaluations.
# - Multiple testing correction: Additional computation for Hochberg adjustment.
# Expected runtime: 2-5 minutes for 500 genes with 3 threads.

cat("\n=== PART 1: GAM Analysis ===\n")

# Flexible parametric approach (GAM)
lm_res <- calculate_lm_interaction(
    ts_se,
    method = "gam",
    paired = TRUE,
    multicorr = "hochberg",
    nthreads = 3,
    verbose = FALSE
)
