# Setup for TSENAT test suite
# ===========================
# This file is automatically sourced by testthat before running tests.
# It loads utility functions and factory functions used across test files.

# Source the test data factory functions
# testthat runs from tests/testthat/ directory, so tests-factory.R should be in current dir
factory_file <- "tests-factory.R"

if (file.exists(factory_file)) {
  source(factory_file, local = FALSE)
} else {
  # Fallback: try relative to package
  fallback_path <- file.path("../testthat", "tests-factory.R")
  if (file.exists(fallback_path)) {
    source(fallback_path, local = FALSE)
  } else {
    warning("Could not find tests-factory.R. Test factory functions will not be available.")
  }
}

# Silence package startup messages during tests
# Note: testthat and SummarizedExperiment are assumed available as test dependencies
# Do not use library() calls - they violate CRAN/Bioconductor guidelines

# ============================================================================
# VIGNETTE DATA LOADING - Load exactly as in TSENAT.Rmd (Lines 310-320)
# ============================================================================
# Load built-in vignette data (readcounts, tpm, effective_length all load together)
data(readcounts, package = "TSENAT", envir = environment())
readcounts <- as.matrix(readcounts)
mode(readcounts) <- "numeric"

# tpm and effective_length are auto-loaded with readcounts data() call
tpm <- as.matrix(tpm)
mode(tpm) <- "numeric"
effective_length <- as.numeric(effective_length)

# Load metadata
metadata_df <- read.table(
  system.file("extdata", "metadata.tsv", package = "TSENAT"),
  header = TRUE, sep = "\t"
)

# Get the GFF3.gz annotation file path
gff3_dataset <- system.file("extdata", "annotation.gff3.gz", package = "TSENAT")

# ============================================================================
# NOTE: TESTTHAT FALSE POSITIVE WARNINGS
# ============================================================================
# testthat's code instrumentation engine produces "stack imbalance" warnings
# for code with nested function calls (like TSENAT::: calls with rep(c(...)))
# These are false positives from testthat's AST rewriting and don't indicate
# actual syntax errors. All tests pass despite these warnings.
# These warnings cannot be suppressed within R code as they're generated
# by the C-level bytecode compiler.
