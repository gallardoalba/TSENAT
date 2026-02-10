# Load package and testthat to ensure test helpers like context() are available
if (!requireNamespace("pkgload", quietly = TRUE)) install.packages("pkgload", repos = 'https://cloud.r-project.org')
if (!requireNamespace("testthat", quietly = TRUE)) install.packages("testthat", repos = 'https://cloud.r-project.org')
suppressPackageStartupMessages(pkgload::load_all())

# Collect test files and compute file coverage
test_files <- list.files("tests/testthat", pattern = "\\.R$", full.names = TRUE)
fc <- covr::file_coverage("R/calc_lm_helpers.R", test_files)
cat("calc_lm_helpers.R coverage:", covr::percent_coverage(fc), "\n\n")
# print a brief summary of uncovered lines
cov_lines <- as.data.frame(fc)
uncovered <- cov_lines[cov_lines$value == 0, ]
if (nrow(uncovered) == 0) {
  cat("No uncovered lines detected\n")
} else {
  cat("Uncovered lines (first 200):\n")
  print(head(uncovered, 200))
}