# Setup file for testthat tests
# This file is automatically loaded before any tests run

# Ensure all necessary packages are available
.onLoad <- function(libname, pkgname) {
  # testthat should already be running
  # Just ensure our functions are sourced
  if (!exists(".tsenat_apply_friedman_test", mode = "function")) {
    source("R/rank_based_methods.R", local = TRUE)
  }
}

# Set seed for reproducibility across test runs
set.seed(42)
