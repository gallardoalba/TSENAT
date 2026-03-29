# Setup for TSENAT test suite
# ===========================
# This file is automatically sourced by testthat before running tests.
# It loads utility functions and factory functions used across test files.

# Source the test data factory functions
# Using relative path from this file's location
source("tests-factory.R", local = FALSE)

# Silence package startup messages during tests
suppressPackageStartupMessages({
  library(testthat)
  library(SummarizedExperiment)
})
