# Test S4 Class Validation for TSENATAnalysis
# Tests parameter validation in TSENATAnalysis constructor and setConfig method
# This ensures invalid configuration parameters are caught early (fail-fast)

context("TSENATAnalysis constructor parameter validation")

test_that("Constructor rejects invalid condition_col", {
    skip_if_not_installed("SummarizedExperiment")
    library(SummarizedExperiment)
    
    se <- SummarizedExperiment(
        assays = list(counts = matrix(1:100, nrow = 10, ncol = 10)),
        colData = data.frame(
            sample_id = sprintf("S%d", 1:10),
            condition = rep(c("A", "B"), 5),
            row.names = sprintf("S%d", 1:10)
        )
    )
    
    expect_error(
        TSENATAnalysis(se, config = list(condition_col = "sample_type")),
        "Invalid condition_col.*not found"
    )
})

test_that("Constructor rejects invalid subject_col when paired", {
    skip_if_not_installed("SummarizedExperiment")
    library(SummarizedExperiment)
    
    se <- SummarizedExperiment(
        assays = list(counts = matrix(1:100, nrow = 10, ncol = 10)),
        colData = data.frame(
            sample_id = sprintf("S%d", 1:10),
            condition = rep(c("A", "B"), 5),
            row.names = sprintf("S%d", 1:10)
        )
    )
    
    expect_error(
        TSENATAnalysis(se, config = list(
            condition_col = "condition",
            paired = TRUE,
            subject_col = "nonexistent_subject"
        )),
        "Invalid subject_col.*not found"
    )
})

test_that("Constructor accepts valid parameters", {
    skip_if_not_installed("SummarizedExperiment")
    library(SummarizedExperiment)
    
    se <- SummarizedExperiment(
        assays = list(counts = matrix(1:100, nrow = 10, ncol = 10)),
        colData = data.frame(
            sample_id = sprintf("S%d", 1:10),
            condition = rep(c("A", "B"), 5),
            row.names = sprintf("S%d", 1:10)
        )
    )
    
    expect_silent(
        obj <- TSENATAnalysis(se, config = list(condition_col = "condition"))
    )
    expect_s4_class(obj, "TSENATAnalysis")
})

# ============================================================================
# setValidity() comprehensive checks
# ============================================================================

context("TSENATAnalysis setValidity() comprehensive checks")

test_that("validObject() catches invalid condition_col", {
    skip_if_not_installed("SummarizedExperiment")
    library(SummarizedExperiment)
    
    se <- SummarizedExperiment(
        assays = list(counts = matrix(1:100, nrow = 10, ncol = 10)),
        colData = data.frame(
            sample_id = sprintf("S%d", 1:10),
            condition = rep(c("A", "B"), 5)
        )
    )
    
    obj <- new("TSENATAnalysis", se = se, config = list(condition_col = "invalid"))
    
    expect_error(
        validObject(obj),
        "condition_col.*not found"
    )
})

test_that("validObject() catches invalid subject_col", {
    skip_if_not_installed("SummarizedExperiment")
    library(SummarizedExperiment)
    
    se <- SummarizedExperiment(
        assays = list(counts = matrix(1:100, nrow = 10, ncol = 10)),
        colData = data.frame(
            sample_id = sprintf("S%d", 1:10),
            condition = rep(c("A", "B"), 5)
        )
    )
    
    obj <- new("TSENATAnalysis", se = se, config = list(
        condition_col = "condition",
        paired = TRUE,
        subject_col = "invalid_subject"
    ))
    
    expect_error(
        validObject(obj),
        "subject_col.*not found"
    )
})

test_that("validObject() allows NULL or missing condition_col", {
    skip_if_not_installed("SummarizedExperiment")
    library(SummarizedExperiment)
    
    se <- SummarizedExperiment(
        assays = list(counts = matrix(1:100, nrow = 10, ncol = 10)),
        colData = data.frame(sample_id = sprintf("S%d", 1:10))
    )
    
    # NULL config
    obj1 <- new("TSENATAnalysis", se = se, config = list())
    expect_true(validObject(obj1))
    
    # NULL condition_col
    obj2 <- new("TSENATAnalysis", se = se, config = list(condition_col = NULL))
    expect_true(validObject(obj2))
})

test_that("validObject() ignores subject_col when paired=FALSE", {
    skip_if_not_installed("SummarizedExperiment")
    library(SummarizedExperiment)
    
    se <- SummarizedExperiment(
        assays = list(counts = matrix(1:100, nrow = 10, ncol = 10)),
        colData = data.frame(
            sample_id = sprintf("S%d", 1:10),
            condition = rep(c("A", "B"), 5)
        )
    )
    
    # subject_col doesn't exist but paired=FALSE, so should be OK
    obj <- new("TSENATAnalysis", se = se, config = list(
        condition_col = "condition",
        paired = FALSE,
        subject_col = "nonexistent"
    ))
    
    expect_true(validObject(obj))
})

# ============================================================================
# setConfig() method with validation
# ============================================================================

context("setConfig() method with validation")

test_that("setConfig() validates config parameter", {
    skip_if_not_installed("SummarizedExperiment")
    library(SummarizedExperiment)
    
    se <- SummarizedExperiment(
        assays = list(counts = matrix(1:100, nrow = 10, ncol = 10)),
        colData = data.frame(
            sample_id = sprintf("S%d", 1:10),
            condition = rep(c("A", "B"), 5)
        )
    )
    
    obj <- TSENATAnalysis(se)
    
    expect_error(
        setConfig(obj, list(condition_col = "invalid")),
        "condition_col.*not found"
    )
})

test_that("setConfig() accepts valid config", {
    skip_if_not_installed("SummarizedExperiment")
    library(SummarizedExperiment)
    
    se <- SummarizedExperiment(
        assays = list(counts = matrix(1:100, nrow = 10, ncol = 10)),
        colData = data.frame(
            sample_id = sprintf("S%d", 1:10),
            condition = rep(c("A", "B"), 5)
        )
    )
    
    obj <- TSENATAnalysis(se)
    
    expect_silent(
        obj_updated <- setConfig(obj, list(condition_col = "condition"))
    )
    expect_equal(obj_updated@config$condition_col, "condition")
})

# ============================================================================
# Integration tests with real vignette data
# ============================================================================

context("Integration tests with real metadata")

test_that("Real vignette metadata works with correct config", {
    skip_if_not_installed("SummarizedExperiment")
    library(SummarizedExperiment)
    
    # Use real metadata from vignette data
    metadata_file <- system.file("extdata", "metadata.tsv", package = "TSENAT")
    
    if (file.exists(metadata_file)) {
        metadata_df <- read.table(metadata_file, header = TRUE, sep = "\t")
        
        # Create SE with real metadata
        se <- SummarizedExperiment(
            assays = list(counts = matrix(rpois(1000, 10), nrow = 100, ncol = 10)),
            colData = metadata_df
        )
        
        # Valid: "condition" is real column
        expect_silent(
            obj <- TSENATAnalysis(se, config = list(condition_col = "condition"))
        )
        expect_s4_class(obj, "TSENATAnalysis")
        
        # Invalid: "sample_type" doesn't exist (this was the original bug)
        expect_error(
            TSENATAnalysis(se, config = list(condition_col = "sample_type")),
            "Invalid condition_col.*not found"
        )
    }
})
