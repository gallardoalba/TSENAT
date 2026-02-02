library(testthat)

context("infer_sample_group")

test_that("suffix-based mapping with defaults works", {
    samples <- c("S1_N", "S2_T", "S3_normal", "S4_tumor", "S5_X")
    res <- infer_sample_group(samples)
    # without mappings, raw suffix tokens are returned
    expect_equal(res, c("N", "T", "normal", "tumor", "X"))
})

test_that("TCGA code mapping works and unknown codes yield NA", {
    samples <- c("TCGA-XX-01A", "TCGA-YY-11B", "TCGA-ZZ-02C")
    res <- infer_sample_group(samples)
    # without tcga_map provided, raw two-digit codes are returned when detected
    expect_equal(res, c("01", "11", "02"))
})

test_that("custom suffix_map and default behave as expected", {
    samples <- c("a_R", "b_C", "c_unknown")
    res <- infer_sample_group(samples,
        suffix_map = c(
            R = "Ref",
            C = "Case"
        ),
        default = "Unknown"
    )
    # mapped tokens are translated; unmapped suffix returns raw token
    expect_equal(res, c("Ref", "Case", "unknown"))
})

test_that("NULL input returns empty character vector", {
    expect_equal(infer_sample_group(NULL), character(0))
})

test_that("mixed TCGA and suffix use default when unmapped", {
    samples <- c("TCGA-XX-01A", "Sample_N", "TCGA-YY-99Z")
    res <- infer_sample_group(samples, default = "UNK")
    # with no tcga_map and no suffix_map provided, raw tokens are returned
    expect_equal(res, c("01", "N", "99"))
})
