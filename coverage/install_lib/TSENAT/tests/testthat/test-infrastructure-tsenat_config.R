context("Orchestration: Configuration Builder")

test_that("TSENAT_config creates config with defaults", {
  cfg <- TSENAT_config()

  expect_s3_class(cfg, "TSENATConfig")
  expect_true(is.list(cfg))
  expect_true("q" %in% names(cfg))
  expect_equal(cfg$fdr_threshold, 0.05)
  expect_equal(cfg$p_threshold, 0.05)
})

test_that("TSENAT_config accepts custom q values", {
  cfg <- TSENAT_config(q = c(0.5, 1.0, 1.5))

  expect_equal(cfg$q, c(0.5, 1.0, 1.5))
})

test_that("TSENAT_config accepts custom thresholds", {
  cfg <- TSENAT_config(p_threshold = 0.01, fdr_threshold = 0.001)

  expect_equal(cfg$p_threshold, 0.01)
  expect_equal(cfg$fdr_threshold, 0.001)
})

test_that("TSENAT_config accepts formula", {
  form <- ~ treatment + batch
  cfg <- TSENAT_config(formula = form)

  expect_equal(cfg$formula, form)
})

test_that("TSENAT_config accepts additional parameters", {
  cfg <- TSENAT_config(custom_param = "value")

  expect_equal(cfg$custom_param, "value")
})

test_that("TSENAT_config has TSENATConfig class", {
  cfg <- TSENAT_config()

  expect_true("TSENATConfig" %in% class(cfg))
  expect_true(is.list(cfg))
})

# ============================================================================
# TEST SUITE: Bioconductor c2e8214 - TSENAT_config Parameter Validation
# ============================================================================
# Tests for match.arg() parameter validation in TSENAT_config function

test_that("TSENAT_config accepts valid bootstrap_method", {
    # This test checks that the parameter is accepted
    # Full functionality testing is done elsewhere
    config <- TSENAT_config(bootstrap_method = "percentile")
    expect_is(config, "TSENATConfig")
    
    config <- TSENAT_config(bootstrap_method = "bca")
    expect_is(config, "TSENATConfig")
})

test_that("TSENAT_config rejects invalid bootstrap_method", {
    expect_error(
        TSENAT_config(bootstrap_method = "invalid_method"),
        regexp = "should be one of"
    )
})

test_that("TSENAT_config accepts valid sait_method", {
    # Test valid sait_method options
    config <- TSENAT_config(sait_method = "gam")
    expect_is(config, "TSENATConfig")
    
    config <- TSENAT_config(sait_method = "lmm")
    expect_is(config, "TSENATConfig")
})

test_that("TSENAT_config rejects invalid sait_method", {
    expect_error(
        TSENAT_config(sait_method = "invalid_sait"),
        regexp = "should be one of"
    )
})

test_that("TSENAT_config accepts valid sait_pcorr", {
    config <- TSENAT_config(sait_pcorr = "BH")
    expect_is(config, "TSENATConfig")
})

test_that("TSENAT_config rejects invalid sait_pcorr", {
    expect_error(
        TSENAT_config(sait_pcorr = "invalid_pcorr"),
        regexp = "should be one of"
    )
})

test_that("TSENAT_config accepts valid assumptions_checks", {
    config <- TSENAT_config(assumptions_checks = "rank")
    expect_is(config, "TSENATConfig")
})

test_that("TSENAT_config rejects invalid assumptions_checks", {
    expect_error(
        TSENAT_config(assumptions_checks = "invalid_check"),
        regexp = "should be one of"
    )
})

test_that("TSENAT_config match.arg parameters provide consistent error messages", {
    # All match.arg calls should produce similar error messages
    
    error1 <- tryCatch(
        TSENAT_config(bootstrap_method = "wrong"),
        error = function(e) e$message
    )
    expect_match(error1, "should be one of")
    
    error2 <- tryCatch(
        TSENAT_config(sait_method = "wrong"),
        error = function(e) e$message
    )
    expect_match(error2, "should be one of")
})

test_that("TSENAT_config default parameters work correctly", {
    # Test that when parameters are omitted, defaults are used correctly
    config_default <- TSENAT_config()
    config_explicit <- TSENAT_config(bootstrap_method = "percentile")
    
    # Both should create valid configs
    expect_is(config_default, "TSENATConfig")
    expect_is(config_explicit, "TSENATConfig")
})

# ============================================================================
# Parameter Validation Function Tests
# ============================================================================

test_that(".validate_tsenat_config_params rejects NULL sample_col", {
    skip_on_bioc()
    
    expect_error(
        TSENAT:::.validate_tsenat_config_params(
            sample_col = NULL,
            condition_col = "condition",
            q = 1.0,
            divergence_ci = 0.95,
            paired = FALSE,
            subject_col = NULL,
            control = NULL
        ),
        "'sample_col' is required"
    )
})

test_that(".validate_tsenat_config_params rejects non-character sample_col", {
    skip_on_bioc()
    
    expect_error(
        TSENAT:::.validate_tsenat_config_params(
            sample_col = 123,
            condition_col = "condition",
            q = 1.0,
            divergence_ci = 0.95,
            paired = FALSE,
            subject_col = NULL,
            control = NULL
        ),
        "'sample_col' is required"
    )
})

test_that(".validate_tsenat_config_params rejects NULL condition_col", {
    skip_on_bioc()
    
    expect_error(
        TSENAT:::.validate_tsenat_config_params(
            sample_col = "sample",
            condition_col = NULL,
            q = 1.0,
            divergence_ci = 0.95,
            paired = FALSE,
            subject_col = NULL,
            control = NULL
        ),
        "'condition_col' is required"
    )
})

test_that(".validate_tsenat_config_params rejects NULL q", {
    skip_on_bioc()
    
    expect_error(
        TSENAT:::.validate_tsenat_config_params(
            sample_col = "sample",
            condition_col = "condition",
            q = NULL,
            divergence_ci = 0.95,
            paired = FALSE,
            subject_col = NULL,
            control = NULL
        ),
        "'q' is required"
    )
})

test_that(".validate_tsenat_config_params rejects out-of-range q values", {
    skip_on_bioc()
    
    expect_error(
        TSENAT:::.validate_tsenat_config_params(
            sample_col = "sample",
            condition_col = "condition",
            q = c(-0.5, 1.0),
            divergence_ci = 0.95,
            paired = FALSE,
            subject_col = NULL,
            control = NULL
        ),
        "'q' must be numeric value"
    )
})

test_that(".validate_tsenat_config_params rejects q > 2", {
    skip_on_bioc()
    
    expect_error(
        TSENAT:::.validate_tsenat_config_params(
            sample_col = "sample",
            condition_col = "condition",
            q = 2.5,
            divergence_ci = 0.95,
            paired = FALSE,
            subject_col = NULL,
            control = NULL
        ),
        "'q' must be numeric value"
    )
})

test_that(".validate_tsenat_config_params rejects invalid divergence_ci", {
    skip_on_bioc()
    
    expect_error(
        TSENAT:::.validate_tsenat_config_params(
            sample_col = "sample",
            condition_col = "condition",
            q = 1.0,
            divergence_ci = 1.5,
            paired = FALSE,
            subject_col = NULL,
            control = NULL
        ),
        "'divergence_ci' must be a probability"
    )
})

test_that(".validate_tsenat_config_params rejects divergence_ci = 0", {
    skip_on_bioc()
    
    expect_error(
        TSENAT:::.validate_tsenat_config_params(
            sample_col = "sample",
            condition_col = "condition",
            q = 1.0,
            divergence_ci = 0,
            paired = FALSE,
            subject_col = NULL,
            control = NULL
        ),
        "'divergence_ci' must be a probability"
    )
})

test_that(".validate_tsenat_config_params requires subject_col for paired=TRUE", {
    skip_on_bioc()
    
    expect_error(
        TSENAT:::.validate_tsenat_config_params(
            sample_col = "sample",
            condition_col = "condition",
            q = 1.0,
            divergence_ci = 0.95,
            paired = TRUE,
            subject_col = NULL,
            control = "A"
        ),
        "Paired design.*requires"
    )
})

test_that(".validate_tsenat_config_params requires control for paired=TRUE", {
    skip_on_bioc()
    
    expect_error(
        TSENAT:::.validate_tsenat_config_params(
            sample_col = "sample",
            condition_col = "condition",
            q = 1.0,
            divergence_ci = 0.95,
            paired = TRUE,
            subject_col = "subject",
            control = NULL
        ),
        "Paired design.*requires"
    )
})

test_that(".validate_tsenat_config_params accepts valid paired configuration", {
    skip_on_bioc()
    
    expect_invisible(
        TSENAT:::.validate_tsenat_config_params(
            sample_col = "sample",
            condition_col = "condition",
            q = 1.0,
            divergence_ci = 0.95,
            paired = TRUE,
            subject_col = "subject",
            control = "A"
        )
    )
})

test_that("TSENAT_config with extra arguments", {
    skip_on_bioc()
    
    cfg <- TSENAT_config(
        q = 1.0,
        custom_param = "custom_value"
    )
    
    expect_true("custom_param" %in% names(cfg))
    expect_equal(cfg$custom_param, "custom_value")
})

test_that("TSENAT_config warns about metadata in config", {
    skip_on_bioc()
    
    expect_warning(
        TSENAT_config(
            q = 1.0,
            metadata = list(some = "data")
        ),
        "metadata.*should not be in config"
    )
})

test_that("TSENAT_config with multiple q-values", {
    skip_on_bioc()
    
    cfg <- TSENAT_config(
        q = c(0.5, 1.0, 1.5, 2.0)
    )
    
    expect_equal(cfg$q, c(0.5, 1.0, 1.5, 2.0))
})

test_that("TSENAT_config with all numeric q boundary values", {
    skip_on_bioc()
    
    cfg0 <- TSENAT_config(q = 0)
    expect_equal(cfg0$q, 0)
    
    cfg2 <- TSENAT_config(q = 2)
    expect_equal(cfg2$q, 2)
})
