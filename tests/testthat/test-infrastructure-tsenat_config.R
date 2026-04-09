context("Orchestration: Configuration Builder")

test_that("TSENAT_config creates config with defaults", {
  cfg <- TSENAT_config()

  expect_s3_class(cfg, "TSENATConfig")
  expect_true(is.list(cfg))
  expect_true("q_values" %in% names(cfg))
  expect_equal(cfg$fdr_threshold, 0.05)
  expect_equal(cfg$p_threshold, 0.05)
})

test_that("TSENAT_config accepts custom q_values", {
  cfg <- TSENAT_config(q_values = c(0.5, 1.0, 1.5))

  expect_equal(cfg$q_values, c(0.5, 1.0, 1.5))
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
