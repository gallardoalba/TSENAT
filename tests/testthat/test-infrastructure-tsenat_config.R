context("Orchestration: Configuration Builder")

test_that("tsenat_config creates config with defaults", {
  cfg <- tsenat_config()

  expect_s3_class(cfg, "TSENATConfig")
  expect_true(is.list(cfg))
  expect_true("q_values" %in% names(cfg))
  expect_true("methods" %in% names(cfg))
  expect_equal(cfg$fdr_threshold, 0.05)
  expect_equal(cfg$p_threshold, 0.05)
})

test_that("tsenat_config accepts custom q_values", {
  cfg <- tsenat_config(q_values = c(0.5, 1.0, 1.5))

  expect_equal(cfg$q_values, c(0.5, 1.0, 1.5))
})

test_that("tsenat_config generates q_values from q_range", {
  cfg <- tsenat_config(q_range = c(1.0, 2.0))

  expect_true(1.0 %in% cfg$q_values)
  expect_true(2.0 %in% cfg$q_values)
  expect_true(all(cfg$q_values >= 1.0 & cfg$q_values <= 2.0))
})

test_that("tsenat_config accepts custom thresholds", {
  cfg <- tsenat_config(p_threshold = 0.01, fdr_threshold = 0.001)

  expect_equal(cfg$p_threshold, 0.01)
  expect_equal(cfg$fdr_threshold, 0.001)
})

test_that("tsenat_config accepts formula", {
  form <- ~ treatment + batch
  cfg <- tsenat_config(formula = form)

  expect_equal(cfg$formula, form)
})

test_that("tsenat_config accepts method selection", {
  cfg <- tsenat_config(methods = c("diversity", "jackknife"))

  expect_equal(cfg$methods, c("diversity", "jackknife"))
})

test_that("tsenat_config rejects invalid methods", {
  expect_error(
    tsenat_config(methods = c("diversity", "invalid_method")),
    "Invalid methods"
  )
})

test_that("tsenat_config accepts plot configuration", {
  cfg <- tsenat_config(
    generate_plots = TRUE,
    plot_types = c("q_curve", "divergence")
  )

  expect_equal(cfg$generate_plots, TRUE)
  expect_equal(cfg$plot_types, c("q_curve", "divergence"))
})

test_that("tsenat_config accepts additional parameters", {
  cfg <- tsenat_config(seed = 42, custom_param = "value")

  expect_equal(cfg$seed, 42)
  expect_equal(cfg$custom_param, "value")
})

test_that("tsenat_config has TSENATConfig class", {
  cfg <- tsenat_config()

  expect_true("TSENATConfig" %in% class(cfg))
  expect_true(is.list(cfg))
})
