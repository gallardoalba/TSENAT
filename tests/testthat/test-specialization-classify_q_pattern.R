library(testthat)

context("classify_q_pattern helper")

# Basic scenarios
per_q1 <- c(q_0.5 = 0.5, q_1 = 0.3, q_2 = 0.1)
per_q2 <- c(q_0.5 = 0.1, q_1 = 0.3, q_2 = 0.5)
per_q3 <- c(q_0.5 = 0.2, q_1 = 0.2, q_2 = 0.2)  # flat values -> balanced

# edge cases
per_q_na <- c(q_0.5 = NA_real_, q_1 = NA_real_)
per_q_short <- c(q_0.5 = 0.2)
per_q_noname <- c(0.1, 0.2, 0.3)  # unnamed vector should return NA

test_that("patterns are classified correctly", {
  expect_equal(classify_q_pattern(per_q1), "RARE_DRIVEN")
  expect_equal(classify_q_pattern(per_q2), "ABUNDANT_DRIVEN")
  expect_equal(classify_q_pattern(per_q3), "BALANCED")
})

test_that("NA or invalid input returns NA", {
  expect_true(is.na(classify_q_pattern(per_q_na)))
  expect_true(is.na(classify_q_pattern(per_q_short)))
  expect_true(is.na(classify_q_pattern(per_q_noname)))
  expect_true(is.na(classify_q_pattern(NULL)))
})
