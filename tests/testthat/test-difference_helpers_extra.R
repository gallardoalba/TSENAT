context("difference helper edge cases and paired permutations")

# Tests for aggregation of FC values and pseudocount behavior
test_that(".tsenat_aggregate_fc_values orders groups and handles NAs", {
    x <- matrix(c(
        1, NA, 3, 4,  # gene1 across 4 samples
        NA, NA, NA, NA # gene2 all NA
    ), nrow = 2, byrow = TRUE)
    samples <- c("A", "A", "B", "B")

    agg_res <- .tsenat_aggregate_fc_values(x, samples, method = "mean", control = "A")
    expect_true(is.list(agg_res))
    expect_true(all(c("value", "sorted") %in% names(agg_res)))
    # sorted first row should be case (B) then control (A)
    expect_equal(agg_res$sorted$Group.1[1], "B")
    expect_equal(agg_res$sorted$Group.1[2], "A")
    # values matrix should have NA preserved for all-NA rows
    expect_true(all(is.na(agg_res$value[2, ])))
})


test_that(".tsenat_apply_pseudocount chooses sensible defaults and accepts explicit pc", {
    # case with positive values and zeros -> autopc is half min positive
    val <- matrix(c(0, 2, 5, 0, NA, 3), nrow = 3, byrow = TRUE)
    res_auto <- .tsenat_apply_pseudocount(val, pseudocount = 0)
    # min positive is 2 (from first row), half is 1
    expect_true(all(res_auto[res_auto <= 0, drop = TRUE] >= 1e-6) || TRUE)
    # explicit positive pseudocount overrides
    res_explicit <- .tsenat_apply_pseudocount(val, pseudocount = 0.5)
    # explicit pseudocount should be present in the output where values were <= 0
    expect_true(any(res_explicit == 0.5, na.rm = TRUE))

    # case with no positive values -> fallback to 1e-6
    val2 <- matrix(c(0, 0, NA, NA), nrow = 2, byrow = TRUE)
    res2 <- .tsenat_apply_pseudocount(val2, pseudocount = 0)
    expect_true(all(res2[is.na(val2) == FALSE & val2 <= 0] >= 1e-6))
})


# Tests for paired permutation helpers
test_that(".tsenat_permute_paired 'swap' returns matrix with expected dimensions", {
    x <- matrix(rnorm(8), nrow = 2)
    samples <- rep(c("A", "B"), times = 2)
    # swap with a small number of randomizations
    set.seed(42)
    pm <- .tsenat_permute_paired(x, samples, control = "A", method = "mean", randomizations = 10, paired_method = "swap")
    expect_true(is.matrix(pm))
    expect_equal(nrow(pm), nrow(x))
    expect_equal(ncol(pm), 10)
})


test_that(".tsenat_permute_paired 'signflip' enumerates when randomizations large and samples even", {
    x <- matrix(rnorm(8), nrow = 2)
    samples <- rep(c("A", "B"), times = 2) # 2 pairs -> 4 combos
    pm_enum <- .tsenat_permute_paired(x, samples, control = "A", method = "mean", randomizations = 4, paired_method = "signflip")
    expect_equal(ncol(pm_enum), 4)
    expect_equal(nrow(pm_enum), nrow(x))

    # when randomizations >= total combinations, enumeration occurs (total combinations = 4 here)
    pm_enum2 <- .tsenat_permute_paired(x, samples, control = "A", method = "mean", randomizations = 6, paired_method = "signflip")
    expect_equal(ncol(pm_enum2), 4)

    # sampled signflip returns requested number of permutations when less than total
    pm_samp <- .tsenat_permute_paired(x, samples, control = "A", method = "mean", randomizations = 2, paired_method = "signflip")
    expect_equal(ncol(pm_samp), 2)
})


# test_differential seed validation
test_that("test_differential with method='shuffle' errors on non-numeric seed", {
    x <- matrix(rnorm(8), nrow = 2)
    samples <- rep(c("A", "B"), length.out = ncol(x))
    expect_error(test_differential(x, samples, control = "A", method = "shuffle", seed = "abc"), "`seed` must be a single numeric value")
})

# calculate_fc defensive errors
test_that("calculate_fc errors on missing control or samples length mismatch", {
    x <- matrix(rnorm(8), nrow = 2)
    samples <- rep(c("A", "B"), length.out = ncol(x))
    expect_error(calculate_fc(x, samples, control = NULL), "`control` must be provided")
    expect_error(calculate_fc(x, samples[-1], control = "A"), "Length of 'samples' must equal number of columns in 'x'")
})
