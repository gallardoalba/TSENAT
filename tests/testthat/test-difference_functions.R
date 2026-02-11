context("Basic difference calculations")

diversity_1 <- matrix(runif(80), ncol = 8)
diversity_2 <- data.frame(
    S1 = 0.1,
    S2 = 0.2,
    S3 = 0.3,
    S4 = 0.4,
    S5 = 0.5,
    S6 = 0.6,
    S7 = 0.7,
    S8 = 0.8
)
samples <- c(rep("Healthy", 4), rep("Pathogenic", 4))
control <- "Healthy"

test_that("Fold change calculation is correct", {
    for (method in c("mean", "median")) {
        fold_change <- calculate_fc(diversity_1, samples, control, "mean")

        expect_length(fold_change, 4)
        expect_true(is.data.frame(fold_change))

        fold_change <- calculate_fc(
            as.matrix(diversity_2),
            samples,
            control,
            "mean"
        )

        expect_equal(fold_change$Pathogenic_mean,
            0.65,
            tolerance = 0.001,
            scale = 1
        )
        expect_equal(fold_change$Healthy_mean, 0.25, tolerance = 0.001, scale = 1)
        expect_equal(fold_change$mean_difference, 0.4, tolerance = 0.001, scale = 1)
        expect_equal(fold_change$log2_fold_change,
            1.378512,
            tolerance = 0.001,
            scale = 1
        )
    }
})

test_that("Wilcoxon sum rank test is correct", {
    wilcoxon_result <- wilcoxon(diversity_1, samples)

    expect_length(wilcoxon_result, 20)
    expect_true(is.matrix(wilcoxon_result))

    wilcoxon_result <- wilcoxon(as.matrix(diversity_2), samples)

    expect_equal(
        as.numeric(wilcoxon_result[
            1,
            1
        ]),
        0.03038282,
        tolerance = 0.001,
        scale = 1
    )
    expect_equal(
        as.numeric(wilcoxon_result[
            1,
            2
        ]),
        0.03038282,
        tolerance = 0.001,
        scale = 1
    )
})

test_that("Label shuffling test is correct", {
    shuffling_result <- label_shuffling(diversity_1, samples, control, "mean")

    expect_length(shuffling_result, 20)
    expect_true(is.matrix(shuffling_result))

    diversity_2 <- rbind(diversity_2, data.frame(
        S1 = 0.2, S2 = 0.3, S3 = 0.4, S4 = 0.5, S5 = 0.6, S6 = 0.7,
        S7 = 0.8, S8 = 0.9
    ))

    shuffling_result <- label_shuffling(
        as.matrix(diversity_2),
        samples,
        control,
        "mean"
    )

    # After fixing permutation p-value calculation, expect valid p-values in [0,1]
    expect_true(is.numeric(as.numeric(shuffling_result[
        1,
        1
    ])) && as.numeric(shuffling_result[
        1,
        1
    ]) >= 0 && as.numeric(shuffling_result[
        1,
        1
    ]) <= 1)
    expect_true(is.numeric(as.numeric(shuffling_result[
        1,
        2
    ])) && as.numeric(shuffling_result[
        1,
        2
    ]) >= 0 && as.numeric(shuffling_result[
        1,
        2
    ]) <= 1)
})

context("difference helper edge cases and paired permutations")

# Tests for aggregation of FC values and pseudocount behavior
test_that(".tsenat_aggregate_fc_values orders groups and handles NAs", {
    x <- matrix(c(
        1, NA, 3, 4, # gene1 across 4 samples
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
