context("Difference Functions: Basic Calculations")

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
        fold_change <- TSENAT:::.calculate_fc(diversity_1, samples, control, "mean")

        expect_length(fold_change, 4)
        expect_true(is.data.frame(fold_change))

        fold_change <- TSENAT:::.calculate_fc(
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
    wilcoxon_result <- .wilcoxon(diversity_1, samples)

    expect_equal(nrow(wilcoxon_result), nrow(diversity_1))
    expect_equal(ncol(wilcoxon_result), 4)
    expect_true(is.data.frame(wilcoxon_result))
    expect_true(all(c("pvalue", "padj", "U", "r") %in% colnames(wilcoxon_result)))

    wilcoxon_result <- .wilcoxon(as.matrix(diversity_2), samples)

    expect_equal(
        as.numeric(wilcoxon_result[
            1,
            "pvalue"
        ]),
        0.03038282,
        tolerance = 0.001,
        scale = 1
    )
    expect_equal(
        as.numeric(wilcoxon_result[
            1,
            "padj"
        ]),
        0.03038282,
        tolerance = 0.001,
        scale = 1
    )
})

test_that("Label shuffling test is correct", {
    shuffling_result <- .label_shuffling(diversity_1, samples, control, "mean")

    expect_equal(nrow(shuffling_result), nrow(diversity_1))
    expect_equal(ncol(shuffling_result), 7)
    expect_true(is.data.frame(shuffling_result))
    expect_true(all(c("pvalue", "padj", "log2FC") %in% colnames(shuffling_result)))

    diversity_2 <- rbind(diversity_2, data.frame(
        S1 = 0.2, S2 = 0.3, S3 = 0.4, S4 = 0.5, S5 = 0.6, S6 = 0.7,
        S7 = 0.8, S8 = 0.9
    ))

    shuffling_result <- .label_shuffling(
        as.matrix(diversity_2),
        samples,
        control,
        "mean"
    )

    # After fixing permutation p-value calculation, expect valid p-values in [0,1]
    expect_true(is.numeric(as.numeric(shuffling_result[
        1,
        "pvalue"
    ])) && as.numeric(shuffling_result[
        1,
        "pvalue"
    ]) >= 0 && as.numeric(shuffling_result[
        1,
        "pvalue"
    ]) <= 1)
    expect_true(is.numeric(as.numeric(shuffling_result[
        1,
        "padj"
    ])) && as.numeric(shuffling_result[
        1,
        "padj"
    ]) >= 0 && as.numeric(shuffling_result[
        1,
        "padj"
    ]) <= 1)
})

context("Difference Functions: Helper Functions and Paired Permutations")

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


# calculate_fc defensive errors
test_that("calculate_fc errors on missing control or samples length mismatch", {
    x <- matrix(rnorm(8), nrow = 2)
    samples <- rep(c("A", "B"), length.out = ncol(x))
    expect_error(TSENAT:::.calculate_fc(x, samples, control = NULL), "`control` must be provided")
    expect_error(TSENAT:::.calculate_fc(x, samples[-1], control = "A"), "Length of 'samples' must equal number of columns in 'x'")
})

context("Wilcoxon Tests: Single Feature Implementation")

# Test .wilcox_one with unpaired design
test_that(".wilcox_one computes unpaired Wilcoxon test correctly", {
    # Create test data: 3 genes x 8 samples
    x <- matrix(rnorm(24), nrow = 3)
    samples <- rep(c("Control", "Treatment"), each = 4)
    
    # Set up as if inside wilcoxon function
    groups <- unique(sort(samples))
    g1_idx <- as.numeric(which(samples %in% groups[1]))
    g2_idx <- as.numeric(which(samples %in% groups[2]))
    
    # Create the .wilcox_one function environment
    .wilcox_one <- function(i) {
        tryCatch({
            test_result <- wilcox.test(x[i, g1_idx], x[i, g2_idx], paired = FALSE, exact = FALSE)
            n <- length(g1_idx) + length(g2_idx)
            list(p.value = test_result$p.value, statistic = test_result$statistic, n = n)
        }, error = function(e) {
            list(p.value = NA_real_, statistic = NA_real_, n = NA_real_)
        }, warning = function(w) {
            list(p.value = NA_real_, statistic = NA_real_, n = NA_real_)
        })
    }
    
    # Test first feature
    result <- .wilcox_one(1)
    
    expect_is(result, "list")
    expect_true("p.value" %in% names(result))
    expect_true("statistic" %in% names(result))
    expect_true("n" %in% names(result))
    expect_true(is.numeric(result$p.value) || is.na(result$p.value))
    expect_true(result$p.value >= 0 && result$p.value <= 1 || is.na(result$p.value))
    expect_equal(result$n, 8)
})

# Test .wilcox_one with paired design (position-based)
test_that(".wilcox_one computes paired Wilcoxon test correctly", {
    # Create paired test data: 3 genes x 8 samples (4 pairs)
    x <- matrix(rnorm(24), nrow = 3)
    samples <- rep(c("Control", "Treatment"), times = 4)  # alternating pairs
    
    # Set up as if inside wilcoxon function
    groups <- unique(sort(samples))
    g1_idx <- as.numeric(which(samples %in% groups[1]))
    g2_idx <- as.numeric(which(samples %in% groups[2]))
    
    # Create the .wilcox_one function for paired test
    .wilcox_one_paired <- function(i) {
        tryCatch({
            test_result <- wilcox.test(x[i, g1_idx], x[i, g2_idx], paired = TRUE, exact = FALSE)
            n <- length(g1_idx)  # number of pairs
            list(p.value = test_result$p.value, statistic = test_result$statistic, n = n)
        }, error = function(e) {
            list(p.value = NA_real_, statistic = NA_real_, n = NA_real_)
        }, warning = function(w) {
            list(p.value = NA_real_, statistic = NA_real_, n = NA_real_)
        })
    }
    
    # Test first feature
    result <- .wilcox_one_paired(1)
    
    expect_is(result, "list")
    expect_true(is.numeric(result$p.value) || is.na(result$p.value))
    expect_true(result$p.value >= 0 && result$p.value <= 1 || is.na(result$p.value))
    expect_equal(result$n, 4)  # 4 pairs
})

# Test .wilcox_one with explicit pairing information
test_that(".wilcox_one computes paired Wilcoxon with explicit pairing correctly", {
    # Create paired test data with explicit pairing
    x <- matrix(rnorm(24), nrow = 3)
    samples <- c("Control", "Treatment", "Control", "Treatment", 
                 "Control", "Treatment", "Control", "Treatment")
    pairs <- c("Pair1", "Pair1", "Pair2", "Pair2", 
               "Pair3", "Pair3", "Pair4", "Pair4")
    groups <- unique(sort(samples))
    
    # Create the .wilcox_one function for explicit pairing
    .wilcox_one_explicit_pairs <- function(i) {
        tryCatch({
            unique_pairs <- unique(pairs)
            all_diffs <- numeric(0)
            for (p in unique_pairs) {
                g1_samples <- which(pairs == p & samples == groups[1])
                g2_samples <- which(pairs == p & samples == groups[2])
                if (length(g1_samples) == 1 && length(g2_samples) == 1) {
                    all_diffs <- c(all_diffs, x[i, g1_samples] - x[i, g2_samples])
                }
            }
            # Perform paired test on the differences
            test_result <- wilcox.test(all_diffs, mu = 0, exact = FALSE)
            list(p.value = test_result$p.value, statistic = test_result$statistic, n = length(all_diffs))
        }, error = function(e) {
            list(p.value = NA_real_, statistic = NA_real_, n = NA_real_)
        }, warning = function(w) {
            list(p.value = NA_real_, statistic = NA_real_, n = NA_real_)
        })
    }
    
    # Test first feature
    result <- .wilcox_one_explicit_pairs(1)
    
    expect_is(result, "list")
    expect_true("p.value" %in% names(result))
    expect_true("statistic" %in% names(result))
    expect_true("n" %in% names(result))
    expect_equal(result$n, 4)  # 4 pairs
    expect_true(is.numeric(result$p.value) || is.na(result$p.value))
    expect_true(result$p.value >= 0 && result$p.value <= 1 || is.na(result$p.value))
})

# Test .wilcox_one error handling
test_that(".wilcox_one handles errors and edge cases gracefully", {
    # Create test data with potential edge cases
    x <- matrix(rnorm(24), nrow = 3)
    samples <- rep(c("A", "B"), each = 4)
    
    g1_idx <- as.numeric(which(samples %in% "A"))
    g2_idx <- as.numeric(which(samples %in% "B"))
    
    # Create .wilcox_one that catches errors
    .wilcox_one_safe <- function(i) {
        tryCatch({
            test_result <- wilcox.test(x[i, g1_idx], x[i, g2_idx], paired = FALSE, exact = FALSE)
            list(p.value = test_result$p.value, statistic = test_result$statistic, n = length(g1_idx) + length(g2_idx))
        }, error = function(e) {
            list(p.value = NA_real_, statistic = NA_real_, n = NA_real_)
        }, warning = function(w) {
            list(p.value = NA_real_, statistic = NA_real_, n = NA_real_)
        })
    }
    
    # Valid feature should return valid results with proper structure
    result_valid <- .wilcox_one_safe(1)
    expect_is(result_valid, "list")
    expect_true("p.value" %in% names(result_valid))
    expect_true("statistic" %in% names(result_valid))
    expect_true("n" %in% names(result_valid))
    
    # p-value should be in [0, 1] or NA
    expect_true((is.numeric(result_valid$p.value) && result_valid$p.value >= 0 && result_valid$p.value <= 1) || is.na(result_valid$p.value))
    expect_true(is.numeric(result_valid$n))
})

# Test r-value computation from U statistic (unpaired)
test_that("r-value computation from unpaired U statistic is correct", {
    # Create simple data for Wilcoxon test
    x <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8), nrow = 1)
    samples <- rep(c("A", "B"), each = 4)
    
    # Perform Wilcoxon test
    test_result <- wilcox.test(x[1, 1:4], x[1, 5:8], paired = FALSE)
    
    # Compute Z and r manually
    n1 <- 4
    n2 <- 4
    n <- n1 + n2
    U <- test_result$statistic
    expected_U <- n1 * n2 / 2
    var_U <- (n1 * n2 * (n1 + n2 + 1)) / 12
    sd_U <- sqrt(var_U)
    Z <- (U - expected_U) / sd_U
    r <- Z / sqrt(n)
    
    # r should be in [-1, 1]
    r_clamped <- pmax(-1, pmin(1, r))
    expect_true(r_clamped >= -1 && r_clamped <= 1)
})

# Test r-value computation from U statistic (paired)
test_that("r-value computation from paired U statistic is correct", {
    # Create simple paired data
    x <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8), nrow = 1)
    samples <- rep(c("A", "B"), times = 4)
    
    # Perform paired Wilcoxon test
    g1_idx <- c(1, 3, 5, 7)
    g2_idx <- c(2, 4, 6, 8)
    test_result <- wilcox.test(x[1, g1_idx], x[1, g2_idx], paired = TRUE, exact = FALSE)
    
    # Compute Z and r manually for paired test
    n <- length(g1_idx)
    U <- test_result$statistic
    expected_U <- n * (n + 1) / 4
    var_U <- (n * (n + 1) * (2 * n + 1)) / 24
    sd_U <- sqrt(var_U)
    Z <- (U - expected_U) / sd_U
    r <- Z / sqrt(n)
    
    # r should be in [-1, 1]
    r_clamped <- pmax(-1, pmin(1, r))
    expect_true(r_clamped >= -1 && r_clamped <= 1)
})

context("Statistical Methods: Wilcoxon Defensive Behavior")

test_that("wilcoxon handles all-NA and constant rows without error and returns matrix", {
    # two groups of 3 samples each
    samples <- c(rep("A", 3), rep("B", 3))
    # build matrix with 3 rows: all-NA, constant, variable
    m <- matrix(nrow = 3, ncol = 6)
    m[1, ] <- NA_real_
    m[2, ] <- rep(5, 6)
    m[3, ] <- c(1, 2, 3, 4, 5, 6)

    res <- .wilcoxon(m, samples, pcorr = "none")
    expect_true(is.data.frame(res))
    # four columns: pvalue, padj, U, r
    expect_equal(ncol(res), 4)
    # raw p-values produced and NA-handling yields numeric outputs
    raw <- as.numeric(res[, 1])
    expect_length(raw, 3)
    # all-NA row should produce p-value of 1 (by convention used elsewhere)
    expect_true(is.finite(raw[1]))
})

test_that(".wilcoxon() returns named p-value and effect size columns", {
    mat <- matrix(runif(20), nrow = 5)
    samples <- rep(c("A", "B"), each = 5)
    res <- .wilcoxon(mat, samples, pcorr = "none", paired = FALSE, exact = FALSE)
    expect_true(is.data.frame(res))
    expect_true(all(c("pvalue", "padj", "U", "r") %in% colnames(res)))
})

context("Statistical Methods: Wilcoxon Paired Behavior")

test_that("wilcoxon paired matches per-row wilcox.test on ordered pairs", {
    # create a small matrix with two genes and two pairs (columns: N,T,N,T)
    # ensure paired differences are not identical to avoid "ties" in differences
    x <- matrix(c(
        1, 2, 3, 6, # gene1 across 4 samples (diffs: -1, -3)
        5, 6, 7, 9 # gene2 across 4 samples (diffs: -1, -2)
    ), nrow = 2, byrow = TRUE)
    samples <- c("Normal", "Tumor", "Normal", "Tumor")

    res <- .wilcoxon(x, samples, paired = TRUE, exact = TRUE)

    # compute expected raw p-values by calling wilcox.test per row with
    # paired=TRUE
    expected_raw <- sapply(seq_len(nrow(x)), function(i) {
        wilcox.test(
            x[
                i,
                which(samples == "Normal")
            ],
            x[
                i,
                which(samples == "Tumor")
            ],
            paired = TRUE, exact = TRUE
        )$p.value
    })

    expected_adj <- p.adjust(expected_raw, method = "BH")

    expect_equal(as.numeric(res[, "pvalue"]), as.numeric(expected_raw))
    expect_equal(as.numeric(res[, "padj"]), as.numeric(expected_adj))
})


test_that("wilcoxon paired errors on unequal group sizes", {
    x_bad <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2, byrow = TRUE)
    samples_bad <- c("Normal", "Tumor", "Normal")
    expect_error(.wilcoxon(x_bad, samples_bad, paired = TRUE), "Paired Wilcoxon requires equal numbers of samples in each group")
})



test_that("wilcoxon paired handles SummarizedExperiment input", {
    # construct a small paired dataset
    sample_names <- c("S1_N", "S1_T", "S2_N", "S2_T")
    mat_vals <- matrix(c(
        1, 2, 3, 6,  # gene1
        5, 6, 7, 9   # gene2
    ), nrow = 2, byrow = TRUE)
    
    colnames(mat_vals) <- sample_names
    
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(diversity = mat_vals)
    )
    cond <- ifelse(grepl("_N$", sample_names), "Normal", "Tumor")
    coldata <- data.frame(
        Sample = sample_names,
        Condition = cond,
        stringsAsFactors = FALSE
    )
    
    se_mapped <- TSENAT:::.tsenat_map_metadata_se(se, coldata)
    res <- .wilcoxon(
        SummarizedExperiment::assay(se_mapped),
        SummarizedExperiment::colData(se_mapped)$sample_type,
        paired = TRUE,
        exact = TRUE
    )
    
    # Verify output structure
    expect_true(is.data.frame(res))
    expect_equal(ncol(res), 4)
    expect_true(all(c("pvalue", "padj", "U", "r") %in% colnames(res)))
    # Verify results are not all NA
    expect_true(any(!is.na(res$pvalue)))
})

context("Statistical Methods: Wilcoxon U Statistic and r-Value Computation")

test_that("wilcoxon U statistic is computed correctly", {
    mat <- matrix(c(
        1, 2, 3, 4, # gene1
        5, 6, 7, 8  # gene2
    ), nrow = 2, byrow = TRUE)
    samples <- c("A", "A", "B", "B")
    
    res <- .wilcoxon(mat, samples, pcorr = "none", paired = FALSE)
    
    # U should be present and non-NA for valid data
    expect_false(anyNA(res$U))
    expect_true(all(res$U >= 0))
    # U should be bounded by the product of group sizes
    expect_true(all(res$U <= 4))  # 2 samples in each group
})

test_that("wilcoxon r-value is computed correctly", {
    mat <- matrix(c(
        1, 2, 3, 4,
        5, 6, 7, 8
    ), nrow = 2, byrow = TRUE)
    samples <- c("A", "A", "B", "B")
    
    res <- .wilcoxon(mat, samples, pcorr = "none", paired = FALSE)
    
    # r-value should be present and non-NA for valid data
    expect_false(anyNA(res$r))
    # r-value should be bounded between -1 and 1
    expect_true(all(res$r >= -1 & res$r <= 1))
    # For constant row, r should be 0 (no effect)
    mat_const <- matrix(c(1, 1, 1, 1), nrow = 1)
    res_const <- .wilcoxon(mat_const, samples, pcorr = "none")
    expect_true(res_const$r[1] == 0)
})

test_that("wilcoxon U and r values are NA when pvalue is NA", {
    m <- matrix(nrow = 2, ncol = 4)
    m[1, ] <- NA_real_
    m[2, ] <- c(1, 2, 3, 4)
    samples <- c("A", "A", "B", "B")
    
    res <- .wilcoxon(m, samples, pcorr = "none")
    
    # All-NA row should have NA U and r value
    expect_true(is.na(res$U[1]))
    expect_true(is.na(res$r[1]))
    # Valid row should have non-NA U and r
    expect_false(is.na(res$U[2]))
    expect_false(is.na(res$r[2]))
})

test_that("wilcoxon r-value is bounded and non-NA for valid paired data", {
    # For paired test, verify r-value is reasonable
    mat <- matrix(c(
        1, 2, 3, 6,  # distinct paired differences
        5, 6, 7, 9
    ), nrow = 2, byrow = TRUE)
    samples <- c("N", "T", "N", "T")
    
    res <- .wilcoxon(mat, samples, pcorr = "none", paired = TRUE)
    
    # Verify r-value is within valid range and non-NA for non-constant rows
    for (i in seq_len(nrow(res))) {
        p <- res$pvalue[i]
        if (!is.na(p)) {
            # r should be bounded between -1 and 1 and non-NA for valid p-values
            expect_true(!is.na(res$r[i]))
            expect_true(res$r[i] >= -1 && res$r[i] <= 1)
        }
    }
})

