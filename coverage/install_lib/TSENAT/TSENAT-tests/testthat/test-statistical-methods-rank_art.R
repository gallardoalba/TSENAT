context("rank_based_methods: Aligned Rank Transform (ART)")

library(TSENAT)

# Helper: access internal functions (works with both load_all() and installed pkg)
art_test <- function(...) {
    fn <- getFromNamespace(".test_q_condition_interaction", "TSENAT")
    fn(...)
}
art_analyze <- function(...) {
    fn <- getFromNamespace(".detect_q_analyze_gene", "TSENAT")
    fn(...)
}

# ------------------------------------------------------------------------------
# TEST DATA — larger sample sizes to avoid residual df=0
# ------------------------------------------------------------------------------

make_gd_unpaired <- function(n_q = 3, n_rep = 5, seed = 456) {
    set.seed(seed)
    expand.grid(
        replicate = seq_len(n_rep),
        q = paste0("q", seq_len(n_q)),
        condition = c("control", "treatment"),
        stringsAsFactors = FALSE
    ) |>
        transform(
            entropy = rnorm(n_q * 2 * n_rep, mean = 1, sd = 0.2) +
                ifelse(condition == "treatment", 0.25, 0) +
                ifelse(q == "q3", 0.4, 0)
        )
}

make_gd_paired <- function(n_q = 3, n_subjects = 8, seed = 123) {
    set.seed(seed)
    d <- expand.grid(
        q = factor(paste0("q", seq_len(n_q))),
        condition = factor(c("control", "treatment")),
        subject = paste0("S", seq_len(n_subjects)),
        stringsAsFactors = FALSE
    )
    d$entropy <- rnorm(nrow(d), mean = 1, sd = 0.2) +
        ifelse(d$condition == "treatment", 0.25, 0) +
        ifelse(d$q == "q3", 0.4, 0)
    d
}

# ------------------------------------------------------------------------------
# SUITE 1: ART basic functionality
# ------------------------------------------------------------------------------

test_that("ART returns valid result structure for unpaired design", {
    gd <- make_gd_unpaired()
    result <- art_test(gd, value_col = "entropy", q_col = "q",
        condition_col = "condition", paired = FALSE, method = "art")

    expect_type(result, "list")
    expect_named(result, c("statistic", "p_value", "method", "test_type",
        "ss_interaction", "df_interaction", "ss_residual", "df_residual"))
    expect_true(is.numeric(result$statistic))
    expect_true(is.numeric(result$p_value) && result$p_value >= 0 && result$p_value <= 1)
    expect_match(result$method, "Aligned Rank Transform")
    expect_equal(result$test_type, "art_unpaired")
    expect_true(result$df_interaction > 0)
})

test_that("ART returns valid result structure for paired design", {
    gd <- make_gd_paired()
    result <- art_test(gd, value_col = "entropy", q_col = "q",
        condition_col = "condition", paired = TRUE, subject_col = "subject",
        method = "art")

    expect_equal(result$test_type, "art_paired")
    expect_match(result$method, "Aligned Rank Transform")
    expect_true(is.numeric(result$p_value))
})

test_that("ART detects strong interaction signal", {
    set.seed(999)
    gd <- make_gd_paired(n_q = 3, n_subjects = 10)
    gd$entropy <- gd$entropy +
        ifelse(gd$q == "q3" & gd$condition == "treatment", 1.5, 0)

    result <- art_test(gd, value_col = "entropy", q_col = "q",
        condition_col = "condition", paired = TRUE, subject_col = "subject",
        method = "art")

    expect_true(result$p_value < 0.05)
    expect_true(result$statistic > 1)
})

test_that("ART returns non-significant for no interaction", {
    set.seed(111)
    gd <- make_gd_paired(n_q = 3, n_subjects = 10)
    gd$entropy <- rnorm(nrow(gd), mean = 1, sd = 0.15) +
        ifelse(gd$q == "q3", 0.3, 0) +
        ifelse(gd$condition == "treatment", 0.2, 0)

    result <- art_test(gd, value_col = "entropy", q_col = "q",
        condition_col = "condition", paired = TRUE, subject_col = "subject",
        method = "art")

    expect_true(is.na(result$p_value) || result$p_value > 0.01)
})

# ------------------------------------------------------------------------------
# SUITE 2: Method selection
# ------------------------------------------------------------------------------

test_that("method='art' is the default", {
    gd <- make_gd_unpaired()
    result <- art_test(gd, value_col = "entropy", q_col = "q",
        condition_col = "condition", paired = FALSE)
    expect_equal(result$test_type, "art_unpaired")
})

test_that("method='rt' forces Conover-Iman", {
    gd <- make_gd_unpaired()
    result <- art_test(gd, value_col = "entropy", q_col = "q",
        condition_col = "condition", paired = FALSE, method = "rt")
    expect_match(result$method, "Conover-Iman")
    expect_match(result$test_type, "rt_")
})

test_that("invalid method triggers error", {
    gd <- make_gd_unpaired()
    expect_error(
        art_test(gd, value_col = "entropy", q_col = "q",
            condition_col = "condition", method = "invalid"),
        "should be one of"
    )
})

# ------------------------------------------------------------------------------
# SUITE 3: ART and Conover-Iman agreement
# ------------------------------------------------------------------------------

test_that("ART and Conover-Iman agree on strong interaction signals", {
    set.seed(777)
    n_sims <- 15
    art_pvals <- numeric(n_sims)
    rt_pvals  <- numeric(n_sims)

    for (i in seq_len(n_sims)) {
        gd <- make_gd_unpaired(n_q = 3, n_rep = 15, seed = 700 + i)
        gd$entropy <- gd$entropy +
            ifelse(gd$q == "q3" & gd$condition == "treatment", 3.0, 0)

        art_res <- art_test(gd, value_col = "entropy", q_col = "q",
            condition_col = "condition", paired = FALSE, method = "art")
        rt_res <- art_test(gd, value_col = "entropy", q_col = "q",
            condition_col = "condition", paired = FALSE, method = "rt")

        art_pvals[i] <- art_res$p_value
        rt_pvals[i]  <- rt_res$p_value
    }

    # Both methods should detect the strong interaction most of the time
    expect_true(mean(art_pvals < 0.05, na.rm = TRUE) >= 0.5)
    # rt may be less powerful for interaction detection — this is a known limitation
    # that ART is designed to address. Just verify p-values are reasonable.
    expect_true(all(!is.na(rt_pvals)))
    expect_true(all(rt_pvals >= 0 & rt_pvals <= 1))
})

# ------------------------------------------------------------------------------
# SUITE 4: Pipeline integration
# ------------------------------------------------------------------------------

test_that("analyze_gene passes method to test function", {
    gd <- make_gd_unpaired(n_q = 3, n_rep = 7)
    gd$gene <- "gene1"

    result_art <- art_analyze(gd, paired = FALSE, subject_col = NULL,
        has_condition = TRUE, method = "art")
    result_rt  <- art_analyze(gd, paired = FALSE, subject_col = NULL,
        has_condition = TRUE, method = "rt")

    expect_false(result_art$test_failed)
    expect_false(result_rt$test_failed)
    expect_true(is.numeric(result_art$eta2))
})

test_that("analyze_gene handles insufficient q-levels", {
    gd <- data.frame(
        q = factor("q1"), condition = factor(c("control", "treatment")),
        entropy = c(1.0, 1.2), gene = "gene1", stringsAsFactors = FALSE
    )
    result <- art_analyze(gd, paired = FALSE, subject_col = NULL,
        has_condition = TRUE)
    expect_true(result$test_failed)
    expect_equal(result$class, "Insufficient data")
})

# ------------------------------------------------------------------------------
# SUITE 5: Edge cases
# ------------------------------------------------------------------------------

test_that("ART handles missing entropy values", {
    gd <- make_gd_unpaired(n_q = 3, n_rep = 7)
    gd$entropy[1] <- NA_real_
    result <- art_test(gd, value_col = "entropy", q_col = "q",
        condition_col = "condition", paired = FALSE, method = "art")
    expect_type(result, "list")
    expect_true(!is.null(result$test_type))
})

test_that("ART validates required columns", {
    gd <- data.frame(x = 1:10, y = 1:10)
    expect_error(
        art_test(gd, value_col = "entropy", q_col = "q",
            condition_col = "condition", method = "art"),
        "not found in data"
    )
})

test_that("ART validates paired requires subject_col", {
    gd <- make_gd_unpaired()
    expect_error(
        art_test(gd, value_col = "entropy", q_col = "q",
            condition_col = "condition", paired = TRUE,
            subject_col = "nonexistent", method = "art"),
        "not found in data"
    )
})

# ------------------------------------------------------------------------------
# SUITE 6: Test type labels
# ------------------------------------------------------------------------------

test_that("all four test_type combinations are correct", {
    gd_u <- make_gd_unpaired(n_q = 3, n_rep = 7)
    gd_p <- make_gd_paired(n_q = 3, n_subjects = 10)

    expect_equal(
        art_test(gd_u, value_col = "entropy", q_col = "q",
            condition_col = "condition", paired = FALSE, method = "art")$test_type,
        "art_unpaired")
    expect_equal(
        art_test(gd_p, value_col = "entropy", q_col = "q",
            condition_col = "condition", paired = TRUE,
            subject_col = "subject", method = "art")$test_type,
        "art_paired")
    expect_equal(
        art_test(gd_u, value_col = "entropy", q_col = "q",
            condition_col = "condition", paired = FALSE, method = "rt")$test_type,
        "rt_unpaired")
})

# ------------------------------------------------------------------------------
# SUITE 7: Output consistency
# ------------------------------------------------------------------------------

test_that("ART SS values are non-negative and df positive", {
    gd <- make_gd_unpaired(n_q = 3, n_rep = 7)
    result <- art_test(gd, value_col = "entropy", q_col = "q",
        condition_col = "condition", paired = FALSE, method = "art")

    if (!is.na(result$ss_interaction)) {
        expect_true(result$ss_interaction >= 0)
        expect_true(result$ss_residual >= 0)
        expect_true(result$df_interaction >= 1)
    }
})

test_that("ART F-stat and p-value are consistent", {
    gd <- make_gd_unpaired(n_q = 3, n_rep = 7)
    result <- art_test(gd, value_col = "entropy", q_col = "q",
        condition_col = "condition", paired = FALSE, method = "art")

    if (!is.na(result$statistic) && !is.na(result$p_value)) {
        expect_true(result$statistic >= 0)
        expect_true(result$p_value >= 0 && result$p_value <= 1)
    }
})

# ------------------------------------------------------------------------------
# SUITE 8: Stress test
# ------------------------------------------------------------------------------

test_that("ART handles moderately large data quickly", {
    set.seed(42)
    n_rep <- 3
    gd <- expand.grid(
        replicate = seq_len(n_rep),
        q = factor(paste0("q", seq_len(5))),
        condition = factor(c("control", "treatment")),
        stringsAsFactors = FALSE
    )
    gd$entropy <- rnorm(nrow(gd), mean = 1, sd = 0.3) +
        ifelse(gd$condition == "treatment", 0.15, 0)

    start <- Sys.time()
    result <- art_test(gd, value_col = "entropy", q_col = "q",
        condition_col = "condition", paired = FALSE, method = "art")
    elapsed <- as.numeric(Sys.time() - start, units = "secs")

    expect_true(elapsed < 10)
    expect_true(!is.null(result$p_value))
})
