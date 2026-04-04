# ================================================================================
# PAIRED RANK TEST CONDITIONAL SELECTION (March 2026)
# ================================================================================
# Conditional logic to choose between paired rank tests based on data
# characteristics: - Standard Friedman test (default) - Aligned Rank Transform
# Friedman (heteroscedastic): Better power with variance drift - Robust
# Friedman median test (extreme skewness): Handles heavy-tailed distributions
# Mathematical basis: Same characteristic detection as unpaired (Breusch-Pagan,
# skewness) but adapted for paired design with blocking structure preservation.
# ================================================================================

#' Select appropriate paired rank test based on data characteristics
#'
#' Implements conditional logic for paired designs to choose between:
#' - Standard Friedman test (default)
#' - Aligned Rank Transform Friedman (heteroscedastic data)
#' - Mood's robust median test adapted for blocking (extreme skewness)
#'
#' @param data Data frame with entropy values and blocking structure
#' @param value_col Column name for values to test (default: 'entropy')
#' @param group_col Column name for treatment/q-parameter (default: 'q')
#' @param subject_col Column name for blocking variable (default:
#' 'paired_samples')
#' @param verbose Logical: print diagnostic information
#'
#' @return List with:
#' - test_selected: Name of selected test ('friedman', 'art_friedman',
#' 'robust_friedman')
#' - characteristics: List of detected characteristics (heteroscedastic,
#' highly_skewed)
#' - reasons: Character vector of reasons for selection
#'
#' @details
#' **Heteroscedasticity detection (Breusch-Pagan adapted for blocks):**
#' 1. Fit additive model: value ~ subject + treatment
#' 2. Extract residuals
#' 3. Test: residuals2 ~ fitted values (parametric heteroscedasticity)
#' 4. Detect if variance differs by treatment (typical in entropy data)
#'
#' **Extreme skewness detection:**
#' - Compute skewness of values pooled across all treatments
#' - Threshold: |skewness| > 2 indicates extreme heavy tails
#' - Robust median test appropriate when normality assumption strongly violated
#'

#' @noRd

# Helper: Safe extraction with error handling
.safe_extract <- function(expr, default = NA) {
    result <- try(expr, silent = TRUE)
    if (inherits(result, "try-error")) default else result
}

#' @noRd
.select_rank_test_paired <- function(data, value_col = "entropy", group_col = "q",
    subject_col = "paired_samples", verbose = FALSE) {

    values <- data[[value_col]]
    groups <- data[[group_col]]
    subjects <- data[[subject_col]]

    # Ensure factors only if not already (avoid unnecessary conversion overhead)
    groups <- if (!is.factor(groups)) factor(groups) else groups
    subjects <- if (!is.factor(subjects)) factor(subjects) else subjects

    characteristics <- list(heteroscedastic = FALSE, boundary_clustered = FALSE,
        highly_skewed = FALSE, n_groups = nlevels(groups), n_subjects = nlevels(subjects),
        n_values = length(values))

    reasons <- character(0)

    # -----------------------------------------------------------------------------
    # 1. HETEROSCEDASTICITY DETECTION (Breusch-Pagan adapted for blocking)
    # -----------------------------------------------------------------------------

    # Fit additive model with blocking
    additive_mod <- try(lm(values ~ subjects + groups), silent = TRUE)

    if (!inherits(additive_mod, "try-error")) {
        residuals <- residuals(additive_mod)
        fitted <- fitted(additive_mod)

        # Breusch-Pagan test: test if residual variance depends on fitted
        # values
        bp_data <- data.frame(residuals_sq = residuals^2, fitted = fitted)

        bp_mod <- try(lm(residuals_sq ~ fitted, data = bp_data), silent = TRUE)

        if (!inherits(bp_mod, "try-error")) {
            bp_anova <- try(anova(bp_mod), silent = TRUE)

            bp_pvalue <- NA_real_
            var_ratio <- NA_real_

            if (!inherits(bp_anova, "try-error") && nrow(bp_anova) >= 1) {
                bp_pvalue <- .safe_extract(as.numeric(bp_anova[1, "Pr(>F)"]))

                if (!is.na(bp_pvalue)) {
                  # Compute variance ratio across treatment groups (not
                  # subjects)
                  group_vars <- tapply(values, groups, var, na.rm = TRUE)
                  finite_vars <- group_vars[is.finite(group_vars)]
                  var_ratio <- if (length(finite_vars) > 1) max(finite_vars)/min(finite_vars) else 1

                  # Heteroscedasticity detected if p < 0.01 AND variance_ratio > 3
                  # (stricter thresholds to reduce false positives on balanced data)
                  if (!is.na(bp_pvalue) && !is.na(var_ratio) && bp_pvalue < 0.01 &&
                    var_ratio > 3) {
                    characteristics$heteroscedastic <- TRUE
                    reasons <- c(reasons, sprintf("Heteroscedasticity: BP_p=%.4f, var_ratio=%.2f (treatment-based)",
                      bp_pvalue, var_ratio))
                  }
                }
            }
        }
    }

    # -----------------------------------------------------------------------------
    # 2. EXTREME SKEWNESS DETECTION
    # -----------------------------------------------------------------------------

    skewness_val <- .compute_skewness(values)

    if (abs(skewness_val) > 2) {
        characteristics$highly_skewed <- TRUE
        reasons <- c(reasons, sprintf("Extreme skewness: |skew|=%.3f", skewness_val))
    }

    # -----------------------------------------------------------------------------
    # TEST SELECTION LOGIC FOR PAIRED DESIGNS
    # -----------------------------------------------------------------------------

    # Priority: Handle extreme skewness first (robust test needed)
    if (characteristics$highly_skewed) {
        test_selected <- "robust_friedman"
        reasons <- c(reasons, "-> Using robust (median-based) Friedman test")
    } else if (characteristics$heteroscedastic) {
        test_selected <- "art_friedman"
        reasons <- c(reasons, "-> Using Aligned Rank Transform + parametric test (paired)")
    } else {
        test_selected <- "friedman"
        reasons <- c(reasons, "-> Using standard Friedman test (no special characteristics)")
    }

    if (verbose && length(reasons) > 0) {
        message("Paired test selection for ", value_col, ":")
        for (r in reasons) message("  ", r)
    }

    return(list(test_selected = test_selected, characteristics = characteristics,
        reasons = reasons))
}

#' Apply Robust Friedman Test (Mood's Median Test adapted for blocks)
#'
#' Applies a robust version of Friedman test using median-based statistics
#' instead of rank sums. More robust to extreme outliers and heavy-tailed
#' distributions than standard Friedman test.
#'
#' **Mathematical basis:**
#' Standard Friedman ranks observations within blocks and tests if rank sums
#' differ across treatments. This can be affected by extreme outliers.
#'
#' Robust version:
#' 1. Compute block-wise medians for each treatment
#' 2. Compute grand median across all observations
#' 3. Test: Are median differences significantly non-zero?
#' 4. Use exact median inference (sign test logic extended to multiple
#' treatments)
#'
#' @param data Data frame with entropy and blocking structure
#' @param value_col Column name for entropy values
#' @param group_col Column name for treatment (q-values)
#' @param subject_col Column name for blocking variable (subjects/pairs)
#'
#' @return List with components:
#' - statistic: Test statistic (chi-squared form, comparable to Friedman)
#' - p_value: P-value from chi-squared distribution
#' - method: 'Robust (Median-based) Friedman Test'
#' - test_type: 'robust_friedman'
#'

#' @noRd
.apply_robust_friedman <- function(data, value_col = "entropy", group_col = "q",
    subject_col = "paired_samples") {

    # Extract values, groups, subjects
    values <- data[[value_col]]
    groups <- as.factor(data[[group_col]])
    subjects <- as.factor(data[[subject_col]])

    n_treatments <- nlevels(groups)

    # OPTIMIZATION: Vectorized median calculation using by() instead of for-loop
    # Compute all block-group medians in one operation (faster than sequential loops)
    block_list <- by(data.frame(values, groups), subjects, function(d) {
        tapply(d$values, d$groups, median, na.rm = TRUE)
    }, simplify = TRUE)
    
    block_ranks <- do.call(rbind, block_list)
    colnames(block_ranks) <- levels(groups)

    # Count median-based differences (robust approach) Use median test: For
    # each block, is a treatment's value above or below grand median?
    grand_median <- median(values, na.rm = TRUE)

    # OPTIMIZATION: Vectorized binary matrix creation Instead of nested loops,
    # use vectorized comparison
    above_median_matrix <- block_ranks > grand_median
    mode(above_median_matrix) <- "numeric"  # Convert TRUE/FALSE to 1/0

    # Compute chi-squared test for independence H0: Probability of being above
    # median is same for all treatments
    # Create contingency table: rows = above/below median, columns = treatments
    # OPTIMIZATION: Directly flatten matrix and create factors instead of intermediate vectors
    contingency_table <- table(
        Above = factor(c(above_median_matrix), levels = c(0, 1), labels = c("Below", "Above")),
        Treatment = rep(colnames(block_ranks), each = nrow(block_ranks))
    )

    # For robustness: Use exact or simulated p-value (Fisher's exact not
    # practical for large tables) Fallback: Chi-squared test
    chi_test <- tryCatch(chisq.test(contingency_table), error = function(e) NULL)

    if (is.null(chi_test)) {
        # If chi-squared fails, use simpler Friedman-like approach Compute sum
        # of squared deviations of treatment medians from overall median
        treat_medians <- .colMedians(block_ranks)
        Q_stat <- sum((treat_medians - median(treat_medians, na.rm = TRUE))^2, na.rm = TRUE)
        p_val <- 1 - pchisq(Q_stat, df = n_treatments - 1)
    } else {
        Q_stat <- chi_test$statistic
        p_val <- chi_test$p.value
    }

    list(statistic = as.numeric(Q_stat), p_value = as.numeric(p_val), method = "Robust (Median-based) Friedman Test",
        test_type = "robust_friedman")
}

# Helper: Compute column medians (simple implementation)

.colMedians <- function(x) {
    apply(x, 2, median, na.rm = TRUE)
}

#' Apply Aligned Rank Transform Friedman Test
#'
#' Applies Aligned Rank Transform (ART) followed by parametric test for paired
#' designs with heteroscedasticity. More powerful than standard Friedman when
#' variances differ across treatment levels.
#'
#' **Mathematical basis (ART, Wobbrock et al., 2011):**
#' 1. Align: Subtract block effect (within-block median) from each observation
#' 2. Rank: Rank aligned values globally
#' 3. Apply parametric test (e.g., repeated measures ANOVA) on ranks
#' Result: Like Friedman but handles heteroscedasticity better
#'
#' @param data Data frame with entropy and blocking structure
#' @param value_col Column name for entropy values
#' @param group_col Column name for treatment (q-values)
#' @param subject_col Column name for blocking variable
#'
#' @return List with components:
#' - statistic: F-statistic (as alternative to Friedman's Q-statistic)
#' - p_value: P-value from F distribution
#' - method: 'Aligned Rank Transform (ART) Friedman Test'
#' - test_type: 'art_friedman'
#'

#' @noRd
.apply_art_friedman <- function(data, value_col = "entropy", group_col = "q", subject_col = "paired_samples") {

    values <- data[[value_col]]
    groups <- as.factor(data[[group_col]])
    subjects <- as.factor(data[[subject_col]])

    # Step 1: Align (remove block/subject effects) OPTIMIZATION: Vectorized
    # median calculation by subject Use ave() to compute subject medians in one
    # operation instead of loop
    subject_medians <- ave(values, subjects, FUN = function(x) median(x, na.rm = TRUE))
    aligned_values <- values - subject_medians

    # Step 2: Rank aligned values globally
    ranked_values <- rank(aligned_values, na.last = "keep")

    # Step 3: Apply parametric test on ranks Use two-way ANOVA (subject + group
    # interaction) on ranks
    rank_data <- data.frame(ranks = ranked_values, subject = subjects, group = groups)

    anova_mod <- try(aov(ranks ~ subject + group, data = rank_data), silent = TRUE)

    if (inherits(anova_mod, "try-error")) {
        return(list(statistic = NA_real_, p_value = NA_real_, method = "ART Friedman (failed)",
            test_type = "art_friedman"))
    }

    # Extract F-statistic and p-value for group effect
    anova_summary <- summary(anova_mod)

    if (length(anova_summary) > 0 && !is.null(anova_summary[[1]])) {
        F_stat <- anova_summary[[1]]["group", "F value"]
        p_val <- anova_summary[[1]]["group", "Pr(>F)"]
    } else {
        F_stat <- NA_real_
        p_val <- NA_real_
    }

    list(statistic = as.numeric(F_stat), p_value = as.numeric(p_val), method = "Aligned Rank Transform (ART) Friedman Test",
        test_type = "art_friedman")
}
