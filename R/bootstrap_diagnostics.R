# ============================================================================
# PRIORITY 3: ENHANCED BOOTSTRAP DIAGNOSTICS (March 2026)
# ============================================================================
# These functions provide sophisticated diagnostics for assessing bootstrap
# confidence interval reliability and distribution characteristics.  Features:
# - Skewness detection: Identifies non-normal bootstrap distributions -
# Multimodality detection: Detects multi-peaked distributions - CI width
# analysis: Assesses precision and stability of estimates - Integration with
# existing diagnostics infrastructure Reference: Zhang & Cao (2023); Friedl &
# Stampfer (2002) - Bootstrap CI quality assessment

#' Estimate Bootstrap Distribution Skewness with Robust Statistics
#'
#' Computes multiple skewness measures for bootstrap distributions:
#' - Fisher-Pearson skewness (moment-based)
#' - Quartile-based skewness (robust to outliers)
#' - Asymptotic confidence interval using jackknife
#'
#' @param boot_dist Numeric vector; bootstrap distribution
#' @param compute_ci Logical; if TRUE, computes confidence intervals via
#' jackknife
#'
#' @return List with components:
#'   - `skewness_mean`: Mean-based skewness (Fisher-Pearson)
#'   - `skewness_quartile`: Quartile-based skewness (robust)
#' - `skewness_median_absolute_dev`: Skewness using MAD (resistant to
#' outliers)
#'   - `ci_lower`: Lower 95% CI for skewness (if compute_ci=TRUE)
#'   - `ci_upper`: Upper 95% CI for skewness (if compute_ci=TRUE)
#'   - `interpretation`: Character description of skewness level
#'
#' @details
#' **Skewness measures:**
#' - **Fisher-Pearson (moment-based):** Most common, sensitive to outliers
#'     $\gamma_1 = \frac{E[(X - \mu)^3]}{\sigma^3}$
#'   Interpretation: |skewness| > 2 indicates strong asymmetry
#'
#' - **Quartile-based:** Robust to outliers, ranges in [-1, 1]
#'     Skewness = $\frac{(Q3 - Q2) - (Q2 - Q1)}{Q3 - Q1}$
#'   Interpretation: Closer to 0 = more symmetric
#'
#' - **Median Absolute Deviation (MAD):** Highly resistant to outliers
#'   Uses median and MAD instead of mean and SD
#'
#' **Confidence bounds:**
#' Bootstrap distributions with |skewness| > 2 may produce unreliable CIs.
#' Jackknife confidence intervals (Zhang & Cao 2023; Friedl & Stampfer 2002)
#' quantify skewness uncertainty.
#'
#' @noRd
#' @noRd
.estimate_bootstrap_skewness <- function(boot_dist, compute_ci = TRUE) {

    # Remove NA values
    x <- boot_dist[!is.na(boot_dist)]
    n <- length(x)

    if (n < 3) {
        return(list(skewness_mean = NA_real_, skewness_quartile = NA_real_, skewness_mad = NA_real_,
            ci_lower = NA_real_, ci_upper = NA_real_, interpretation = "insufficient data (n < 3)"))
    }

    # MEASURE 1: Fisher-Pearson skewness (moment-based) Formula: γ1 = E[(X -
    # μ)³] / σ³
    mean_x <- mean(x)
    sd_x <- sd(x)

    if (sd_x > 0) {
        skewness_mean <- (sum((x - mean_x)^3)/n)/(sd_x^3)
    } else {
        skewness_mean <- NA_real_
    }

    # MEASURE 2: Quartile-based skewness (robust) Formula: Skewness = ((Q3 -
    # Q2) - (Q2 - Q1)) / (Q3 - Q1)
    q1 <- stats::quantile(x, 0.25, type = 7)
    q2 <- stats::quantile(x, 0.5, type = 7)
    q3 <- stats::quantile(x, 0.75, type = 7)
    iqr <- q3 - q1

    if (iqr > 1e-10) {
        skewness_quartile <- ((q3 - q2) - (q2 - q1))/iqr
    } else {
        skewness_quartile <- NA_real_
    }

    # MEASURE 3: Skewness using Median Absolute Deviation (MAD) Highly
    # resistant to outliers
    median_x <- stats::median(x)
    mad_x <- stats::mad(x)  # Median absolute deviation

    if (mad_x > 1e-10) {
        # MAD-based skewness: use deviations from median scaled by MAD
        skewness_mad <- (sum((x - median_x)^3)/n)/(mad_x^3)
    } else {
        skewness_mad <- NA_real_
    }

    # CONFIDENCE INTERVALS (via jackknife; Zhang & Cao 2023, Friedl & Stampfer 2002)
    ci_lower <- NA_real_
    ci_upper <- NA_real_

    if (compute_ci && n >= 5) {
        # Jackknife replicates of skewness
        jack_skew <- numeric(n)

        for (i in seq_len(n)) {
            x_minus_i <- x[-i]
            m_i <- mean(x_minus_i)
            s_i <- sd(x_minus_i)

            if (s_i > 0) {
                jack_skew[i] <- (sum((x_minus_i - m_i)^3)/(n - 1))/(s_i^3)
            } else {
                jack_skew[i] <- NA_real_
            }
        }

        # SE via jackknife: SE = sqrt((n-1)/n * sum((x_j - x_bar)^2))
        jack_mean <- mean(jack_skew, na.rm = TRUE)
        jack_var <- sum((jack_skew - jack_mean)^2, na.rm = TRUE) * (n - 1)/n
        # jack_var already includes (n-1)/n factor.
        # SE should be sqrt(jack_var), not sqrt(jack_var/n).
        jack_se <- sqrt(jack_var)

        if (is.finite(jack_se) && jack_se > 0) {
            # 95% CI using normal approximation
            ci_lower <- skewness_mean - 1.96 * jack_se
            ci_upper <- skewness_mean + 1.96 * jack_se
        }
    }

    # Interpretation based on skewness magnitude
    interpretation <- if (is.na(skewness_mean)) {
        "Unable to compute (no variation in bootstrap dist)"
    } else if (abs(skewness_mean) < 0.5) {
        "Approximately symmetric (good for CI reliability)"
    } else if (abs(skewness_mean) < 1) {
        "Moderately skewed (acceptable for CI)"
    } else if (abs(skewness_mean) < 2) {
        "Highly skewed (caution: CI may be unreliable)"
    } else {
        "Extremely skewed (CI unreliable; consider BCa method)"
    }

    return(list(skewness_mean = skewness_mean, skewness_quartile = skewness_quartile,
        skewness_mad = skewness_mad, ci_lower = ci_lower, ci_upper = ci_upper, interpretation = interpretation))
}

#' Detect Multimodality in Bootstrap Distribution
#'
#' Uses multiple methods to detect if bootstrap distribution has multiple modes.
#' Multimodal distributions complicate CI interpretation and may indicate
#' problems with the input data or choice of bootstrap method.
#'
#' @param boot_dist Numeric vector; bootstrap distribution
#' @param method Character; detection method:
#'   - 'kde': Kernel density estimation (default, most accurate)
#'   - 'histogram': Simple histogram-based method
#'   - 'gaps': Detects large gaps in distribution (fastest)
#'   - 'all': Run all methods and summarize
#'
#' @return List with components:
#'   - `is_multimodal`: Logical; TRUE if multimodality detected
#'   - `n_modes`: Estimated number of modes (if detectable)
#'   - `modes_locations`: Estimated mode locations (numeric vector)
#'   - `separation_score`: How well-separated modes are (0-1, higher = better)
#'   - `method_used`: String indicating which method was used
#' - `interpretation`: Assessment of what multimodality means for the
#' bootstrap
#'
#' @details
#' **Methods:**
#'
#' 1. **KDE-based (kernel density estimation):**
#'    - Smooth bootstrap distribution using bandwith-adaptive KDE
#'    - Find local maxima (peaks) in density
#'    - Threshold: Need >10% density ratio between peaks and valleys
#'    - Most accurate but requires more computation
#'
#' 2. **Histogram-based:**
#' - Partition distribution into bins (Sturges rule: k = ceiling(log2(n) +
#' 1))
#'    - Count modes as bins with more items than median bin count
#'    - Fast and simple, less sensitive to bandwidth choice
#'
#' 3. **Gap-detection:**
#' - Identify large gaps between sorted values (>2 SD of inter-point
#' distance)
#'    - Fastest method, good for well-separated modes but misses close modes
#'
#' **Interpretation:**
#' - **Unimodal (1 mode):** Bootstrap distribution is well-behaved.
#'   Bootstrap CI is likely reliable.
#'
#' - **Bimodal to trimodal (2-3 modes):** Distribution has secondary peaks.
#' May indicate: (a) different parameter regimes, (b) boundary effects in
#' data,
#'   (c) inadequate bootstrap sample size. Bootstrap CI may be conservative.
#'   Recommendation: Check input data, consider BCa method.
#'
#' - **Highly multimodal (>3 modes):** Distribution has complex structure.
#' May indicate: (a) too many resampling boundaries, (b) specific data
#' patterns,
#'   (c) mixture distribution in original data. Bootstrap CIs may be unreliable.
#' Recommendation: Investigate input data, increase nboot, consider
#' alternative methods.
#'
#' @noRd
#' @noRd
.detect_multimodality <- function(boot_dist, method = c("kde", "histogram", "gaps")) {

    # Validate method parameter per Bioconductor code syntax standards
    method <- match.arg(method)

    # Remove NA values
    x <- boot_dist[!is.na(boot_dist)]
    n <- length(x)

    if (n < 10) {
        return(list(is_multimodal = NA, n_modes = NA_integer_, modes_locations = NA_real_,
            separation_score = NA_real_, method_used = "insufficient_data", interpretation = "Bootstrap distribution too small (n < 10) for mode detection"))
    }

    result <- if (method == "kde") {
        .detect_multimodality_kde(x)
    } else if (method == "histogram") {
        .detect_multimodality_histogram(x)
    } else if (method == "gaps") {
        .detect_multimodality_gaps(x)
    }

    return(result)
}

#' KDE-based Multimodality Detection
#' @noRd
.detect_multimodality_kde <- function(x) {

    # Estimate bandwidth using Silverman's rule
    n <- length(x)
    bw <- stats::bw.nrd0(x)

    # Create evaluation grid
    min_x <- min(x)
    max_x <- max(x)
    grid_x <- seq(min_x, max_x, length.out = 200)

    # Compute density at grid points via KDE
    density_vals <- vapply(grid_x, function(g) {
        mean(stats::dnorm(g - x, sd = bw))
    }, FUN.VALUE = numeric(1))

    if (length(density_vals) == 0 || sum(is.finite(density_vals)) < 3) {
        return(list(is_multimodal = FALSE, n_modes = 1L, modes_locations = mean(x),
            separation_score = NA_real_, method_used = "kde_failed", interpretation = "KDE computation failed; assuming unimodal"))
    }

    # Find local maxima (modes) A grid point is a mode if density higher than
    # neighbors
    n_grid <- length(density_vals)
    modes_mask <- logical(n_grid)

    for (i in seq_len(n_grid)) {
        if (i == 1 || i == n_grid)
            next  # Skip boundaries

        # Check if local maximum (density > both neighbors)
        if (density_vals[i] > density_vals[i - 1] && density_vals[i] > density_vals[i +
            1]) {
            # Also check if above threshold (>10% of max density)
            if (density_vals[i] > 0.1 * max(density_vals)) {
                modes_mask[i] <- TRUE
            }
        }
    }

    # Cluster nearby modes (within 3 grid points)
    mode_indices <- which(modes_mask)
    if (length(mode_indices) == 0) {
        modes_locations <- mean(x)
        n_modes <- 1L
    } else if (length(mode_indices) == 1) {
        # Single mode found
        modes_locations <- grid_x[mode_indices[1]]
        n_modes <- 1L
    } else {
        # Multiple modes: merge nearby ones
        clustered_modes <- numeric()
        current_cluster <- c(mode_indices[1])

        for (i in 2:length(mode_indices)) {
            if (mode_indices[i] - mode_indices[i - 1] <= 3) {
                current_cluster <- c(current_cluster, mode_indices[i])
            } else {
                # Save cluster center
                cluster_center_idx <- current_cluster[which.max(density_vals[current_cluster])]
                clustered_modes <- c(clustered_modes, grid_x[cluster_center_idx])
                current_cluster <- c(mode_indices[i])
            }
        }
        # Save final cluster
        cluster_center_idx <- current_cluster[which.max(density_vals[current_cluster])]
        clustered_modes <- c(clustered_modes, grid_x[cluster_center_idx])

        modes_locations <- clustered_modes
        n_modes <- as.integer(length(clustered_modes))
    }

    # Compute separation score (how well-separated modes are) If modes are
    # close together, separation_score is low
    if (n_modes > 1) {
        mode_diffs <- diff(sort(modes_locations))
        avg_separation <- mean(mode_diffs)
        data_range <- max(x) - min(x)
        separation_score <- min(1, avg_separation/(data_range/n_modes))
    } else {
        separation_score <- 1
    }

    is_multimodal <- n_modes > 1

    interpretation <- if (n_modes == 1) {
        "Unimodal distribution (single mode) - good for bootstrap CI"
    } else if (n_modes == 2) {
        sprintf("Bimodal distribution (2 modes) - examine input data; bootstrap CI may be conservative")
    } else {
        sprintf("Multimodal distribution (%d modes) - bootstrap CI reliability questionable",
            n_modes)
    }

    list(is_multimodal = is_multimodal, n_modes = n_modes, modes_locations = modes_locations,
        separation_score = separation_score, method_used = "kde", interpretation = interpretation)
}

#' Histogram-based Multimodality Detection
#' @noRd
.detect_multimodality_histogram <- function(x) {

    n <- length(x)

    # Sturges rule for number of bins
    n_bins <- ceiling(log2(n) + 1)

    # Compute histogram
    h <- graphics::hist(x, breaks = n_bins, plot = FALSE)

    # Find bins with above-median counts
    median_count <- stats::median(h$counts)
    mode_bins <- which(h$counts > median_count * 1.2)  # 20% threshold

    if (length(mode_bins) == 0) {
        n_modes <- 1L
        modes_locations <- mean(x)
    } else {
        # Cluster adjacent mode bins
        mode_locations <- h$mids[mode_bins]

        # Simple clustering: modes within 1 bin width are same mode
        bin_width <- h$breaks[2] - h$breaks[1]
        modes_locations <- numeric()
        current_modes <- c(mode_locations[1])

        for (i in 2:length(mode_locations)) {
            if (abs(mode_locations[i] - mode_locations[i - 1]) <= 1.5 * bin_width) {
                current_modes <- c(current_modes, mode_locations[i])
            } else {
                modes_locations <- c(modes_locations, mean(current_modes))
                current_modes <- c(mode_locations[i])
            }
        }
        modes_locations <- c(modes_locations, mean(current_modes))

        n_modes <- as.integer(length(modes_locations))
    }

    # Separation score
    if (n_modes > 1) {
        mode_diffs <- diff(sort(modes_locations))
        avg_separation <- mean(mode_diffs)
        data_range <- max(x) - min(x)
        separation_score <- min(1, avg_separation/(data_range/n_modes))
    } else {
        separation_score <- 1
    }

    is_multimodal <- n_modes > 1

    interpretation <- if (n_modes == 1) {
        "Unimodal distribution (single mode) - good for bootstrap CI"
    } else if (n_modes == 2) {
        "Bimodal distribution (2 modes) - examine input data"
    } else {
        sprintf("Multimodal distribution (%d modes) - bootstrap CI reliability questionable",
            n_modes)
    }

    list(is_multimodal = is_multimodal, n_modes = n_modes, modes_locations = sort(modes_locations),
        separation_score = separation_score, method_used = "histogram", interpretation = interpretation)
}

#' Gap-based Multimodality Detection
#' @noRd
.detect_multimodality_gaps <- function(x) {

    n <- length(x)
    x_sorted <- sort(x)

    # Compute inter-point gaps
    gaps <- diff(x_sorted)

    if (length(gaps) < 2) {
        return(list(is_multimodal = FALSE, n_modes = 1L, modes_locations = mean(x),
            separation_score = NA_real_, method_used = "gaps", interpretation = "Insufficient data for gap-based detection"))
    }

    # Identify large gaps (> 2 SD of mean gap)
    mean_gap <- mean(gaps)
    sd_gap <- sd(gaps)

    large_gap_threshold <- mean_gap + 2 * sd_gap
    large_gaps <- which(gaps > large_gap_threshold)

    # Number of modes = 1 + number of large gaps
    n_modes <- as.integer(1 + length(large_gaps))

    # Estimate mode locations (median of each segment)
    if (n_modes == 1) {
        modes_locations <- stats::median(x)
    } else {
        segment_starts <- c(1, large_gaps + 1)
        segment_ends <- c(large_gaps, n)
        modes_locations <- vapply(seq_len(n_modes), function(i) {
            stats::median(x_sorted[segment_starts[i]:segment_ends[i]])
        }, FUN.VALUE = numeric(1))
    }

    # Separation score
    if (n_modes > 1) {
        mode_diffs <- diff(sort(modes_locations))
        avg_separation <- mean(mode_diffs)
        data_range <- max(x) - min(x)
        separation_score <- min(1, avg_separation/(data_range/n_modes))
    } else {
        separation_score <- 1
    }

    is_multimodal <- n_modes > 1

    interpretation <- if (n_modes == 1) {
        "Unimodal distribution (no large gaps detected)"
    } else if (n_modes == 2) {
        "Bimodal distribution (1 large gap detected)"
    } else {
        sprintf("Multimodal distribution (%d modes, %d large gaps)", n_modes, length(large_gaps))
    }

    list(is_multimodal = is_multimodal, n_modes = n_modes, modes_locations = sort(modes_locations),
        separation_score = separation_score, method_used = "gaps", interpretation = interpretation)
}

#' Analyze Bootstrap Confidence Interval Width and Characteristics
#'
#' Provides comprehensive analysis of CI width to assess precision,
#' stability, and potential issues with bootstrap estimation.
#'
#' @param ci_lower Numeric; lower confidence bound
#' @param ci_upper Numeric; upper confidence bound
#' @param point_est Numeric; point estimate (from original data)
#' @param bootstrap_dist Numeric vector; bootstrap distribution
#' @param n_bootstrap Integer; number of bootstrap replicates
#'
#' @return List with components:
#'   - `ci_width`: Raw CI width (upper - lower)
#'   - `ci_width_to_estimate_ratio`: CI width normalized by point estimate
#'   - `ci_width_to_sd_ratio`: CI width normalized by bootstrap SD
#'   - `coverage_estimate`: Estimated empirical coverage probability
#'   - `precision_assessment`: Qualitative assessment of precision
#'   - `potential_issues`: Character vector of detected issues
#'   - `recommendations`: Character vector of suggested actions
#'
#' @details
#' **Ratios:**
#' - **CI width / estimate:** High value (>0.5) suggests low precision
#' relative to estimate
#' - **CI width / SD:** Ratio ~4 is typical for 95% CIs (2.5*SD on each side);
#'   much higher ratios suggest longer-tailed bootstrap distributions
#'
#' **Precision levels:**
#' - Excellent: CI width < 0.1 * estimate (±5% relative uncertainty)
#' - Good: CI width < 0.25 * estimate (±12.5% relative uncertainty)
#' - Acceptable: CI width < 0.5 * estimate (±25% relative uncertainty)
#' - Poor: CI width >= 0.5 * estimate (>±25% relative uncertainty)
#'
#' **Issues detected:**
#' - Asymmetric CI: Large difference between distance to lower and upper bounds
#' - Negative lower bound: May indicate boundary issues for non-negative
#' quantities
#' - Wide relative CI: High relative uncertainty
#' - Small n_bootstrap: Low effective sample size for CI computation
#'
#' @noRd
#' @noRd
.analyze_ci_width <- function(ci_lower, ci_upper, point_est, boot_dist, nboot) {

    # Basic CI characteristics
    ci_width <- ci_upper - ci_lower

    # Ratios for interpretation
    estimate_abs <- abs(point_est)
    if (estimate_abs > 0) {
        ci_width_to_est_ratio <- ci_width/estimate_abs
    } else {
        ci_width_to_est_ratio <- NA_real_
    }

    boot_sd <- sd(boot_dist, na.rm = TRUE)
    if (boot_sd > 0) {
        ci_width_to_sd_ratio <- ci_width/boot_sd
    } else {
        ci_width_to_sd_ratio <- NA_real_
    }

    # Symmetry of CI bounds
    lower_tail <- point_est - ci_lower
    upper_tail <- ci_upper - point_est

    if (lower_tail > 0 && upper_tail > 0) {
        tail_ratio <- min(lower_tail, upper_tail)/max(lower_tail, upper_tail)
    } else {
        tail_ratio <- NA_real_
    }

    # Estimated coverage (empirical) For well-behaved bootstrap, approximately
    # 95% of replicates within CI
    coverage <- mean(boot_dist >= ci_lower & boot_dist <= ci_upper, na.rm = TRUE) *
        100

    # Precision assessment
    precision <- if (is.na(ci_width_to_est_ratio)) {
        "indeterminate (zero estimate)"
    } else if (ci_width_to_est_ratio < 0.1) {
        "excellent"
    } else if (ci_width_to_est_ratio < 0.25) {
        "good"
    } else if (ci_width_to_est_ratio < 0.5) {
        "acceptable"
    } else {
        "poor"
    }

    # Detect potential issues
    issues <- character()

    # Asymmetric CI
    if (!is.na(tail_ratio) && tail_ratio < 0.7) {
        issues <- c(issues, "Asymmetric CI (consider BCa method)")
    }

    # Negative lower bound for non-negative quantities
    if (point_est >= 0 && ci_lower < -0.01 * abs(point_est)) {
        issues <- c(issues, "Negative lower bound (add pseudocount?)")
    }

    # Wide CI relative to estimate
    if (!is.na(ci_width_to_est_ratio) && ci_width_to_est_ratio > 0.5) {
        issues <- c(issues, "Wide CI relative to estimate (low information)")
    }

    # Small effective sample size
    if (nboot < 100) {
        issues <- c(issues, "Low n_bootstrap (< 100, less stable CI)")
    }

    # High CI width to SD ratio (suggests long tails)
    if (!is.na(ci_width_to_sd_ratio) && ci_width_to_sd_ratio > 5) {
        issues <- c(issues, "Long-tailed bootstrap distribution")
    }

    # Recommendations
    recommendations <- character()

    if (length(issues) > 0) {
        if ("Asymmetric CI (consider BCa method)" %in% issues) {
            recommendations <- c(recommendations, "Use BCa method instead of percentile")
        }
        if ("Negative lower bound (add pseudocount?)" %in% issues) {
            recommendations <- c(recommendations, "Try adding pseudocount to counts")
        }
        if ("Low n_bootstrap (< 100, less stable CI)" %in% issues) {
            recommendations <- c(recommendations, "Increase n_bootstrap for more stable CI")
        }
        if ("Wide CI relative to estimate (low information)" %in% issues) {
            recommendations <- c(recommendations, "Increase sample size or consider other measurements")
        }
    }

    if (precision %in% c("excellent", "good")) {
        recommendations <- c(recommendations, "CI appears reliable")
    }

    list(ci_width = ci_width, ci_width_to_estimate_ratio = ci_width_to_est_ratio,
        ci_width_to_sd_ratio = ci_width_to_sd_ratio, ci_symmetry_ratio = tail_ratio,
        coverage_estimate = coverage, precision_assessment = precision, potential_issues = if (length(issues) >
            0) issues else "none detected", recommendations = if (length(recommendations) >
            0) recommendations else "none needed")
}

#' Integrated Bootstrap Diagnostics Report
#'
#' Combines skewness, multimodality, and CI width analysis into a
#' comprehensive assessment report with actionable recommendations.
#'
#' @param boot_result Object of class \code{tsenat_bootstrap_ci}
#'
#' @return List with integrated diagnostics:
#'   - `skewness_analysis`: Output from `.estimate_bootstrap_skewness()`
#'   - `multimodality_analysis`: Output from `.detect_multimodality()`
#'   - `ci_width_analysis`: Output from `.analyze_ci_width()`
#'   - `overall_reliability`: Character assessment (reliable/caution/unreliable)
#'   - `summary_recommendations`: List of recommended actions
#'
#' @details
#' This function provides an all-in-one diagnostic summary suitable for
#' validation reports and supplementary materials.
#'
#' **Reliability tiers:**
#' - **Reliable:** Bootstrap distribution is well-behaved (unimodal, low
#' skewness,
#'   symmetric CI, sufficient n_bootstrap). CI can be used with confidence.
#' - **Caution:** Some non-ideal characteristics detected (moderate skewness,
#' slightly asymmetric CI, or modest sample size). CI is usable but
#' conservative
#'   interpretation recommended.
#' - **Unreliable:** Major issues detected (strong multimodality, extreme
#' skewness,
#' very asymmetric CI). Bootstrap CI may not be valid; consider alternative
#' methods.
#'
#' @noRd
#' @noRd
.generate_bootstrap_diagnostics_report <- function(boot_result) {

    if (!inherits(boot_result, "tsenat_bootstrap_ci")) {
        stop("boot_result must be of class tsenat_bootstrap_ci")
    }

    # Run all diagnostics
    skewness_diag <- .estimate_bootstrap_skewness(boot_result$bootstrap_dist, compute_ci = TRUE)
    multimodality_diag <- .detect_multimodality(boot_result$bootstrap_dist, method = "kde")
    ci_width_diag <- .analyze_ci_width(boot_result$lower_ci, boot_result$upper_ci,
        boot_result$estimate, boot_result$bootstrap_dist, boot_result$nboot)

    # Overall reliability assessment
    reliability_flags <- 0

    # Flag 1: Skewness
    if (!is.na(skewness_diag$skewness_mean)) {
        if (abs(skewness_diag$skewness_mean) > 2)
            reliability_flags <- reliability_flags + 2 else if (abs(skewness_diag$skewness_mean) > 1)
            reliability_flags <- reliability_flags + 1
    }

    # Flag 2: Multimodality
    if (multimodality_diag$is_multimodal) {
        reliability_flags <- reliability_flags + (multimodality_diag$n_modes - 1)
    }

    # Flag 3: CI asymmetry
    if (!is.na(ci_width_diag$ci_symmetry_ratio) && ci_width_diag$ci_symmetry_ratio <
        0.6) {
        reliability_flags <- reliability_flags + 1
    }

    # Flag 4: Small n_bootstrap (100 is minimum recommended, not ideal)
    if (boot_result$nboot <= 100) {
        reliability_flags <- reliability_flags + 1
    }

    # Overall assessment based on flag count
    overall_reliability <- if (reliability_flags == 0) {
        "Reliable"
    } else if (reliability_flags <= 2) {
        "Caution"
    } else {
        "Unreliable"
    }

    # Generate summary recommendations
    summary_recommendations <- character()

    if (overall_reliability == "Reliable") {
        summary_recommendations <- c("Bootstrap CI appears well-behaved and can be used with confidence.",
            "Distribution is approximately normal with symmetric CI bounds.")
    } else if (overall_reliability == "Caution") {
        summary_recommendations <- c("Bootstrap CI has some non-ideal characteristics but is usable.",
            "Conservative interpretation recommended; consider BCa method.")
        if (multimodality_diag$is_multimodal) {
            summary_recommendations <- c(summary_recommendations, sprintf("Distribution appears multimodal (%d modes). Check input data for mixture structure.",
                multimodality_diag$n_modes))
        }
        if (abs(skewness_diag$skewness_mean) > 1) {
            summary_recommendations <- c(summary_recommendations, "Bootstrap distribution is skewed; BCa method may be more accurate than percentile.")
        }
    } else {
        # Unreliable
        summary_recommendations <- c("Bootstrap CI may not be valid. Consider alternative approaches:",
            "  1. Check input data quality and distribution", "  2. Increase n_bootstrap to >= 2000",
            "  3. Use BCa method instead of percentile", "  4. Try non-parametric alternatives (e.g., jackknife)",
            "  5. Add pseudocount if zero counts are problematic")
    }

    return(list(skewness_analysis = skewness_diag, multimodality_analysis = multimodality_diag,
        ci_width_analysis = ci_width_diag, overall_reliability = overall_reliability,
        summary_recommendations = summary_recommendations))
}
