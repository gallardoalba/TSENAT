# =========================================================================
# PRIVATE HELPER FUNCTIONS (must come before main roxygen block)
# =========================================================================

# -------------------------------------------------------------------------
# PRIVATE HELPER: Auto-detect group column and control group
# -------------------------------------------------------------------------

#' Auto-Detect Group Column and Control Group from SummarizedExperiment colData
#'
#' Automatically searches for a column containing group/condition assignments
#' and identifies the reference (control) group.
#'
#' **Detection Strategy for group_col:**
#' Searches colData in precedence order for common grouping column names:
#' 1. 'group' (TSENAT default)
#' 2. 'condition', 'treatment', 'sample_type' (common alternatives)
#'
#' **Detection Strategy for control_group:**
#' Once group column is found, identifies the reference group by looking for
#' common control/reference names: 'Normal', 'Control', 'WT' (wild-type),
#' 'Reference', 'Baseline', 'Wild-type', 'wild_type', 'wt'.
#' If no standard name matches, an error is raised because the control group
#' is a scientific decision that cannot be reliably guessed --- the user must
#' specify it explicitly via the \code{control_group} parameter.
#'
#' @param se SummarizedExperiment object with sample metadata in colData
#'
#' @return List with elements:
#'   \describe{
#'     \item{group_col}{Name of the colData column used for  grouping,  or 
#' NA_character_ if  none detected}
#'     \item{control_group}{Name of the control/reference group,  or 
#' NA_character_ if  none detected}
#'     \item{groups}{Character vector of all unique groups found}
#'     \item{sample_counts}{Named integer vector of sample counts per group (names:
#'  group names)}
#'   }
#'

#' @noRd
.auto_detect_groups <- function(se) {

    cd <- SummarizedExperiment::colData(se)
    cd_colnames <- colnames(cd)

    # Candidate column names (in priority order) for group/condition.
    # We intentionally restrict this list to likely biological/experimental
    # grouping variables, not technical factors such as batch identifiers.
    group_col_candidates <- c("sample_type", "group", "condition", "treatment",
        "phenotype")

    group_col <- NA_character_

    # Find first matching group column
    for (col_name in group_col_candidates) {
        if (col_name %in% cd_colnames) {
            group_col <- col_name
            break
        }
    }

    # If no group column found, return NAs
    if (is.na(group_col)) {
        return(list(group_col = NA_character_, control_group = NA_character_, groups = character(0),
            sample_counts = integer(0)))
    }

    # Get group values and counts
    group_vec <- as.character(cd[[group_col]])
    unique_groups <- unique(group_vec)
    group_counts <- table(group_vec)

    # Candidate names for control/reference group (in priority order).
    # Both capitalized and lowercase variants are included because real-world
    # metadata uses inconsistent casing (e.g., "normal" in vignette, "Control"
    # in TCGA-style datasets).
    control_candidates <- c("Normal", "normal", "Control", "control", "WT", "wt",
        "Reference", "reference", "Baseline", "baseline", "Wild-type", "wild_type")

    control_group <- NA_character_

    # Try to match control candidates
    for (control_name in control_candidates) {
        if (control_name %in% unique_groups) {
            control_group <- control_name
            break
        }
    }

    # If no standard control label is found, we cannot reliably guess the
    # reference group.  The control group is a scientific decision, not a
    # computational convenience --- choosing the smallest group, the first
    # alphabetically, or any other heuristic can silently produce incorrect
    # results.  Per Bioconductor reproducibility guidelines, we error out
    # with a clear message that lists the available groups so the user can
    # make an explicit, documented choice.
    if (is.na(control_group)) {
        if (length(unique_groups) >= 2) {
            stop(
                "Could not auto-detect control_group. ",
                "Available groups: ", paste(sQuote(unique_groups), collapse = ", "), ". ",
                "Please specify 'control_group' explicitly (e.g., control_group = \"",
                unique_groups[1], "\").",
                call. = FALSE
            )
        } else if (length(unique_groups) == 1) {
            control_group <- unique_groups[1]
        }
    }

    return(list(group_col = group_col, control_group = control_group, groups = unique_groups,
        sample_counts = as.vector(group_counts)))
}


# -------------------------------------------------------------------------
# PRIVATE HELPER: Auto-detect paired samples from metadata
# -------------------------------------------------------------------------

#' Detect Paired Sample Structure from SummarizedExperiment colData
#'
#' Automatically searches for a column containing paired sample identifiers
#' (e.g., 'pair_id', 'paired_samples', 'patient_id', 'subject_id').
#' This enables pair-respecting bootstrap resampling in divergence calculations.
#'
#' **Detection Strategy:**
#' Searches colData in precedence order for common pairing column names:
#' 1. 'paired_samples' (TSENAT default, matches readcounts metadata)
#' 2. 'pair_id', 'pair_samples', 'subject_id', 'patient_id' (common
#' alternatives)
#'
#' Returns a mapping from sample names to pair identifiers, or NULL if no
#' pairing column is found. A valid pairing column has:
#' - Non-NA values for all samples
#' - At least 2 samples per pair
#' - Deterministic structure (e.g., all A's paired with another A sample nearby)
#'
#' Pairing auto-detection follows standard paired-design conventions:
#' subject/patient identifiers are shared across conditions.
#'
#' @param se SummarizedExperiment object with sample metadata in colData
#'
#' @return List with elements:
#'   \describe{
#'     \item{pair_ids}{Character vector (names:  sample names,  values:
#'  pair identifiers)
#'       or NULL if no pairing detected}
#'     \item{column_name}{Name of the colData column used,  or 
#' NA_character_ if  none}
#'     \item{num_pairs}{Number of unique pairs (0 if none detected)}
#'     \item{samples_per_pair}{Vector of samples per pair (names:  pair IDs,
#'  values:  counts)}
#'   }
#'
#' @note Paired samples detected from any of: 'paired_samples', 'pair_id',
#' 'pair_samples',
#' 'subject_id', 'patient_id'. Returns NULL if none present or validation
#' fails.
#'

#' @noRd
.detect_pair_ids <- function(se) {

    cd <- SummarizedExperiment::colData(se)
    sample_names <- colnames(se)

    # Candidate column names (in priority order)
    candidate_cols <- c("paired_samples", "pair_id", "pair_samples", "subject_id",
        "patient_id")

    for (col_name in candidate_cols) {
        if (col_name %in% colnames(cd)) {
            pair_col <- cd[[col_name]]

            # Validate: must be non-NA for all samples
            if (any(is.na(pair_col))) {
                next  # Skip if any NAs
            }

            # Valid pairing structure found Keep pair_ids as character (names
            # are character in the data)
            pair_ids <- setNames(as.character(pair_col), sample_names)
            unique_pairs <- unique(pair_ids)
            samples_per_pair <- table(pair_ids)

            # Structural validation: a pairing column must
            # actually PAIR observations. A column in which every identifier
            # appears once (e.g., a plain sample/subject ID) is not a pair
            # structure; skip it and try the next candidate.
            if (sum(samples_per_pair >= 2) == 0) {
                next
            }

            return(list(pair_ids = pair_ids, column_name = col_name, num_pairs = length(unique_pairs),
                samples_per_pair = samples_per_pair))
        }
    }

    # No pairing detected
    return(list(pair_ids = NULL, column_name = NA_character_, num_pairs = 0, samples_per_pair = numeric(0)))
}


#' Validate the paired-design invariant: at most 1 control + 1 treatment per pair
#'
#' AUDIT R5: pair-respecting bootstrap resampling requires a well-defined
#' resampling unit. A pair that mixes both conditions with REPEATED
#' observations in one of them (e.g., control + control + treatment) is
#' rejected: its resampling unit is ambiguous. Incomplete/single-condition
#' pairs are tolerated because the bootstrap machinery handles unmatched
#' samples as unpaired units. Validate this BEFORE bootstrap; stop with an
#' actionable message instead of silently skipping replicates.
#'
#' @param se SummarizedExperiment with sample metadata in colData
#' @param pair_ids Named character vector (names = sample names, values = pair IDs)
#' @param group_col character; colData column with group/condition labels
#' @param control_group character; control group label
#'
#' @return invisible(TRUE) if the invariant holds; errors otherwise
#' @noRd
.validate_pair_structure <- function(se, pair_ids, group_col, control_group) {
    cd <- SummarizedExperiment::colData(se)
    if (is.null(pair_ids) || length(pair_ids) == 0) {
        return(invisible(TRUE))
    }
    if (is.null(group_col) || !group_col %in% colnames(cd)) {
        return(invisible(TRUE))  # No group info: cannot validate, stay permissive
    }

    groups <- as.character(cd[[group_col]])
    names(groups) <- colnames(se)

    bad_pairs <- character(0)
    for (pid in unique(as.character(pair_ids))) {
        idx <- which(as.character(pair_ids) == pid)
        smp <- names(pair_ids)[idx]
        if (length(smp) == 0) next
        g <- groups[smp]
        n_ctrl <- sum(g == control_group, na.rm = TRUE)
        n_trt <- sum(g != control_group, na.rm = TRUE)
        # Reject only pairs that mix both conditions AND contain repeated
        # observations within a condition (ill-defined resampling unit).
        # Incomplete or single-condition pairs are handled as unpaired units
        # by the bootstrap machinery.
        if (n_ctrl >= 1 && n_trt >= 1 && (n_ctrl > 1 || n_trt > 1)) {
            bad_pairs <- c(bad_pairs, pid)
        }
    }

    if (length(bad_pairs) > 0) {
        stop("[.validate_pair_structure] Paired-design invariant violated: pair(s) ",
            paste(sQuote(unique(bad_pairs)), collapse = ", "),
            " contain repeated samples from the same condition (each pair may contain at most one control and one ",
            "treatment sample; the pair is the resampling unit). ",
            "Fix the pairing column in colData or set bootstrap = FALSE / paired = FALSE.",
            call. = FALSE)
    }

    invisible(TRUE)
}


#' Resample Data Respecting Paired Structure
#'
#' When resampling paired data, both members of a pair are selected or discarded
#' together. This preserves within-pair correlations critical for
#' statistical validity
#' in matched designs (Efron & Tibshirani 1993).
#'
#' For k pairs:
#' - Draw k pair indices uniformly with replacement: pair_idx ~ U(1:k)
#' - For each drawn pair i, include both samples from pair i
#' - This maintains pairing structure across bootstrap replicates
#'
#' **Statistical Justification:**
#' - Ramsay (2005), Springer Series in Statistics: Bootstrap for confidence intervals requires preserving data structure
#' - Paired resampling is required for matched samples; within-pair correlation invalidates
#' independent resampling
#'
#' @param control_samples Vector of counts for control group
#' @param treatment_samples Vector of counts for treatment group
#' @param pair_ids Named character vector: names = sample names, values =
#' pair IDs
#' @param group_col Character vector: group assignment (one per sample)
#' @param control_group Character: name of control group
#'
#' @return List with elements:
#'   \describe{
#'     \item{control_resampled}{Resampled control group counts}
#'     \item{treatment_resampled}{Resampled treatment group counts}
#'   }
#'

#' @noRd
.jis_resample_paired_data <- function(control_samples, treatment_samples, pair_ids,
    group_col, control_group) {

    # Map sample names to indices
    all_samples <- c(names(control_samples), names(treatment_samples))
    all_groups <- c(rep(control_group, length(control_samples)), rep(setdiff(unique(group_col),
        control_group), length(treatment_samples)))

    # Get unique pairs involved
    pairs_in_data <- unique(pair_ids[all_samples])
    num_pairs <- length(pairs_in_data)

    # Resample pairs with replacement
    resampled_pairs <- sample(seq_len(num_pairs), size = num_pairs, replace = TRUE)
    resampled_pair_ids <- pairs_in_data[resampled_pairs]

    # Collect samples for each resampled pair (pre-allocate for O(n) performance)
    n_resampled <- length(resampled_pair_ids) * 2L  # upper bound: 2 samples per pair
    resampled_control <- numeric(n_resampled)
    resampled_treatment <- numeric(n_resampled)
    c_idx <- 1L
    t_idx <- 1L

    for (pair_id in resampled_pair_ids) {
        # Get both samples from this pair
        pair_samples <- names(pair_ids)[pair_ids == pair_id]

        for (sample_name in pair_samples) {
            if (sample_name %in% names(control_samples)) {
                resampled_control[c_idx] <- control_samples[sample_name]
                c_idx <- c_idx + 1L
            } else if (sample_name %in% names(treatment_samples)) {
                resampled_treatment[t_idx] <- treatment_samples[sample_name]
                t_idx <- t_idx + 1L
            }
        }
    }
    # Trim to actual used length
    resampled_control <- resampled_control[seq_len(c_idx - 1L)]
    resampled_treatment <- resampled_treatment[seq_len(t_idx - 1L)]

    # BUGFIX #3: Validate balanced groups after paired resampling Ensures
    # control and treatment have equal sizes (required for divergence
    # computation)
    if (length(resampled_control) != length(resampled_treatment)) {
        warning("Paired bootstrap produced unequal group sizes (", length(resampled_control),
            " control vs ", length(resampled_treatment), " treatment). ",
            "Skipping this bootstrap replicate. Check for unbalanced or incomplete pairs in input data.")
        return(list(control_resampled = numeric(0), treatment_resampled = numeric(0),
            failed = TRUE))
    }

    return(list(control_resampled = resampled_control, treatment_resampled = resampled_treatment,
        failed = FALSE))
}





# -------------------------------------------------------------------------
# HELPER FUNCTION: Bootstrap CI Computation
# -------------------------------------------------------------------------

#' Calculate Tsallis Divergence Bootstrap Confidence Intervals
#'
#' Internal helper function for compute_divergence_bootstrap.
#' Computes bootstrap confidence intervals for Tsallis divergence.
#'

#' @noRd
.prepare_paired_bootstrap_data <- function(x, y, pair_ids) {
    # Extract separate pair_ids for x and y, ready for flexible C++ bootstrap
    # OPTIMIZATION (March 2026): Use enhanced C++ supporting mixed
    # paired/unpaired This version handles complete pairs, incomplete pairs,
    # and unpaired samples No R fallback - all cases use C++ (~10x speedup)
    # Returns: List with components: $x_pair_ids: Pair IDs for x samples (0 =
    # unpaired) $y_pair_ids: Pair IDs for y samples (0 = unpaired) $valid: TRUE
    # if extraction successful, FALSE otherwise

    tryCatch({
        # Step 1: Get sample names from x and y
        x_names <- names(x)
        y_names <- names(y)

        # Handle NAs in pair_ids before any operations Accept both integer and
        # numeric vectors with names (required for matching)
        if (is.null(names(pair_ids))) {
            # pair_ids must have names to match with x and y samples
            warning("pair_ids has no names --- cannot perform paired resampling. Falling back to unpaired bootstrap.")
            return(list(valid = FALSE))
        }

        if (!(is.numeric(pair_ids) || is.integer(pair_ids) || is.character(pair_ids))) {
            warning("pair_ids must be a named integer/numeric/character vector. Falling back to unpaired bootstrap.")
            return(list(valid = FALSE))
        }

        # Save original names before any conversion
        pair_ids_names <- names(pair_ids)

        # If pair_ids are character, convert to numeric indices while
        # preserving names
        if (is.character(pair_ids)) {
            unique_pair_values <- unique(pair_ids)
            pair_id_map <- setNames(seq_along(unique_pair_values), unique_pair_values)
            pair_ids_numeric <- as.numeric(pair_id_map[pair_ids])
            pair_ids <- setNames(pair_ids_numeric, pair_ids_names)
        } else {
            # Ensure pair_ids is numeric for consistent handling
            pair_ids <- as.numeric(pair_ids)
            names(pair_ids) <- pair_ids_names  # Restore names that may be lost in conversion
        }

        # Replace NAs with 0 for unpaired samples
        pair_ids[is.na(pair_ids)] <- 0

        if (is.null(x_names) || is.null(y_names)) {
            warning("x and y must have names for paired bootstrap pairing.")
            return(list(valid = FALSE))
        }

        if (length(x_names) == 0 || length(y_names) == 0) {
            warning("x and y have empty names.")
            return(list(valid = FALSE))
        }

        # Step 2: Extract pair_ids for x samples
        if (!all(x_names %in% names(pair_ids))) {
            warning("Not all x samples found in pair_ids. Cannot prepare paired bootstrap.")
            return(list(valid = FALSE))
        }

        x_pair_ids <- pair_ids[x_names]

        # Step 3: Extract pair_ids for y samples
        if (!all(y_names %in% names(pair_ids))) {
            warning("Not all y samples found in pair_ids. Cannot prepare paired bootstrap.")
            return(list(valid = FALSE))
        }

        y_pair_ids <- pair_ids[y_names]

        # Step 4: Identify pairing structure for summary
        x_paired_mask <- x_pair_ids > 0
        y_paired_mask <- y_pair_ids > 0

        x_paired_count <- sum(x_paired_mask)
        y_paired_count <- sum(y_paired_mask)
        x_unpaired_count <- sum(!x_paired_mask)
        y_unpaired_count <- sum(!y_paired_mask)

        # Log pairing structure
        if (length(pair_ids) > 0) {
            # Count complete pairs (pair_id in both x and y)
            x_pair_set <- setNames(x_pair_ids[x_paired_mask], NULL)
            y_pair_set <- setNames(y_pair_ids[y_paired_mask], NULL)
            complete_pair_ids <- intersect(unique(x_pair_set[x_pair_set > 0]), unique(y_pair_set[y_pair_set >
                0]))

            if (length(complete_pair_ids) > 0) {
                message(sprintf("Paired bootstrap structure: %d complete pairs, %d unpaired x, %d unpaired y",
                  length(complete_pair_ids), x_unpaired_count, y_unpaired_count),
                  domain = NA)
            }
        }

        # Return extracted pair_ids as integers (no NAs, safe for C++)
        x_result <- as.numeric(pair_ids[x_names])
        y_result <- as.numeric(pair_ids[y_names])
        x_result[is.na(x_result)] <- 0
        y_result[is.na(y_result)] <- 0

        return(list(x_pair_ids = x_result, y_pair_ids = y_result, valid = TRUE))

    }, error = function(e) {
        warning("Error preparing paired bootstrap data: ", e$message)
        return(list(valid = FALSE))
    })
}

#' @noRd

.calculate_divergence_bootstrap <- function(x, y, q = 1, nboot = 1000, ci = 0.95,
    method = "percentile", log_base = exp(1), pseudocount = 0.5, gene_name = NA_character_,
    verbose = FALSE, paired = FALSE, pair_ids = NULL) {

    # Seed handling left to caller for Bioconductor compliance

    # Compute point estimate
    point_est <- .tsallis_divergence_scalar(x, y, q, pseudocount, log_base)

    # Bootstrap confidence interval
    if (nboot > 0) {
        # Determine if using pair_ids-based pairing
        use_complex_paired_bootstrap <- isTRUE(paired) && !is.null(pair_ids)

        if (use_complex_paired_bootstrap) {
            # OPTIMIZATION (March 2026): Use C++ flexible paired bootstrap
            # Supports mixed paired/unpaired data, handles sparse pairings
            pair_data <- .prepare_paired_bootstrap_data(x, y, pair_ids)

            if (isTRUE(pair_data$valid)) {
                # Use enhanced C++ implementation for all pairing scenarios No
                # fallback - all cases handled by C++ (~10x speedup)
                bootstrap_dist <- tryCatch({
                  divergence_bootstrap_flexible_cpp_wrapper(x = as.numeric(x), y = as.numeric(y),
                    x_pair_ids = pair_data$x_pair_ids, y_pair_ids = pair_data$y_pair_ids,
                    nboot = as.integer(nboot), q = q, pseudocount = pseudocount,
                    log_base = log_base)
                }, error = function(e) {
                  stop("C++ flexible paired bootstrap failed: ", e$message)
                })
            } else {
                stop("Failed to prepare paired bootstrap data from pair_ids")
            }
        } else {
            # Use C++ accelerated version for independent bootstrap (10-15x
            # faster)
            bootstrap_dist <- divergence_bootstrap_compute_cpp_wrapper(x = x, y = y,
                q = q, nboot = as.integer(nboot), paired = FALSE, pseudocount = pseudocount,
                log_base = log_base)
        }

        alpha <- (1 - ci)/2

        if (method == "percentile") {
            lower_ci <- stats::quantile(bootstrap_dist, probs = alpha, names = FALSE)
            upper_ci <- stats::quantile(bootstrap_dist, probs = 1 - alpha, names = FALSE)
        } else if (method == "bca") {
            # BCa not appropriate for divergence (requires two-sample
            # jackknife). Fall back to percentile method with warning.
            warning("BCa bootstrap not available for divergence (requires ",
                "two-sample jackknife acceleration). Falling back to percentile method.",
                call. = FALSE)
            lower_ci <- stats::quantile(bootstrap_dist, probs = alpha, names = FALSE)
            upper_ci <- stats::quantile(bootstrap_dist, probs = 1 - alpha, names = FALSE)
        }
    } else {
        lower_ci <- NA_real_
        upper_ci <- NA_real_
    }

    return(list(estimate = point_est, lower_ci = lower_ci, upper_ci = upper_ci, q = q,
        nboot = nboot,
        method = if (nboot > 0) {
            if (method == "bca") "percentile" else method
        } else NA_character_))
}


#' Compute Tsallis Divergence Between Two Count Vectors
#'
#' Computes scalar Tsallis divergence D_q(p || q) using the Furuichi formula.
#'

#' @noRd
.tsallis_divergence_scalar <- function(x, y, q_val, pseudocount = 0.5, log_base = exp(1)) {
    # Validate input vectors
    if (length(x) == 0 || length(y) == 0) {
        return(NA_real_)
    }

    # BUGFIX: Ensure x and y have equal length (required for divergence)
    if (length(x) != length(y)) {
        # This can happen if paired bootstrap resampling produces unequal group
        # sizes Return NA rather than crashing
        return(NA_real_)
    }

    # Normalize to probabilities
    p <- (x + pseudocount)/(sum(x) + length(x) * pseudocount)
    r <- (y + pseudocount)/(sum(y) + length(y) * pseudocount)

    if (any(is.na(p)) || any(is.na(r))) {
        return(NA_real_)
    }

    # q = 0 limit (convention 0^0 = 0) is evaluated on the UNCLAMPED,
    # PRE-PSEUDOCOUNT support so the limit keeps its support-difference
    # meaning: D_0(p||r) = 1 - sum_{i: P_i > 0} R_i (the R-mass on P's zero
    # support). With pseudocount > 0 the regularized vectors are all-positive
    # and would make D_0 identically 0; computing it on the raw counts avoids
    # that degeneracy (review_divergence.md, Option A).
    q_tol <- 1e-10
    if (q_val == 0) {
        if (sum(x, na.rm = TRUE) <= 0 || sum(y, na.rm = TRUE) <= 0) {
            return(NA_real_)
        }
        p0 <- x/sum(x)
        r0 <- y/sum(y)
        div <- 1 - sum(r0[p0 > 0], na.rm = TRUE)
        if (div < 0 && div > -1e-12) {
            div <- 0
        }
        return(div)
    }

    # Only apply min-probability clamping when
    # pseudocount is zero. When pseudocount > 0, it already handles zero
    # probabilities --- applying both is a double-correction that distorts
    # divergence values. When pseudocount == 0, min_prob acts as a safety
    # net against log(0) and power-of-zero numerical issues.
    if (pseudocount < 1e-10) {
        min_prob <- 1e-10
        p[p < min_prob] <- min_prob
        r[r < min_prob] <- min_prob

        # Re-normalize to maintain probability constraint (sum = 1)
        p <- p/sum(p)
        r <- r/sum(r)
    }

    # Compute Tsallis divergence (Furuichi 2006; van Erven & Harremoës 2014):
    #   D_q(p||r) = (sum_i p_i^q * r_i^(1-q) - 1) / (q - 1)   for q > 0, q != 1
    #   D_1(p||r) = sum_i p_i * log(p_i / r_i)               (KL limit, q = 1)
    #   D_0(p||r) = 1 - sum_{i: p_i > 0} r_i                 (q -> 0 limit with
    #               the convention 0^0 = 0, i.e. the Q-mass on P's zero
    #               support; under pseudocount > 0 every bin is positive and
    #               D_0 = 0).
    # Support violations (p_i > 0 with r_i = 0) make D_q = +Inf for q >= 1 in
    # the exact theory; the pseudocount/min-probability regularization keeps
    # the estimate finite here.
    # The formula is non-negative for valid distributions; abs() is NOT used,
    # only tiny negative roundoff is clamped to zero.
    q_tol <- 1e-10

    if (q_val > 0 && abs(q_val - 1) < q_tol) {
        # KL divergence (q = 1, exact up to numerical tolerance 1e-10):
        # lim_{q->1} D_q = sum(p * log(p/r))
        div <- sum(p * log(p/r), na.rm = TRUE)
    } else if (q_val > 0) {
        # Standard Tsallis divergence: D_q(p||r) = (sum(p^q * r^(1-q)) - 1)/(q - 1)
        p_power <- p^q_val
        r_power <- r^(1 - q_val)

        # Check for numerical issues (inf, nan, underflow)
        if (any(is.nan(p_power)) || any(is.infinite(p_power)) || any(is.nan(r_power)) ||
            any(is.infinite(r_power))) {
            # Log-space computation for numerical stability when q is far from
            # 1
            log_p_power <- q_val * log(p + 1e-10)
            log_r_power <- (1 - q_val) * log(r + 1e-10)
            sum_term <- sum(exp(log_p_power + log_r_power), na.rm = TRUE)
        } else {
            sum_term <- sum(p_power * r_power, na.rm = TRUE)
        }

        div <- (sum_term - 1)/(q_val - 1)
    } else {
        # Invalid q value
        return(NA_real_)
    }

    if (is.nan(div) || !is.finite(div)) {
        return(NA_real_)
    }

    # Numerical safety: the formula is non-negative for valid distributions.
    # Clamp only tiny negative roundoff to zero; larger negatives would
    # indicate a real anomaly and are left visible.
    if (div < 0 && div > -1e-12) {
        div <- 0
    }

    # Log_base normalization only applies to the
    # q\u21921 (KL divergence) limit. The Tsallis divergence for q\u22601 is scale-invariant
    # and does not involve a logarithm base.
    if (q_val > 0 && abs(q_val - 1) < q_tol && log_base != exp(1)) {
        div <- div/log(log_base)
    }

    return(div)
}

#' Vectorized Tsallis Divergence Calculation (OPTIMIZED)
#'
#' Computes Tsallis divergence for multiple q-values simultaneously using
#' vectorized operations.
#' This provides 2-3x speedup compared to sequential q-value loops by:
#' 1. Computing p^q and r^(1-q) matrices once
#' 2. Reusing these matrices for all q-values
#' 3. Vectorizing the sum and divergence formula computation
#'
#' @param x,y Count vectors (samples for each isoform)
#' @param q_vals Numeric vector of q-values to compute divergence for
#' @param pseudocount Pseudo-count for probability normalization (default: 0.5)
#' @param log_base Log base for divergence scaling (default: e)
#'
#' @return Numeric vector of divergence values, one per q-value
#'
#' @noRd
.tsallis_divergence_vector <- function(x, y, q_vals, pseudocount = 0.5, log_base = exp(1)) {
    # Input validation (identical to the previous pure-R implementation)
    if (length(x) == 0 || length(y) == 0) {
        return(rep(NA_real_, length(q_vals)))
    }

    # CRITICAL: Ensure x and y have equal length (required for divergence) This
    # can happen if paired bootstrap or group extraction produces unequal
    # lengths
    if (length(x) != length(y)) {
        return(rep(NA_real_, length(q_vals)))
    }

    if (any(is.na(x)) || any(is.na(y))) {
        return(rep(NA_real_, length(q_vals)))
    }

    # Multi-q C++ kernel: ONE normalization pass for all q, exact mirror of
    # the R semantics (pseudocount normalization, q=0 limit on unclamped
    # probabilities, min_prob clamp only for pseudocount==0, KL limit with
    # log_base correction, log-space fallback, roundoff clamp, non-finite
    # to NA).
    tsallis_divergence_vector_cpp(as.numeric(x), as.numeric(y), as.numeric(q_vals),
        pseudocount, log_base)
}

