# ============================================================================
# HELPER FUNCTIONS FOR BOOTSTRAP CI ANALYSIS (10 total)
# ============================================================================

#' Internal: Auto-select nboot based on input size and method

#' @noRd
.bootstrap_auto_select_nboot <- function(n_genes, use_bca, nthreads) {
    .suggest_nboot(n_genes, use_bca = use_bca, nthreads = nthreads)
}

#' Internal: Validate all bootstrap input parameters

#' @noRd
.bootstrap_validate_inputs <- function(x, q, nboot, ci, paired, show_messages = FALSE) {
    if (!is.numeric(x) || any(x < 0, na.rm = TRUE)) {
        stop("x must be a vector of non-negative numeric values.")
    }
    if (!is.numeric(q) || any(q < 0)) {
        stop("q must be non-negative numeric value(s) (q >= 0).")
    }

    # PRIORITY 2: Extreme q value validation (Issue #8) q=0 is valid (species
    # richness), but q values in (0, 0.001) or > 100 need caution
    if (any(q > 0 & q < 0.001)) {
        warning("Detected very small q values in (0, 0.001). Numerical behavior not well-tested. ",
            "Consider using q >= 0.001 or q = 0 for species richness.")
    }
    if (any(q > 100)) {
        warning("Detected extreme q values > 100 (very high weighting). Behavior not well-tested. ",
            "Consider using q <= 100.")
    }

    if (!is.numeric(nboot) || nboot < 1) {
        stop("nboot must be a numeric value >= 1.")
    }
    if (nboot < 10) {
        if (show_messages)
            warning("nboot = ", nboot, " is below minimum recommended (10). ", "CI bounds may be unreliable. Consider nboot >= 50 for production use.")
    } else if (nboot < 100 && !isTRUE(getOption("TSENAT.suppress_nboot_warning")) &&
        show_messages) {
        warning("nboot = ", nboot, " is below recommended minimum (100). ", "Consider nboot >= 100 for stable CI estimates (jackknife resampling literature, e.g., Zhang & Cao 2023).")
    }
    if (!is.numeric(ci) || ci <= 0 || ci >= 1) {
        stop("ci must be a probability in (0, 1).")
    }
    if (!is.logical(paired) || length(paired) != 1) {
        stop("'paired' must be a single logical value (TRUE or FALSE)")
    }
    if (isTRUE(paired) && (length(x)%%2 != 0)) {
        stop("For paired=TRUE, data must have even length")
    }

    total_count <- sum(x, na.rm = TRUE)
    if (total_count < 10 && show_messages) {
        warning("Total count (", total_count, ") below recommended minimum (10-20).\n",
            "Bootstrap estimates may be unreliable (jackknife resampling literature, e.g., Zhang & Cao 2023).")
    }
}

#' Internal: Enhanced validation for bootstrap data quality and edge cases

#' @noRd
.validate_bootstrap_data <- function(x, effective_length = NULL, pseudocount = 0) {
    # Check 1: Empty input
    if (length(x) == 0) {
        stop("Input x must be a non-empty vector. Received empty vector.")
    }

    # Check 2: All zeros with no pseudocount
    total_count <- sum(x, na.rm = TRUE)
    if (total_count == 0 && pseudocount == 0) {
        stop("All counts are zero and pseudocount = 0. Bootstrap entropy is undefined.\n",
            "Either: (1) provide non-zero counts, (2) set pseudocount > 0, or (3) check data quality.")
    }

    # Check 3: All zeros but pseudocount provided (warning only)
    if (total_count == 0 && pseudocount > 0) {
        warning("All counts are zero. Adding pseudocount = ", pseudocount, " for calculation. Results represent artificial distribution.")
    }

    # Check 4: Effective length validation
    if (!is.null(effective_length)) {
        if (length(effective_length) != length(x)) {
            stop("Length mismatch: effective_length (length = ", length(effective_length),
                ") must match x (length = ", length(x), ")")
        }

        # Check for negative or zero effective_length
        if (any(effective_length <= 0, na.rm = TRUE)) {
            n_bad <- sum(effective_length <= 0, na.rm = TRUE)
            warning("Found ", n_bad, " position(s) with effective_length <= 0. ",
                "These will be set to NA in normalization.")
        }
    }

    # Check 5: Single isoform (entropy = 0)
    if (length(x) == 1) {
        warning("Single isoform detected (n = 1). Entropy will be 0 with zero-width CI [0, 0].")
    }

    # Check 6: Very high proportion of zeros
    n_zeros <- sum(x == 0)
    frac_zeros <- n_zeros/length(x)
    if (frac_zeros >= 0.9) {
        warning("High proportion of zeros (", round(frac_zeros * 100, 1), "%). Bootstrap distribution may be concentrated in few categories.")
    }

    # Validation passed
    invisible(TRUE)
}

#' Internal: Process matrix input with parallelization

#' @noRd
.bootstrap_process_matrix <- function(x, q, norm, nboot, ci, method, log_base, pseudocount,
    what, gene_name, verbose, include_diagnostics, use_job, nthreads, paired,
    resample_by = "read", counts_matrix = NULL) {
    if (!is.numeric(nthreads) || nthreads < 1)
        stop("'nthreads' must be positive")
    nthreads <- as.integer(nthreads)

    gene_names <- rownames(x) %||% paste0("Gene_", seq_len(nrow(x)))

    # DIAGNOSTIC: Check input matrix
    if (nrow(x) == 0) {
        stop("Bootstrap matrix has 0 rows. This typically means all genes were filtered out.",
            call. = FALSE)
    }

    # Use .bplapply for cross-platform parallel support (Windows compatible)
    results_list <- .bplapply(seq_len(nrow(x)), function(i) {
        .calculate_tsallis_entropy_bootstrap(x = x[i, ], se = NULL, res = NULL,
            top_n = 1, q = q, norm = norm, nboot = nboot, ci = ci, method = method,
            log_base = log_base, pseudocount = pseudocount, what = what, gene_name = gene_names[i],
            verbose = FALSE, include_diagnostics = include_diagnostics, use_job = use_job,
            nthreads = 1, paired = paired, resample_by = resample_by,
            counts_matrix = counts_matrix)
    }, nthreads = nthreads)

    # DIAGNOSTIC: Check output list
    if (length(results_list) != nrow(x)) {
        stop("Bootstrap results list length (", length(results_list), ") does not match input matrix rows (",
            nrow(x), "). This indicates bootstrap computation failed for some genes.",
            call. = FALSE)
    }

    if (length(results_list) > 0) {
        names(results_list) <- gene_names
    } else {
        stop("Bootstrap results_list is empty. Check that the input matrix has rows.",
            call. = FALSE)
    }
    structure(results_list, class = c("tsenat_bootstrap_ci_list", "list"))
}

#' Internal: Extract gene counts from SummarizedExperiment

#' @noRd
.bootstrap_extract_gene <- function(se, target_gene) {
    gene_tx_idx <- which(rownames(se) == target_gene)

    if (length(gene_tx_idx) == 0) {
        rd <- SummarizedExperiment::rowData(se)
        if (!is.null(rd)) {
            if ("gene_name" %in% colnames(rd)) {
                gene_tx_idx <- which(rd$gene_name == target_gene)
            } else if ("gene_id" %in% colnames(rd)) {
                gene_tx_idx <- which(rd$gene_id == target_gene)
            }
        }
    }

    if (length(gene_tx_idx) == 0)
        stop("Gene '", target_gene, "' not found in 'se'")

    as.numeric(colSums(as.matrix(SummarizedExperiment::assay(se, "counts")[gene_tx_idx,
        , drop = FALSE])))
}

#' Internal: Process SummarizedExperiment and results data.frame

#' @noRd
.bootstrap_process_se <- function(se, res, top_n, q, norm, nboot, ci, method, log_base,
    pseudocount, what, gene_name, verbose, include_diagnostics, use_job, paired) {
    if (!methods::is(se, "SummarizedExperiment"))
        stop("'se' must be SummarizedExperiment")
    if (!is.data.frame(res))
        stop("'res' must be data.frame")

    if (!("gene_id" %in% colnames(res))) {
        stop("'res' data.frame must have a 'gene_id' column containing gene identifiers")
    }
    res_genes <- res$gene_id
    top_genes <- head(res_genes, top_n)

    # For paired data, use lower threshold (pairs have less total count)
    min_count_threshold <- if (isTRUE(paired))
        5 else 10
    valid_genes <- character()

    for (gene in head(res_genes, top_n * 2)) {
        if (length(valid_genes) >= top_n)
            break
        gene_counts <- tryCatch(.bootstrap_extract_gene(se, gene), error = function(e) NULL)
        if (!is.null(gene_counts) && sum(gene_counts, na.rm = TRUE) >= min_count_threshold) {
            valid_genes <- c(valid_genes, gene)
        }
    }

    if (length(valid_genes) == 0) {
        warning("No genes with sufficient counts for bootstrap.")
        return(NULL)
    }

    top_genes <- valid_genes[seq_len(min(top_n, length(valid_genes)))]

    if (length(top_genes) > 1) {
        results_list <- lapply(seq_along(top_genes), function(i) {
            .calculate_tsallis_entropy_bootstrap(x = NULL, se = se, res = data.frame(gene_id = top_genes[i],
                row.names = i), top_n = 1, q = q, norm = norm, nboot = nboot, ci = ci,
                method = method, log_base = log_base, pseudocount = pseudocount,
                what = what, gene_name = top_genes[i], verbose = FALSE, include_diagnostics = include_diagnostics,
                use_job = use_job, nthreads = 1, paired = paired)
        })
        names(results_list) <- top_genes
        return(structure(results_list, class = c("tsenat_bootstrap_ci_list", "list")))
    }

    target_gene <- top_genes[1]
    if (is.null(gene_name))
        gene_name <- target_gene
    gene_counts <- .bootstrap_extract_gene(se, target_gene)

    .calculate_tsallis_entropy_bootstrap(x = gene_counts, q = q, norm = norm, nboot = nboot,
        ci = ci, method = method, log_base = log_base, pseudocount = pseudocount,
        what = what, gene_name = gene_name, verbose = verbose, include_diagnostics = include_diagnostics,
        use_job = use_job, paired = paired)
}

#' Internal: Process multiple q values

#' @noRd
.bootstrap_process_multiple_q <- function(x, q, norm, nboot, ci, method, log_base,
    pseudocount, what, gene_name, verbose, include_diagnostics, use_job, paired,
    effective_length = NULL, min_valid_frac = 0.75, resample_by = "read",
    counts_matrix = NULL) {
    resample_by <- match.arg(resample_by)

    # FAST PATH (optimization audit): standard read-level bootstrap on a
    # vector. ONE C++ call resamples once per iteration and evaluates ALL q,
    # instead of Q separate bootstrap runs. QC regeneration and degenerate
    # cases fall back to the original per-q pipeline for identical behavior.
    if (!paired && identical(resample_by, "read") && is.null(counts_matrix) &&
        identical(method, "percentile")) {
        fast <- .bootstrap_process_multiple_q_fast(x, q, norm, nboot, ci, log_base,
            pseudocount, what, include_diagnostics, use_job, effective_length,
            min_valid_frac)
        if (!is.null(fast)) {
            names(fast) <- paste0("q=", q)
            return(structure(fast, class = c("tsenat_bootstrap_ci_list", "list")))
        }
    }

    results_list <- lapply(q, function(q_val) {
        .calculate_tsallis_entropy_bootstrap(x = x, se = NULL, res = NULL, top_n = 1,
            q = q_val, norm = norm, nboot = nboot, ci = ci, method = method, log_base = log_base,
            pseudocount = pseudocount, what = what, gene_name = NULL, verbose = FALSE,
            include_diagnostics = include_diagnostics, use_job = use_job, paired = paired,
            effective_length = effective_length, min_valid_frac = min_valid_frac,
            resample_by = resample_by, counts_matrix = counts_matrix)
    })
    names(results_list) <- paste0("q=", q)
    structure(results_list, class = c("tsenat_bootstrap_ci_list", "list"))
}

#' Fast multi-q bootstrap: one resample -> all q (percentile, unpaired, read)
#'
#' Mirrors the data transformation and result assembly of the per-q pipeline.
#' Returns NULL when the fast path is not applicable (caller falls back).
#'
#' @noRd
.bootstrap_process_multiple_q_fast <- function(x, q, norm, nboot, ci, log_base,
    pseudocount, what, include_diagnostics, use_job, effective_length, min_valid_frac) {
    q <- as.numeric(q)
    n_q <- length(q)

    # Same resampling-input transformation as .bootstrap_compute_ci:
    # p_hat must equal the point-estimate proportions exactly (see
    # .prepare_bootstrap_resample in R/bootstrap.R).
    prep <- .prepare_bootstrap_resample(x, effective_length, pseudocount)
    x_for_calc <- prep$x
    pseudocount_eff <- prep$pseudocount

    # ONE C++ call: nboot x n_q matrix (same resample for all q)
    dist_matrix <- if (what == "S") {
        bootstrap_compute_multi_q_cpp_wrapper(x = x_for_calc, q = q, normalize = norm,
            nboot = nboot, log_base = log_base, pseudocount = pseudocount_eff)
    } else {
        mat <- bootstrap_compute_multi_q_cpp_wrapper(x = x_for_calc, q = q, normalize = FALSE,
            nboot = nboot, log_base = log_base, pseudocount = pseudocount_eff)
        # Hill number conversion per q (same as the single-q path)
        for (j in seq_len(n_q)) {
            qv <- q[j]
            colj <- mat[, j]
            if (abs(qv - 1) < 1e-06) {
                mat[, j] <- exp(colj)
            } else {
                base <- 1 - (qv - 1) * colj
                n_neg <- sum(base < 0, na.rm = TRUE)
                if (n_neg > 0) {
                  warning("Hill number conversion produced ", n_neg, " negative base values for q=",
                    qv, ". Clamping to small positive value (1e-10).", call. = FALSE)
                  base <- pmax(base, 1e-10)
                }
                mat[, j] <- base^(1/(1 - qv))
            }
        }
        mat
    }

    results <- vector("list", n_q)
    for (j in seq_len(n_q)) {
        q_val <- q[j]
        dist_j <- dist_matrix[, j]
        n_valid <- sum(!is.na(dist_j))

        # Degenerate/QC fallback: identical to the legacy per-q pipeline
        # (includes regeneration attempts and all-Na warnings)
        if (n_valid == 0 || (n_valid/nboot) < min_valid_frac) {
            results[[j]] <- .calculate_tsallis_entropy_bootstrap(x = x, se = NULL,
                res = NULL, top_n = 1, q = q_val, norm = norm, nboot = nboot,
                ci = ci, method = "percentile", log_base = log_base, pseudocount = pseudocount,
                what = what, gene_name = NULL, verbose = FALSE, include_diagnostics = include_diagnostics,
                use_job = use_job, paired = FALSE, effective_length = effective_length,
                min_valid_frac = min_valid_frac, resample_by = "read", counts_matrix = NULL)
            next
        }

        # Point estimate from the estimator itself on the raw input (same
        # transformation as the resampling probabilities, by construction).
        point_est <- .calculate_tsallis_entropy(x, q = q_val, norm = norm,
            what = what, log_base = log_base, pseudocount = pseudocount,
            effective_length = effective_length)
        ci_result <- .ci_percentile(dist_j, ci = ci)
        diag_list <- .bootstrap_compute_diag(point_est, dist_j, use_job, paired = FALSE,
            x = x_for_calc, q = q_val, norm = norm, nboot = nboot, ci = ci,
            method = "percentile", log_base = log_base, pseudocount = pseudocount,
            what = what, accel_factor = NA_real_)
        results[[j]] <- .bootstrap_assemble_result(point_est, ci_result, dist_j,
            ci, "percentile", nboot, diag_list, include_diagnostics, use_job)
    }
    results
}



#' Internal: Bootstrap resampling with quality control enforcement
#'
#' Wraps .bootstrap_resample_optimized() and enforces min_valid_frac by
#' regenerating
#' invalid (NA/NaN) replicates until the quality threshold is met.
#'
#' @noRd
.bootstrap_resample_with_quality_control <- function(x, q, norm, nboot, log_base,
    pseudocount, what, paired = FALSE, effective_length = NULL, min_valid_frac = 0.75,
    resample_by = c("read", "replicate"), counts_matrix = NULL) {
    resample_by <- match.arg(resample_by)
    # Generate initial bootstrap replicates
    bootstrap_dist <- .bootstrap_resample_optimized(x, q = q, norm = norm, nboot = nboot,
        log_base = log_base, pseudocount = pseudocount, what = what, paired = paired,
        effective_length = effective_length, resample_by = resample_by,
        counts_matrix = counts_matrix)

    # Count invalid replicates and check against threshold
    n_invalid <- sum(is.na(bootstrap_dist) | is.nan(bootstrap_dist))
    n_valid <- nboot - n_invalid
    valid_frac <- n_valid/nboot

    # If all replicates are valid, return early
    if (n_invalid == 0) {
        return(bootstrap_dist)
    }

    # Regenerate invalid replicates until min_valid_frac is met Database
    # validation (Ramsay (2005), Springer Series in Statistics, 2005):
    # Bootstrap must operate on raw data with consistency checks
    max_attempts <- 10
    attempt <- 1
    regenerated_total <- 0

    while (valid_frac < min_valid_frac && attempt <= max_attempts) {
        # Find indices of invalid replicates
        invalid_idx <- which(is.na(bootstrap_dist) | is.nan(bootstrap_dist))
        n_to_regenerate <- length(invalid_idx)

        if (n_to_regenerate == 0)
            break  # All replicates are valid

        # Regenerate only the invalid replicates
        replacement_dist <- .bootstrap_resample_optimized(x, q = q, norm = norm,
            nboot = n_to_regenerate, log_base = log_base, pseudocount = pseudocount,
            what = what, paired = paired, effective_length = effective_length,
            resample_by = resample_by, counts_matrix = counts_matrix)

        # Replace invalid replicates with regenerated ones
        bootstrap_dist[invalid_idx] <- replacement_dist
        regenerated_total <- regenerated_total + n_to_regenerate

        # Recount and check
        n_invalid <- sum(is.na(bootstrap_dist) | is.nan(bootstrap_dist))
        n_valid <- nboot - n_invalid
        valid_frac <- n_valid/nboot

        if (valid_frac >= min_valid_frac) {
            # Quality threshold met
            if (regenerated_total > 0) {
                message(sprintf("[Bootstrap QC] Regenerated %d replicates across %d attempt(s). Final valid_frac: %.1f%%",
                  regenerated_total, attempt, valid_frac * 100))
            }
            return(bootstrap_dist)
        }

        attempt <- attempt + 1
    }

    # If we exit the loop without meeting threshold, check severity and handle appropriately
    if (valid_frac < min_valid_frac) {
        # Create informative warning with diagnostic information
        warning_msg <- sprintf(
            "Bootstrap regeneration could not achieve min_valid_frac=%.0f%% (got %.1f%% after %d attempts, %d replicates regenerated). ",
            min_valid_frac * 100, valid_frac * 100, max_attempts, regenerated_total)
        
        # Escalate response based on severity
        if (valid_frac < 0.5) {
            # Critical: Less than 50% valid - data quality is severely compromised
            stop(warning_msg, 
                 "CRITICAL: Less than 50% valid bootstrap replicates. ",
                 "This indicates severe data quality issues. Check for: ",
                 "(1) All-zero counts, (2) Extreme sparsity, (3) Invalid effective_length values. ",
                 "Consider filtering genes with adequate read depth.")
        } else if (valid_frac < 0.75) {
            # Warning: 50-75% valid - proceed with caution
            warning(warning_msg,
                   "CAUTION: CI may be unreliable with only ", 
                   sprintf("%.1f%%", valid_frac*100), " valid replicates. ",
                   "Consider reviewing data quality and results interpretation.")
        }
    }

    bootstrap_dist
}

#' Internal: Compute bootstrap CI

#' @noRd
.bootstrap_compute_ci <- function(x, q, norm, nboot, ci, method, log_base, pseudocount,
    what, paired = FALSE, effective_length = NULL, min_valid_frac = 0.75,
    resample_by = c("read", "replicate"), counts_matrix = NULL) {
    resample_by <- match.arg(resample_by)
    # BOOTSTRAP INVARIANT (auditx follow-up, 2026-08): the bootstrap must
    # resample from EXACTLY the point-estimate proportions
    # T(x, l, c) = (x/l + c) / sum(x/l + c). The point estimate is computed by
    # the estimator itself on the raw input; the resampling input is prepared
    # by .prepare_bootstrap_resample(), which embeds the pseudocount on the
    # effective-abundance scale and rescales to the original depth (the scale
    # factor cancels in the probabilities, so p_hat = T(x, l, c) exactly).
    # Previously the point estimate was computed on depth-rescaled values with
    # the pseudocount added afterwards, which does not factor with the scale
    # factor and therefore disagreed with the assay's stored estimate when
    # c > 0.
    point_est <- .calculate_tsallis_entropy(x, q = q, norm = norm, what = what,
        log_base = log_base, pseudocount = pseudocount, effective_length = effective_length)

    prep <- .prepare_bootstrap_resample(x, effective_length, pseudocount)
    x_for_calc <- prep$x
    pseudocount_eff <- prep$pseudocount

    # Replicate-level input (transcripts × samples) receives the same T:
    # divide by effective length and embed the pseudocount per value.
    if (!is.null(counts_matrix) && is.matrix(counts_matrix) && !is.null(effective_length) &&
        length(effective_length) == nrow(counts_matrix)) {
        counts_matrix <- sweep(counts_matrix, 1, effective_length, "/") + pseudocount
        pseudocount_eff <- 0
    }

    # Use quality-controlled bootstrap resampling on the prepared data
    # Pass effective_length=NULL since the transformation is already applied
    # Enforces min_valid_frac by regenerating invalid replicates
    bootstrap_dist <- .bootstrap_resample_with_quality_control(x_for_calc, q = q,
        norm = norm, nboot = nboot, log_base = log_base, pseudocount = pseudocount_eff,
        what = what, paired = paired, effective_length = NULL, min_valid_frac = min_valid_frac,
        resample_by = resample_by, counts_matrix = counts_matrix)

    # After quality control, check if CI is computable
    n_valid_values <- sum(!is.na(bootstrap_dist))

    # If all bootstrap replicates are invalid, return NAs for CI
    if (n_valid_values == 0) {
        warning("All bootstrap replicates produced NA/NaN. Returning NA for confidence intervals.")
        return(list(point_est = point_est, bootstrap_dist = bootstrap_dist, ci_result = list(lower = NA_real_,
            upper = NA_real_), accel_factor = NA_real_))
    }

    if (method == "percentile") {
        ci_result <- .ci_percentile(bootstrap_dist, ci = ci)
        accel_factor <- NA_real_
    } else {
        # Pass the prepared resampling input and its embedded pseudocount; the
        # BCa jackknife inside .ci_bca evaluates the same transformation T.
        ci_result <- .ci_bca(x_for_calc, bootstrap_dist, q = q, norm = norm, ci = ci, log_base = log_base,
            pseudocount = pseudocount_eff, what = what, point_est = point_est)
        accel_factor <- if (!is.null(ci_result$a))
            ci_result$a else NA_real_
    }

    list(point_est = point_est, bootstrap_dist = bootstrap_dist, ci_result = ci_result,
        accel_factor = accel_factor)
}

#' Internal: Compute bootstrap diagnostics and JOB

#' @noRd
.bootstrap_compute_diag <- function(point_est, bootstrap_dist, use_job, paired = FALSE,
    x, q, norm, nboot, ci, method, log_base, pseudocount, what, accel_factor = NA_real_) {
    diagnostics <- list(effective_sample_size = .compute_effective_n(bootstrap_dist),
        skewness = .compute_skewness(bootstrap_dist), bias = point_est - median(bootstrap_dist,
            na.rm = TRUE), acceleration_factor = accel_factor)

    job_stability <- NULL
    if (isTRUE(use_job) && isTRUE(paired)) {
        warning("JOB not supported with paired=TRUE.")
    } else if (isTRUE(use_job) && length(x) >= 3) {
        job_stability <- .compute_job(x, q = q, norm = norm, nboot = nboot, ci = ci,
            method = method, log_base = log_base, pseudocount = pseudocount, what = what,
            paired = paired)
    } else if (isTRUE(use_job) && length(x) < 3) {
        warning("JOB requires n >= 3. Skipping.")
    }

    list(diagnostics = diagnostics, job_stability = job_stability)
}

#' Internal: Assemble final bootstrap result

#' @noRd
.bootstrap_assemble_result <- function(point_est, ci_result, bootstrap_dist, ci,
    method, nboot, diag_list, include_diagnostics, use_job) {
    result <- list(estimate = as.numeric(point_est), lower_ci = ci_result$lower,
        upper_ci = ci_result$upper, ci_level = ci, method = method, nboot = nboot,
        bootstrap_dist = bootstrap_dist)

    if (include_diagnostics && !is.null(diag_list$diagnostics)) {
        result$diagnostics <- diag_list$diagnostics
    }

    if (isTRUE(use_job) && !is.null(diag_list$job_stability)) {
        result$job_stability <- diag_list$job_stability[c("ci_lower_stable", "ci_upper_stable",
            "ci_width_variation", "bound_variability", "n_outlier_bounds")]
    }

    structure(result, class = "tsenat_bootstrap_ci")
}

#' Internal: Print bootstrap results to console

#' @noRd
.bootstrap_print_results <- function(result, gene_name, ci, verbose) {
    if (!is.null(gene_name) && verbose) {
        message("\nBootstrap Confidence Intervals: ", gene_name)
        message("Point estimate: ", sprintf("%.6f", result$estimate))
        message(sprintf("%d%% CI: [%.6f, %.6f]", as.integer(ci * 100), result$lower_ci,
            result$upper_ci))
        message("CI width: ", sprintf("%.6f", result$upper_ci - result$lower_ci))
        message("Interpretation: We are ", sprintf("%.0f%%", ci * 100), " confident the true Tsallis entropy lies within this range.")
    }
}
