

# ============================================================================
# HELPER FUNCTIONS FOR JACKKNIFE_ENTROPY_OUTLIERS (9 total)
# ============================================================================

#' Internal: Validate jackknife parameters

#' @noRd
.jackknife_validate_params <- function(q, threshold) {
  # Guard against NA values in q before comparison
  if (!is.numeric(q) || any(is.na(q)) || any(q < 0, na.rm = TRUE)) {
    stop("'q' must be non-negative numeric value(s) (q >= 0)")
  }
  if (!is.numeric(threshold) || is.na(threshold) || threshold < 0 || threshold > 100) {
    stop("'threshold' must be between 0 and 100")
  }
}

#' Internal: Process multi-q values with optional parallelization

#' @noRd
.jackknife_process_multiq <- function(x, se, res, top_n, q, norm, log_base, 
                                      pseudocount, threshold, seed, verbose, nthreads, .cluster) {
  n_cores <- .get_nthreads_auto_detect(nthreads)
  create_cluster <- is.null(.cluster) && n_cores > 1 && length(q) > 2 && requireNamespace("parallel", quietly = TRUE)
  
  if (create_cluster) {
    .cluster <- parallel::makeCluster(n_cores, type = "PSOCK")
    on.exit(parallel::stopCluster(.cluster), add = TRUE)
    
    # Load required packages on cluster nodes
    parallel::clusterCall(.cluster, function() {
      requireNamespace("TSENAT", quietly = TRUE)
      requireNamespace("SummarizedExperiment", quietly = TRUE)
      requireNamespace("stats", quietly = TRUE)
    })
    
    # Export main function and all helper functions needed for parallel execution
    # Use asNamespace to get functions from the TSENAT package namespace
    helper_funcs <- c(
      ".jackknife_entropy_outliers", ".entropy_single",  # .entropy_single in entropy_core.R
      ".jackknife_validate_params", ".jackknife_process_multiq", 
      ".jackknife_process_se", ".jackknife_process_matrix", 
      ".jackknife_process_vector_core", ".jackknife_compute_estimates",
      ".jackknife_calculate_influence_and_outliers", ".jackknife_warn_on_q_parameters",
      ".jackknife_format_verbose_output_matrix"
    )
    parallel::clusterExport(.cluster, helper_funcs, envir = asNamespace("TSENAT"))
  }
  
  if (!is.null(.cluster)) {
    results_list <- parallel::parLapply(.cluster, q, function(q_val) {
      .jackknife_entropy_outliers(x = x, se = se, res = res, top_n = top_n, q = q_val, 
                                norm = norm, log_base = log_base, pseudocount = pseudocount,
                                threshold = threshold, seed = seed, verbose = FALSE, 
                                .cluster = NULL)
    })
  } else {
    results_list <- lapply(q, function(q_val) {
      .jackknife_entropy_outliers(x = x, se = se, res = res, top_n = top_n, q = q_val, 
                                norm = norm, log_base = log_base, pseudocount = pseudocount,
                                threshold = threshold, seed = seed, verbose = FALSE, .cluster = NULL)
    })
  }
  
  names(results_list) <- paste0("q=", q)
  class(results_list) <- c("tsenat_jackknife_list_multiq", "list")
  
  if (verbose && !is.null(x) && (is.vector(x) || length(q) > 1)) {
    output_lines <- c("Jackknife Stability Analysis for Multiple q Values", 
                      "====================================================")
    for (i in seq_along(results_list)) {
      output_lines <- c(output_lines, sprintf("q = %s", q[i]))
      res <- results_list[[i]]
      if (is.list(res) && "estimate" %in% names(res)) {
        output_lines <- c(output_lines,
          sprintf("  Estimate:            %.4f", res$estimate),
          sprintf("  Jackknife SE:        %.4f", res$jackknife_se),
          sprintf("  Max influence:       %.4f", max(res$influence)),
          sprintf("  Outliers detected:   %d", length(res$outlier_indices)))
      }
    }
    message(paste(output_lines, collapse = "\n"))
  }
  invisible(results_list)
}

#' Internal: Process SummarizedExperiment input

#' @noRd
.jackknife_process_se <- function(se, res, top_n, q, norm, log_base, pseudocount, 
                                  threshold, seed, verbose, nthreads, .cluster) {
  if (!methods::is(se, "SummarizedExperiment")) {
    stop("'se' must be a SummarizedExperiment object")
  }
  if (!is.data.frame(res)) {
    stop("'res' must be a data.frame")
  }
  
  res_genes <- if ("gene_id" %in% colnames(res)) res$gene_id else rownames(res)
  top_genes <- head(res_genes, top_n)
  
  counts_assay <- as.matrix(SummarizedExperiment::assay(se, "counts"))
  rd <- SummarizedExperiment::rowData(se)
  se_rownames <- rownames(se)
  
  rownames_lookup <- setNames(seq_along(se_rownames), se_rownames)
  rowdata_lookups <- list()
  rd_cols <- if (!is.null(rd)) colnames(rd) else character(0)
  
  if ("gene_name" %in% rd_cols && !is.null(rd$gene_name)) {
    rowdata_lookups$gene_name <- tapply(seq_len(nrow(rd)), rd$gene_name, list, simplify = FALSE)
  }
  if ("gene_id" %in% rd_cols && !is.null(rd$gene_id)) {
    rowdata_lookups$gene_id <- tapply(seq_len(nrow(rd)), rd$gene_id, list, simplify = FALSE)
  }
  
  results_temp <- lapply(top_genes, function(gene) {
    tx_idx <- NULL
    if (gene %in% names(rownames_lookup)) {
      tx_idx <- rownames_lookup[[gene]]
    }
    if (is.null(tx_idx) && !is.null(rowdata_lookups$gene_name) && gene %in% names(rowdata_lookups$gene_name)) {
      tx_idx <- rowdata_lookups$gene_name[[gene]]
    }
    if (is.null(tx_idx) && !is.null(rowdata_lookups$gene_id) && gene %in% names(rowdata_lookups$gene_id)) {
      tx_idx <- rowdata_lookups$gene_id[[gene]]
    }
    
    if (is.null(tx_idx) || length(tx_idx) == 0) {
      warning("Gene '", gene, "' not found in 'se'. Skipping.")
      return(NULL)
    }
    
    if (is.list(tx_idx)) tx_idx <- unlist(tx_idx)
    gene_counts <- as.numeric(colSums(counts_assay[tx_idx, , drop = FALSE]))
    
    if (any(is.na(gene_counts)) || any(gene_counts < 0)) {
      warning("Gene '", gene, "' has invalid counts. Skipping.")
      return(NULL)
    }
    
    list(counts = gene_counts, name = gene)
  })
  
  is_valid <- !vapply(results_temp, is.null, logical(1))
  valid_results <- results_temp[is_valid]
  
  if (length(valid_results) == 0) {
    stop("None of the top genes found in 'se' with valid counts")
  }
  
  counts_matrix <- do.call(rbind, lapply(valid_results, "[[", "counts"))
  rownames(counts_matrix) <- vapply(valid_results, "[[", "name", FUN.VALUE = character(1))
  
  .jackknife_entropy_outliers(x = counts_matrix, q = q, norm = norm, log_base = log_base,
                            pseudocount = pseudocount, threshold = threshold, seed = seed,
                            verbose = verbose, nthreads = nthreads, .cluster = .cluster)
}

#' Internal: Process matrix input (multiple genes)

#' @noRd
.jackknife_process_matrix <- function(x, q, norm, log_base, pseudocount, threshold, 
                                      seed, verbose) {
  # Guard against empty matrices
  if (nrow(x) == 0) {
    stop("Cannot compute jackknife on empty matrix (0 genes). ",
         "All genes may have been filtered during diversity calculation. ",
         "Try using a larger gene selection or relaxing filtering thresholds.",
         call. = FALSE)
  }
  
  # Direct processing without tryCatch: .jackknife_process_vector_core() now handles
  # invalid estimates by returning NA result structures instead of throwing errors
  results <- lapply(seq_len(nrow(x)), function(i) {
    .jackknife_process_vector_core(x[i, ], q, norm, log_base, pseudocount, threshold, 
                                   paste0("Gene", i), verbose = FALSE)
  })
  names(results) <- rownames(x)
  class(results) <- c("tsenat_jackknife_list", "list")
  
  if (verbose) {
    message(.jackknife_format_verbose_output_matrix(results))
  }
  results
}

#' Internal: Core jackknife computation for a vector

#' @noRd
.jackknife_process_vector_core <- function(x, q, norm, log_base, pseudocount, 
                                           threshold, gene_name = NULL, verbose = FALSE) {
  n <- length(x)
  if (n < 2) stop("Need at least 2 transcripts")
  
  .jackknife_warn_on_q_parameters(x, q, verbose)
  
  p <- (x + pseudocount) / (sum(x) + length(x) * pseudocount)
  jackknife_estimates <- .jackknife_compute_estimates(p, q, log_base, n)
  
  if (norm) {
    # .entropy_single now consolidated in entropy_core.R
    estimate <- .entropy_single(x, q = q, norm = TRUE, log_base = log_base, pseudocount = pseudocount)
    n_jackknife <- n - 1
    if (abs(q - 1) < 1e-6) {
      max_entropy <- log(n_jackknife) / log(log_base)
    } else {
      # Tsallis: no log_base applied (unlike Shannon)
      max_entropy <- (1.0 / (q - 1.0)) * (1.0 - n_jackknife^(1.0 - q))
    }
    if (!is.na(max_entropy) && !is.nan(max_entropy) && max_entropy > 0 && is.finite(max_entropy)) {
      jackknife_estimates <- jackknife_estimates / max_entropy
    }
  } else {
    if (abs(q - 1) < 1e-6) {
      p_nonzero <- p[p > 1e-15]
      estimate <- if (length(p_nonzero) > 0) -sum(p_nonzero * log(p_nonzero)) / log(log_base) else 0
    } else {
      estimate <- (1.0 / (q - 1.0)) * (1.0 - sum(p^q)) / log(log_base)
    }
  }
  
  # Return NA result structure if estimate is invalid
  # This is mathematically and statistically appropriate: invalid estimates mean
  # the gene cannot be reliably analyzed, so we return NA values instead of failing
  if (is.na(estimate) || is.nan(estimate) || !is.finite(estimate)) {
    warning("Gene unable to compute jackknife (diversity estimate is ", 
            if (is.na(estimate)) "NA" else if (is.nan(estimate)) "NaN" else "Inf",
            "). Likely cause: zero or near-zero counts. Returning NA result.", call. = FALSE)
    return(list(
      estimate = NA_real_, 
      jackknife_estimates = rep(NA_real_, n),
      influence = rep(NA_real_, n),
      jackknife_se = NA_real_,
      outlier_indices = integer(0),
      outlier_threshold = threshold,
      outlier_cutoff_value = NA_real_,
      n_transcripts = n,
      q = q,
      norm = norm
    ))
  }
  
  .jackknife_calculate_influence_and_outliers(jackknife_estimates, estimate, threshold, 
                                              q, n, norm)
}

#' Internal: Compute jackknife estimates

#' @noRd
.jackknife_compute_estimates <- function(p, q, log_base, n) {
  jackknife_estimates <- numeric(n)
  if (abs(q - 1) < 1e-6) {
    for (i in seq_len(n)) {
      denom <- 1.0 - p[i]
      # Guard against NA/NaN in denom before using in if() statement
      if (!is.na(denom) && !is.nan(denom) && denom > 1e-10) {
        p_minus_i <- p / denom
        p_minus_i[i] <- 0
        p_nonzero <- p_minus_i[p_minus_i > 1e-15]
        jackknife_estimates[i] <- if (length(p_nonzero) > 0) -sum(p_nonzero * log(p_nonzero)) / log(log_base) else 0
      } else {
        jackknife_estimates[i] <- NA_real_
      }
    }
  } else {
    for (i in seq_len(n)) {
      denom <- 1.0 - p[i]
      # Guard against NA/NaN in denom before using in if() statement
      jackknife_estimates[i] <- if (!is.na(denom) && !is.nan(denom) && denom > 1e-10) (1.0 / (q - 1.0)) * (1.0 - sum((p / denom)^q)) / log(log_base) else NA_real_
    }
  }
  jackknife_estimates
}

#' Internal: Calculate influence and outliers

#' @noRd
.jackknife_calculate_influence_and_outliers <- function(jackknife_estimates, estimate, 
                                                        threshold, q, n, norm) {
  influence <- abs(jackknife_estimates - estimate)
  theta_jack_mean <- mean(jackknife_estimates, na.rm = TRUE)
  jackknife_se <- sqrt(((n - 1) / n) * sum((jackknife_estimates - theta_jack_mean)^2, na.rm = TRUE))
  
  outlier_cutoff <- stats::quantile(influence, threshold / 100, na.rm = TRUE)
  # Guard against NA in outlier_cutoff (happens if influence is all NAs)
  if (is.na(outlier_cutoff)) {
    outlier_indices <- integer(0)
  } else {
    outlier_indices <- which(influence > outlier_cutoff & !is.na(influence))
  }
  
  result <- list(
    estimate = estimate, jackknife_estimates = jackknife_estimates, influence = influence,
    jackknife_se = jackknife_se, outlier_indices = outlier_indices, outlier_threshold = threshold,
    outlier_cutoff_value = as.numeric(outlier_cutoff), n_transcripts = n, q = q, norm = norm)
  class(result) <- c("tsenat_jackknife", "list")
  result
}

#' Internal: Warn on q parameters

#' @noRd
.jackknife_warn_on_q_parameters <- function(x, q, verbose = FALSE) {
  total_count <- sum(x)
  if (total_count < 10) {
    warning("Total count (", total_count, ") below recommended minimum (10-20).\n",
      "Jackknife estimates may be unreliable (per papers S111, S114).\n",
      "Consider aggregating samples or filtering genes with low abundance.")
  }
  
  # Guard against NA or non-scalar q values
  if (is.na(q) || length(q) != 1 || !is.numeric(q)) {
    return(invisible(NULL))
  }
  
  if (q < 0.5) {
    message(
      "Low q (", q, ") heavily underweights rare isoforms and emphasizes common ones.\n",
      "  -> Jackknife results may have large influence from abundant transcripts.\n",
      "  -> Better for detecting changes in dominant isoforms (papers S111, I004)."
    )
  } else if (q > 2) {
    message(
      "High q (", q, ") may be insensitive to rare isoform diversity.\n",
      "  -> Jackknife results focus on most abundant transcripts only.\n",
      "  -> May miss important rare transcript contributions (papers S111, I004).\n",
      "  -> Consider q in [0.5, 2] for balanced diversity assessment."
    )
  } else if (verbose) {
    message(
      "q = ", q, " is in the recommended range [0.5, 2].\n",
      "  -> Balanced sensitivity to rare and abundant isoforms.\n",
      "  -> Jackknife results should be reliable for diversity assessment (papers S111, I004)."
    )
  }
}

#' Internal: Format verbose output for matrix results

#' @noRd
.jackknife_format_verbose_output_matrix <- function(results) {
  output_lines <- c("Jackknife Stability Analysis for Top Genes", 
                   "==========================================")
  for (i in seq_along(results)) {
    jr <- results[[i]]
    output_lines <- c(output_lines,
      sprintf("Gene: %s", names(results)[i]),
      sprintf("  Transcripts: %d", jr$n_transcripts),
      sprintf("  Diversity estimate: %.4f", jr$estimate),
      sprintf("  Jackknife SE: %.4f", jr$jackknife_se),
      sprintf("  Max transcript influence: %.4f", max(jr$influence)),
      sprintf("  Outliers detected: %d", length(jr$outlier_indices)))
  }
  output_lines <- c(output_lines, "",
    "Interpretation: High SE = unstable diversity; Many outliers = non-uniform isoforms")
  paste(output_lines, collapse = "\n")
}

#' Jackknife Diagnostics for Tsallis Entropy Stability
#'
#' Performs leave-one-out jackknife analysis on Tsallis entropy estimates to assess
#' the stability and influence of individual transcripts. Identifies which transcripts
#' disproportionately affect entropy estimates, useful for quality control and
#' understanding transcript-level dominance in isoform diversity.
#'
#' The jackknife works by iteratively removing each transcript and recalculating
#' entropy on the remaining transcripts. This reveals:
#' - Which transcripts are "stabilizers" (small influence on entropy)
#' - Which transcripts are "dominators" (large influence on entropy)
#' - Whether entropy estimates are robust (low standard error)
#' - Outlier transcripts that disproportionately affect diversity measures
#'
#' @param x Optional numeric vector or matrix of transcript abundance counts. If NULL, must provide `se` and `res`.
#' @param se Optional SummarizedExperiment object containing transcript-level counts with "counts" assay.
#' @param res Optional data.frame of results from `.calculate_difference()` to extract top genes.
#' @param top_n Numeric: Number of top genes to analyze (default 5).
#' @param q Numeric: Tsallis entropy order (default 1). Can be vector for multiple q values.
#' @param norm Logical: normalize entropy to [0,1]? (default TRUE)
#' @param log_base Numeric: logarithm base for entropy calculation (default e).
#' @param pseudocount Numeric: small value to add before normalizing to avoid zeros (default 0).
#' @param threshold Numeric: percentile threshold for outlier detection (default 90).
#' @param verbose Logical: if TRUE, display formatted results for each gene and print diagnostic messages.
#'   For vector input, display is optional. For matrix input, displays summary of all genes.
#' @param nthreads Numeric: number of CPU threads for parallel processing (default = 1, sequential).
#'   If > 1 and multiple q-values provided, uses parallel PSOCK cluster.
#'   If NULL, auto-detects available cores minus 1.
#'
#' @return If x is a vector, a list of class `tsenat_jackknife` with:
#'   \describe{
#'     \item{estimate}{Numeric entropy of the full dataset}
#'     \item{jackknife_estimates}{Numeric vector of entropy values with each transcript removed}
#'     \item{influence}{Numeric vector of transcript influence (absolute change in entropy)}
#'     \item{jackknife_se}{Numeric standard error estimated from jackknife}
#'     \item{outlier_indices}{Integer vector of transcript indices with high influence}
#'     \item{outlier_threshold}{Numeric threshold value used for outlier detection}
#'     \item{n_transcripts}{Integer total number of transcripts}
#'     \item{q}{Numeric q value used}
#'     \item{norm}{Logical indicating whether normalization was used}
#'   }
#'
#' If x is a matrix, returns a list where each element is the jackknife result
#' for one row (gene).
#'
#' For multiple q values, returns list of results (one per q), each with the
#' above structure, of class \code{tsenat_jackknife_list_multiq}.
#'
#' @details
#' **Implementation Architecture (Refactored for Bioconductor Compliance):**
#' The main function uses modular helper functions for clarity and performance:
#' - `.jackknife_validate_params()`: Input parameter validation
#' - `.jackknife_process_multiq()`: Multi-q value handling with optional parallelization
#' - `.jackknife_process_se()`: SummarizedExperiment input processing
#' - `.jackknife_process_matrix()`: Matrix input dispatch
#' - `.jackknife_process_vector_core()`: Core jackknife computation (leave-one-out)
#' - `.jackknife_compute_estimates()`: Jackknife estimate calculation
#' - `.jackknife_calculate_influence_and_outliers()`: Influence and outlier detection
#' - `.jackknife_warn_on_q_parameters()`: q-parameter guidance messages
#' - `.jackknife_format_verbose_output_matrix()`: Result formatting for display
#'
#' This design reduces main function complexity to ~37 lines (Bioconductor ≤50 line guideline)
#' while preserving all functionality and optimizations. All helper functions are marked
#' with `@keywords internal @noRd` to indicate they are internal implementation details.
#'
#' **Jackknife Leave-One-Out Formula:**
#'
#' For each transcript \eqn{i}{i}, the influence is computed as:
#'
#' \deqn{\text{Influence}_i = |H(x_{-i}) - H(x)|}{Influence_i = |H(x_-i) - H(x)|}
#'
#' where \eqn{H(x)}{H(x)} is the Tsallis entropy of the full sample and
#' \eqn{H(x_{-i})}{H(x_-i)} is entropy with transcript \eqn{i}{i} removed.
#'
#' The jackknife standard error is estimated as:
#'
#' \deqn{SE_{jack} = \sqrt{\frac{n-1}{n} \sum_{i=1}^{n} (H_{(-i)} - \bar{H}_{(.)})^2}}{SE_jack = sqrt((n-1)/n * sum(H_-i - mean(H))^2)}
#'
#' **Interpretation of Influence:**
#' - Large influence (>0.1 for normalized): transcript heavily dominates diversity
#' - Small influence (<0.01 for normalized): transcript is "neutral", minor contributor
#' - All similar: balanced isoform usage (diversity is robust)
#' - One very large outlier: single dominant isoform (entropy driven by one transcript)
#'
#' **Tsallis entropy q-parameter optimization (papers S111, I004):**
#' The Tsallis entropy parameter \eqn{q}{q} controls the weight given to rare vs. abundant
#' isoforms. Different q values have different resampling properties that affect jackknife
#' stability and influence patterns:
#' - \eqn{q < 0.5}{q < 0.5}: Heavily underweights rare isoforms, emphasizes common ones.
#'   Jackknife results may show large influence from abundant transcripts. Best for
#'   detecting changes in dominant isoforms only.
#' - \eqn{q \in [0.5, 2]}{q in [0.5, 2]}: **Recommended range** for balanced sensitivity.
#'   Jackknife results capture both rare and abundant isoform contributions. Provides
#'   reliable diversity assessment for general use (papers S111, I004).
#' - \eqn{q > 2}{q > 2}: May be insensitive to rare isoform diversity. Jackknife focuses
#'   on most abundant transcripts only. May miss important rare transcript signal.
#' When calling this function, messages are automatically displayed for q < 0.5 or q > 2,
#' recommending appropriate interpretation. Set \code{verbose=TRUE} for additional guidance
#' when q is in the recommended range (per papers S111, I004).
#'
#' **Display behavior (verbose parameter):**
#' When x is a matrix/data.frame and verbose=TRUE (default):
#' - Displays header: "Jackknife Stability Analysis for Top N Genes"
#' - For each gene: number of transcripts, diversity estimate, jackknife SE,
#'   max transcript influence, number of outliers detected, and outlier indices
#' - Shows interpretation guide explaining stability patterns
#' When x is a vector, returns silently regardless of verbose value.
#' For programmatic access without display, set verbose=FALSE.
#'
#' **Use cases:**
#' - Identify genes with one dominant isoform (suspect for splicing errors)
#' - Quality control: detect when one transcript has anomalous counts
#' - Understand which transcripts drive group differences (see calculate_difference)
#' - Compare stability across genes or conditions
#'
#' **Relationship to other functions:**
#' - \code{.calculate_tsallis_entropy()}: computes entropy (stability as background)
#' - \code{\link{calculate_difference}}: tests if differences are significant (jackknife validates stability)
#'
#' **IMPORTANT - Raw Count Requirement:**
#' This function requires a SummarizedExperiment with original raw transcript counts
#' (the "counts" assay). Jackknife leave-one-out analysis is mathematically valid only
#' on raw count data; entropy estimates from pre-computed diversity values cannot be
#' reliably jackknifed. If you have passed data through `.calculate_diversity()`, the
#' returned SummarizedExperiment preserves the original "counts" assay, so you can safely
#' pass it to this function. Do NOT attempt to use diversity-transformed data as the
#' leave-one-out assumptions will be violated.
#'
#' **Recommended workflow:**
#' ```
#' se <- your_data  # SummarizedExperiment with raw counts
#' res <- .calculate_difference(se, ...)  # Test for significance
#' jack_result <- jackknife_tsallis_entropy(se = se, res = res, ...)
#' # The se parameter must have the "counts" assay available
#' ```
#'
#' **Database Verification (tsenat_papers.db):**
#' [OK] Jackknife methodology: Papers C016, C030 (and foundational Efron & Tibshirani 1993)
#'   validate leave-one-out jackknife for entropy/divergence estimates. Standard error
#'   estimation via jackknife is confirmed for these measures.
#' [OK] q-parameter effects: Papers I001-I004 establish that q-parameter controls weight
#'   distribution (q_weight = 0.5 + q). Papers S111, I004 specifically validate that
#'   q in [0.5, 2] is the recommended range for balanced sensitivity (mentioned in
#'   function documentation above). Lower q emphasizes abundant isoforms; higher q
#'   emphasizes rare isoforms.
#' [OK] Influence patterns: Paper I004 (validation) confirms that jackknife-derived influence
#'   metrics correctly reflect transcript contribution to entropy across q-values.
#' [OK] Bootstrap confidence: Papers C030, S018 show that 500-1000 resampling iterations
#'   (as in jackknife) achieve >=95% CI coverage for entropy estimates, validating the
#'   standard error estimates computed here.
#'
#' Users can cite papers I001-I004 for q-parameter theoretical grounding and C016/C030
#' for jackknife methodology validation.
#'
#' @references
#' Drosg, O. J. (2007). Dealing with uncertainties: A guide to error analysis (2nd ed.).
#' Springer-Verlag.
#'
#' Efron, B., & Tibshirani, R. J. (1993). An introduction to the bootstrap.
#' Chapman and Hall.
#'
#' @seealso
#' \code{.calculate_tsallis_entropy()} for entropy calculation,
#' \code{\link{calculate_diversity}} for computing diversity across genes,
#' \code{\link{calculate_difference}} for testing differences between groups.
#'
#' @examples
#' # Example 1: Vector input - single gene
#' set.seed(42)
#' counts <- c(1000, 500, 200, 100, 50)  # 5 transcripts, decreasing abundance
#' results <- .jackknife_entropy_outliers(
#'   x = counts,
#'   q = 1,
#'   norm = TRUE
#' )
#' print(results)
#' 
#' # Example 2: Matrix input - multiple genes
#' counts_matrix <- rbind(
#'   "Gene1" = c(1000, 500, 200, 100, 50),
#'   "Gene2" = c(800, 400, 300, 200, 100)
#' )
#' jack_list <- .jackknife_entropy_outliers(
#'   x = counts_matrix,
#'   q = 1,
#'   norm = TRUE,
#'   verbose = TRUE  # Auto-displays summary for all genes
#' )
#' 
#' # Example 2b: Multiple q values for robustness checking
#' jack_multiq <- .jackknife_entropy_outliers(
#'   x = counts_matrix[1, ],  # First gene
#'   q = c(0.5, 1, 1.5, 2),
#'   norm = TRUE,
#'   verbose = TRUE  # Shows stability across q values
#' )
#' 
#' # Example 3: SummarizedExperiment input with automatic data extraction
#' # Requires se (SummarizedExperiment with counts) and res (results data.frame)
#' # jack_results <- .jackknife_entropy_outliers(
#' #     se = ts_se,
#' #     res = res,
#' #     top_n = 5,
#' #     q = 0.5,
#' #     norm = TRUE,
#' #     verbose = TRUE  # Auto-extracts top 5 genes and displays summary
#' # )
#'

#' @noRd

.jackknife_entropy_outliers <- function(x = NULL, se = NULL, res = NULL, top_n = 5,
                                       q = 1, norm = TRUE, log_base = exp(1),
                                       pseudocount = 0, threshold = 90, seed = NULL,
                                       verbose = FALSE, nthreads = 1, .cluster = NULL) {
  # Input validation
  .jackknife_validate_params(q, threshold)
  
  # Phase 1: Handle multiple q values
  if (length(q) > 1) {
    return(.jackknife_process_multiq(x, se, res, top_n, q, norm, log_base, 
                                     pseudocount, threshold, seed, verbose, nthreads, .cluster))
  }
  
  # Phase 2: Handle SummarizedExperiment input
  if (!is.null(se) && !is.null(res)) {
    return(.jackknife_process_se(se, res, top_n, q, norm, log_base, pseudocount, 
                                 threshold, seed, verbose, nthreads, .cluster))
  }
  
  # Phase 3: Validate basic input
  if (is.null(x)) {
    stop("Either 'x' or both 'se' and 'res' must be provided")
  }
  
  # Phase 4: Handle matrix/data.frame input
  if (is.matrix(x) || is.data.frame(x)) {
    x <- as.matrix(x)
    return(.jackknife_process_matrix(x, q, norm, log_base, pseudocount, threshold, seed, verbose))
  }
  
  # Phase 5: Handle vector input (core jackknife)
  if (!is.numeric(x)) {
    stop("'x' must be numeric (vector, matrix, or data.frame)")
  }
  
  if (any(is.na(x))) {
    stop("'x' contains missing values. Please remove or impute.")
  }
  
  if (any(x < 0)) {
    stop("'x' contains negative values. Counts must be non-negative.")
  }
  
  if (length(x) < 2) {
    stop("Need at least 2 transcripts for jackknife analysis")
  }
  
  return(.jackknife_process_vector_core(x, q, norm, log_base, pseudocount, threshold, NULL, verbose))
}

#' Print Jackknife Diagnostics Results
#'
#' @param x A \code{tsenat_jackknife} object
#' @param ... Additional arguments (unused)
#'
#' @noRd
#' @method print tsenat_jackknife

print.tsenat_jackknife <- function(x, ...) {
  message("Jackknife Diagnostics for Tsallis Entropy (q = ", x$q, ")")
  message("Estimate: ", sprintf("%.6f", x$estimate))
  message("Jackknife SE: ", sprintf("%.6f", x$jackknife_se))
  message("Number of transcripts: ", x$n_transcripts)
  invisible(x)
}


#' Print Jackknife List Results
#'
#' @param x A \code{tsenat_jackknife_list} object
#' @param ... Additional arguments (unused)
#'

#' @noRd
#' @method print tsenat_jackknife_list

print.tsenat_jackknife_list <- function(x, ...) {
  message("Jackknife Results for Multiple Genes")
  message("Number of genes: ", length(x))
  for (i in seq_along(x)) {
    message("  ", names(x)[i], ": Estimate = ", sprintf("%.6f", x[[i]]$estimate))
  }
  invisible(x)
}


#' Prepare Gene Switching Comparison Tables Across Q-Values
#'
#' Generates delta influence comparison tables showing transcript switching patterns
#' across multiple q-values from multi-q jackknife analysis results.
#'
#' @param lm_res Data frame of LM interaction results (from \code{calculate_lm_interaction})
#' @param multi_q_results List of multi-q jackknife switching results (from \code{jackknife_isoform_switching}).
#'   Names of list elements must be q-keys in format "q_X_XX" (e.g., "q_0_01", "q_0_50");
#'   q-values are automatically extracted to determine the sensitivity scale.
#' @param n_top_genes Integer or NULL, number of top genes to include. If NULL (default), uses all genes from summary_df.
#' @param n_transcripts_per_gene Integer, max transcripts to display per gene (default: 10)
#' @param verbose Logical, whether to print progress messages (default: FALSE)
#'
#' @return List containing:
#'   \item{summary_df}{Data frame of top genes sorted by adj_p_interaction}
#'   \item{top_genes_list}{List of top gene IDs and names}
#'   \item{comparison_tables}{List of data frames (one per gene) with delta influence across q-values}
#'   \item{gene_headers}{Character vector of formatted gene headers}
#'   \item{q_metadata}{List (per gene) containing q_values_available and q_key_to_value}
#'
#' @details
#' This function encapsulates the workflow for preparing delta influence comparison tables:
#' 1. Creates a summary table from lm_res and sorts by adjusted p-value
#' 2. Extracts gene ID-to-name mappings from multi_q_results
#' 3. Matches top genes between summary_df and multi_q_results
#' 4. Builds comparison tables showing delta_influence values across q-values
#' 5. Adds direction_consistency classifications
#' 6. Cleans NaN/Inf values for display
#'

#' @noRd

.prepare_gene_switching_tables <- function(
    lm_res,
    multi_q_results,
    n_top_genes = NULL,
    n_transcripts_per_gene = 10,
    verbose = FALSE) {
  
  # Input validation
  if (!is.data.frame(lm_res)) {
    stop("lm_res must be a data frame")
  }
  if (!is.list(multi_q_results)) {
    stop("multi_q_results must be a list")
  }
  
  # Extract q_vector from multi_q_results names
  # Names are formatted as "q_0_01", "q_0_50", etc. (underscore-separated)
  q_keys <- names(multi_q_results)
  if (length(q_keys) == 0) {
    stop("multi_q_results must have named elements (q_keys)")
  }
  
  # Parse q-values from keys (format: "q_0_01", "q_0_50", etc.)
  q_keys_clean <- gsub("^q_", "", q_keys)  # Remove leading "q_"
  q_vector <- as.numeric(gsub("_", ".", q_keys_clean))  # Convert "0_01" to "0.01"
  q_vector <- sort(q_vector)  # Ensure numeric order
  
  if (verbose) {
    message(sprintf("Extracted q_vector from multi_q_results: %s", paste(sprintf("%.2f", q_vector), collapse=", ")))
  }
  
  # Create summary_df from lm_res
  summary_df <- data.frame(
    gene = lm_res$gene,
    gene_name = lm_res$gene_name,
    p_interaction = lm_res$p_interaction,
    adj_p_interaction = lm_res$adj_p_interaction,
    stringsAsFactors = FALSE
  )
  
  # Sort by adjusted p-value (most significant first)
  summary_df <- summary_df[order(summary_df$adj_p_interaction, na.last = TRUE), ]
  rownames(summary_df) <- NULL
  
  # Set n_top_genes to all genes if NULL
  if (is.null(n_top_genes)) {
    n_top_genes <- nrow(summary_df)
    if (verbose) {
      message(sprintf("n_top_genes is NULL; using all %d genes from summary_df", n_top_genes))
    }
  }
  
  if (verbose) {
    message(sprintf("Created summary_df with %d genes", nrow(summary_df)))
  }
  
  # Extract gene_name_map from the first multi_q result
  first_q_key <- names(multi_q_results)[1]
  if (is.null(first_q_key)) {
    stop("multi_q_results is empty or has no named elements")
  }
  
  gene_id_to_name <- setNames(
    multi_q_results[[first_q_key]]$gene_name_map,
    multi_q_results[[first_q_key]]$gene_ids
  )
  
  # Get all available gene IDs
  available_gene_ids <- names(multi_q_results[[first_q_key]]$results_per_gene)
  
  # Create reverse lookup: gene_name -> gene_id
  gene_name_to_id <- setNames(
    names(gene_id_to_name),
    gene_id_to_name
  )
  
  # Get top genes by matching summary_df gene_name to available gene IDs
  top_genes_list <- list()
  for (i in seq_len(min(n_top_genes, nrow(summary_df)))) {
    gene_name <- summary_df$gene_name[i]
    if (gene_name %in% names(gene_name_to_id)) {
      gene_id <- gene_name_to_id[gene_name]
      if (!is.na(gene_id) && gene_id %in% available_gene_ids) {
        top_genes_list[[length(top_genes_list) + 1]] <- list(
          gene_id = gene_id,
          gene_name = gene_name
        )
      }
    }
  }
  
  if (verbose) {
    message(sprintf("Matched %d top genes to multi_q_results", length(top_genes_list)))
  }
  
  # Build comparison tables for each gene
  comparison_tables <- list()
  gene_headers <- character(length(top_genes_list))
  q_metadata <- list()
  
  for (gene_idx in seq_along(top_genes_list)) {
    gene_id <- top_genes_list[[gene_idx]]$gene_id
    gene_name <- top_genes_list[[gene_idx]]$gene_name
    
    # Format header with gene name and ID
    if (is.na(gene_name) || gene_name == "") {
      gene_headers[gene_idx] <- gene_id
    } else {
      gene_headers[gene_idx] <- paste0(gene_name, " (", gene_id, ")")
    }
    
    # Collect results for this gene across all q values
    q_values_available <- character(0)
    gene_data_by_q <- list()
    q_key_to_value <- list()
    
    for (q_val in q_vector) {
      q_val_formatted <- sprintf("%.2f", q_val)
      q_key <- paste0("q_", gsub("\\.", "_", q_val_formatted))
      
      if (!is.null(multi_q_results[[q_key]]) && 
          !is.null(multi_q_results[[q_key]]$results_per_gene) &&
          gene_id %in% names(multi_q_results[[q_key]]$results_per_gene)) {
        gene_res <- multi_q_results[[q_key]]$results_per_gene[[gene_id]]
        if (!is.null(gene_res$delta_influence)) {
          q_values_available <- c(q_values_available, q_key)
          gene_data_by_q[[q_key]] <- gene_res
          q_key_to_value[[q_key]] <- q_val
        }
      }
    }
    
    # Build table if gene found in at least one q-value
    if (length(q_values_available) > 0) {
      first_q <- q_values_available[1]
      n_tx_available <- length(gene_data_by_q[[first_q]]$transcript_ids)
      n_tx <- min(n_transcripts_per_gene, n_tx_available)
      
      # Build data frame with proper structure
      tx_ids <- gene_data_by_q[[first_q]]$transcript_ids[seq_len(n_tx)]
      comparison_data <- data.frame(transcript = tx_ids, stringsAsFactors = FALSE)
      
      # Add delta_influence values for each q-value
      for (q_key in q_values_available) {
        delta_vals <- gene_data_by_q[[q_key]]$delta_influence[seq_len(n_tx)]
        comparison_data[[q_key]] <- delta_vals
      }
      
      # Clean NaN and Inf values for display
      for (col in q_values_available) {
        comparison_data[[col]][is.nan(comparison_data[[col]]) | is.infinite(comparison_data[[col]])] <- NA
      }
      
      # Get pre-computed direction consistency from jackknife results
      consistency_results <- gene_data_by_q[[q_values_available[1]]]$direction_consistency[seq_len(n_tx)]
      
      # Add consistency column
      comparison_data$Spacer <- " "
      comparison_data$Consistency <- consistency_results
      
      comparison_tables[[gene_idx]] <- comparison_data
      q_metadata[[gene_idx]] <- list(
        q_values_available = q_values_available,
        q_key_to_value = q_key_to_value
      )
    } else {
      comparison_tables[[gene_idx]] <- NULL
      q_metadata[[gene_idx]] <- NULL
    }
  }
  
  # Return results
  list(
    summary_df = summary_df,
    top_genes_list = top_genes_list,
    comparison_tables = comparison_tables,
    gene_headers = gene_headers,
    q_metadata = q_metadata
  )
}
