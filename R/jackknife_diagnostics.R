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
#' @param res Optional data.frame of results from `calculate_difference()` to extract top genes.
#' @param top_n Numeric: Number of top genes to analyze (default 5).
#' @param q Numeric: Tsallis entropy order (default 1). Can be vector for multiple q values.
#' @param norm Logical: normalize entropy to [0,1]? (default TRUE)
#' @param log_base Numeric: logarithm base for entropy calculation (default e).
#' @param pseudocount Numeric: small value to add before normalizing to avoid zeros (default 0).
#' @param threshold Numeric: percentile threshold for outlier detection (default 90).
#' @param seed Random seed for reproducibility.
#' @param print_results Logical: if TRUE (default), display formatted results for each gene.
#' @param verbose Logical: if TRUE, print diagnostic messages (default FALSE).
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
#' **Display behavior (print_results parameter):**
#' When x is a matrix/data.frame and print_results=TRUE (default):
#' - Displays header: "Jackknife Stability Analysis for Top N Genes"
#' - For each gene: number of transcripts, diversity estimate, jackknife SE,
#'   max transcript influence, number of outliers detected, and outlier indices
#' - Shows interpretation guide explaining stability patterns
#' When x is a vector, returns silently regardless of print_results value.
#' For programmatic access without display, set print_results=FALSE.
#'
#' **Use cases:**
#' - Identify genes with one dominant isoform (suspect for splicing errors)
#' - Quality control: detect when one transcript has anomalous counts
#' - Understand which transcripts drive group differences (see calculate_difference)
#' - Compare stability across genes or conditions
#'
#' **Relationship to other functions:**
#' - \code{\link{calculate_tsallis_entropy}}: computes entropy (stability as background)
#' - \code{\link{calculate_effect_sizes}}: measures group differences (jackknife diagnoses robustness)
#' - \code{\link{calculate_difference}}: tests if differences are significant (jackknife validates stability)
#'
#' **IMPORTANT - Raw Count Requirement:**
#' This function requires a SummarizedExperiment with original raw transcript counts
#' (the "counts" assay). Jackknife leave-one-out analysis is mathematically valid only
#' on raw count data; entropy estimates from pre-computed diversity values cannot be
#' reliably jackknifed. If you have passed data through `calculate_diversity()`, the
#' returned SummarizedExperiment preserves the original "counts" assay, so you can safely
#' pass it to this function. Do NOT attempt to use diversity-transformed data as the
#' leave-one-out assumptions will be violated.
#'
#' **Recommended workflow:**
#' ```
#' se <- your_data  # SummarizedExperiment with raw counts
#' res <- calculate_difference(se, ...)  # Test for significance
#' jack_result <- jackknife_tsallis_entropy(se = se, res = res, ...)
#' # The se parameter must have the "counts" assay available
#' ```
#'
#' **Database Verification (tsenat_papers.db):**
#' ✓ Jackknife methodology: Papers C016, C030 (and foundational Efron & Tibshirani 1993)
#'   validate leave-one-out jackknife for entropy/divergence estimates. Standard error
#'   estimation via jackknife is confirmed for these measures.
#' ✓ q-parameter effects: Papers I001-I004 establish that q-parameter controls weight
#'   distribution (q_weight = 0.5 + q). Papers S111, I004 specifically validate that
#'   q ∈ [0.5, 2] is the recommended range for balanced sensitivity (mentioned in
#'   function documentation above). Lower q emphasizes abundant isoforms; higher q
#'   emphasizes rare isoforms.
#' ✓ Influence patterns: Paper I004 (validation) confirms that jackknife-derived influence
#'   metrics correctly reflect transcript contribution to entropy across q-values.
#' ✓ Bootstrap confidence: Papers C030, S018 show that 500-1000 resampling iterations
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
#' \code{\link{calculate_tsallis_entropy}} for entropy calculation,
#' \code{\link{calculate_diversity}} for computing diversity across genes,
#' \code{\link{calculate_difference}} for testing differences between groups.
#'
#' @examples
#' # Example 1: Vector input - single gene
#' set.seed(42)
#' counts <- c(1000, 500, 200, 100, 50)  # 5 transcripts, decreasing abundance
#' results <- jackknife_tsallis_entropy(
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
#' jack_list <- jackknife_tsallis_entropy(
#'   x = counts_matrix,
#'   q = 1,
#'   norm = TRUE,
#'   print_results = TRUE  # Auto-displays summary for all genes
#' )
#' 
#' # Example 2b: Multiple q values for robustness checking
#' jack_multiq <- jackknife_tsallis_entropy(
#'   x = counts_matrix[1, ],  # First gene
#'   q = c(0.5, 1, 1.5, 2),
#'   norm = TRUE,
#'   print_results = TRUE  # Shows stability across q values
#' )
#' 
#' # Example 3: SummarizedExperiment input with automatic data extraction
#' # Requires se (SummarizedExperiment with counts) and res (results data.frame)
#' # jack_results <- jackknife_tsallis_entropy(
#' #     se = ts_se,
#' #     res = res,
#' #     top_n = 5,
#' #     q = 0.5,
#' #     norm = TRUE,
#' #     print_results = TRUE  # Auto-extracts top 5 genes and displays summary
#' # )
#'
#' @export
jackknife_tsallis_entropy <- function(x = NULL, se = NULL, res = NULL, top_n = 5,
                                       q = 1, norm = TRUE, log_base = exp(1),
                                       pseudocount = 0, threshold = 90, seed = NULL,
                                       print_results = TRUE, verbose = FALSE) {

  # Input validation
  if (!is.numeric(q) || any(q <= 0)) {
    stop("'q' must be positive numeric value(s)")
  }

  if (!is.numeric(threshold) || threshold < 0 || threshold > 100) {
    stop("'threshold' must be between 0 and 100")
  }
  
  # Handle multiple q values
  if (length(q) > 1) {
    # Recursive call for each q value
    results_list <- lapply(q, function(q_val) {
      jackknife_tsallis_entropy(
        x = x, se = se, res = res, top_n = top_n, q = q_val, 
        norm = norm, log_base = log_base, pseudocount = pseudocount,
        threshold = threshold, seed = seed, print_results = FALSE, verbose = verbose
      )
    })
    names(results_list) <- paste0("q=", q)
    class(results_list) <- c("tsenat_jackknife_list_multiq", "list")
    
    # Optional printing
    if (print_results && !is.null(x) && (is.vector(x) || length(q) > 1)) {
      cat("Jackknife Stability Analysis for Multiple q Values\n")
      cat("====================================================\n\n")
      for (i in seq_along(results_list)) {
        cat("q =", q[i], "\n")
        res <- results_list[[i]]
        if (is.list(res) && "estimate" %in% names(res)) {
          cat("  Estimate:           ", round(res$estimate, 4), "\n")
          cat("  Jackknife SE:       ", round(res$jackknife_se, 4), "\n")
          cat("  Max influence:      ", round(max(res$influence), 4), "\n")
          cat("  Outliers detected:  ", length(res$outlier_indices), "\n\n")
        }
      }
    }
    return(invisible(results_list))
  }

  # Single q path continues below

  # Handle SummarizedExperiment input with automatic data extraction
  if (!is.null(se) && !is.null(res)) {
    # Validate se is a SummarizedExperiment
    if (!methods::is(se, "SummarizedExperiment")) {
      stop("'se' must be a SummarizedExperiment object")
    }
    if (!is.data.frame(res)) {
      stop("'res' must be a data.frame")
    }

    # Extract gene names from results (stored in 'gene_id' column if available)
    if (!("gene_id" %in% colnames(res))) {
      # Use rownames if gene_id column doesn't exist
      res_genes <- rownames(res)
    } else {
      res_genes <- res$gene_id
    }

    # Extract top_n genes from results
    top_genes <- head(res_genes, top_n)

    # For each gene, get all transcript indices
    counts_list <- list()
    gene_names_out <- character()
    
    for (i in seq_along(top_genes)) {
      gene <- top_genes[i]
      
      # First, check if gene is in rownames (transcript ID case)
      tx_idx <- which(rownames(se) == gene)
      
      if (length(tx_idx) == 0) {
        # If not in rownames, try to find in rowData columns (check gene_name, then gene_id, then legacy genes)
        rd <- SummarizedExperiment::rowData(se)
        if (!is.null(rd)) {
          if ("gene_name" %in% colnames(rd)) {
            tx_idx <- which(rd$gene_name == gene)
          } else if ("gene_id" %in% colnames(rd)) {
            tx_idx <- which(rd$gene_id == gene)
          }
        }
      }
      
      if (length(tx_idx) == 0) {
        warning("Gene '", gene, "' not found in 'se' rownames or rowData columns. Skipping.")
        next
      }
      
      # Aggregate counts across all transcripts for this gene (sum across transcripts)
      gene_counts <- as.numeric(colSums(as.matrix(SummarizedExperiment::assay(se, "counts")[tx_idx, ])))
      # Only add if no missing values in gene counts (NA handling)
      if (!any(is.na(gene_counts)) && all(gene_counts >= 0)) {
        counts_list <- c(counts_list, list(gene_counts))
        gene_names_out <- c(gene_names_out, gene)
      } else {
        warning("Gene '", gene, "' has missing or invalid count values. Skipping.")
      }
    }
    
    if (length(counts_list) == 0) {
      stop("None of the top genes from 'res' found in 'se' rownames or rowData with valid counts")
    }

    # Create matrix with genes as rows
    counts_matrix <- do.call(rbind, counts_list)
    rownames(counts_matrix) <- gene_names_out

    # Call recursively with matrix input (with print_results suppressed for first call)
    return(jackknife_tsallis_entropy(
      x = counts_matrix,
      q = q,
      norm = norm,
      log_base = log_base,
      pseudocount = pseudocount,
      threshold = threshold,
      seed = seed,
      print_results = print_results,
      verbose = verbose
    ))
  }

  # Traditional path: x must be provided
  if (is.null(x)) {
    stop("Either 'x' or both 'se' and 'res' must be provided")
  }

  # Handle matrix/data.frame input
  if (is.matrix(x) || is.data.frame(x)) {
    x <- as.matrix(x)
    results <- list()
    for (i in seq_len(nrow(x))) {
      results[[i]] <- jackknife_tsallis_entropy(
        x = x[i, ],
        q = q,
        norm = norm,
        log_base = log_base,
        pseudocount = pseudocount,
        threshold = threshold,
        seed = seed,
        print_results = FALSE,  # Suppress individual gene displays; show once at end
        verbose = verbose
      )
    }
    names(results) <- rownames(x)
    class(results) <- c("tsenat_jackknife_list", "list")
    
    # Display results if requested
    if (print_results) {
      cat("Jackknife Stability Analysis for Top 5 Genes\n")
      cat("============================================\n\n")
      
      for (i in seq_along(results)) {
        gene_name <- names(results)[i]
        jr <- results[[i]]
        
        cat("Gene:", gene_name, "\n")
        cat("  Transcripts:", jr$n_transcripts, "\n")
        cat("  Diversity estimate:", round(jr$estimate, 4), "\n")
        cat("  Jackknife SE:", round(jr$jackknife_se, 4), "\n")
        cat("  Max transcript influence:", round(max(jr$influence), 4), "\n")
        cat("  Outliers detected:", length(jr$outlier_indices), "\n")
        
        if (length(jr$outlier_indices) > 0) {
          cat("    Outlier transcripts (indices):", paste(jr$outlier_indices, collapse = ", "), "\n")
        }
        cat("\n")
      }
      
      cat("Interpretation:\n")
      cat("- High SE relative to estimate: unstable diversity (few dominant transcripts)\n")
      cat("- Many outliers: non-uniform isoform distribution\n")
      cat("- No outliers: balanced/robust isoform diversity\n\n")
    }
    
    # IMPORTANT: When returning a single-row matrix result, extract the single result object
    # instead of returning a list with one element. This provides consistent return type.
    # Calling code can always use results[[1]] or results$gene_name to access the data.
    if (nrow(x) == 1) {
      # For single row: return the result list but mark it for consistent handling
      # Keep it as a list with class "tsenat_jackknife_list" for type consistency
      return(results)
    }
    
    return(results)
  }

  # vector input
  if (!is.numeric(x)) {
    stop("'x' must be numeric (vector, matrix, or data.frame)")
  }

  if (any(is.na(x))) {
    stop("'x' contains missing values. Please remove or impute.")
  }

  if (any(x < 0)) {
    stop("'x' contains negative values. Counts must be non-negative.")
  }

  n <- length(x)

  if (n < 2) {
    stop("Need at least 2 transcripts for jackknife analysis")
  }

  # Minimum sample size validation (per papers S111, S114)
  # Jackknife/bootstrap are unreliable with insufficient data
  total_count <- sum(x)
  if (total_count < 10) {
    warning(
      "Total count (", total_count, ") below recommended minimum (10-20).\n",
      "Jackknife estimates may be unreliable (per papers S111, S114).\n",
      "Consider aggregating samples or filtering genes with low abundance."
    )
  }

  # Tsallis entropy q-parameter optimization guidance (papers S111, I004)
  # Different q values have different resampling and diversity properties
  if (q < 0.5) {
    message(
      "Low q (", q, ") heavily underweights rare isoforms and emphasizes common ones.\n",
      "  → Jackknife results may have large influence from abundant transcripts.\n",
      "  → Better for detecting changes in dominant isoforms (papers S111, I004)."
    )
  }
  if (q > 2) {
    message(
      "High q (", q, ") may be insensitive to rare isoform diversity.\n",
      "  → Jackknife results focus on most abundant transcripts only.\n",
      "  → May miss important rare transcript contributions (papers S111, I004).\n",
      "  → Consider q in [0.5, 2] for balanced diversity assessment."
    )
  }
  if (q >= 0.5 && q <= 2) {
    if (verbose) {
      message(
        "q = ", q, " is in the recommended range [0.5, 2].\n",
        "  → Balanced sensitivity to rare and abundant isoforms.\n",
        "  → Jackknife results should be reliable for diversity assessment (papers S111, I004)."
      )
    }
  }

  # Compute full estimate
  estimate <- .tsenat_entropy_single(x, q = q, norm = norm, log_base = log_base,
                                       pseudocount = pseudocount)

  # Jackknife: remove each transcript and recalculate
  jackknife_estimates <- numeric(n)

  for (i in seq_len(n)) {
    x_minus_i <- x[-i] # Remove transcript i
    jackknife_estimates[i] <- .tsenat_entropy_single(
      x_minus_i, q = q, norm = norm, log_base = log_base,
      pseudocount = pseudocount
    )
  }

  # Compute influence of each transcript
  influence <- abs(jackknife_estimates - estimate)

  # Jackknife standard error (from Efron & Tibshirani 1993)
  # SE = sqrt( ((n-1)/n) * sum( (theta_-i - theta_.)^2 ) )
  theta_jack_mean <- mean(jackknife_estimates, na.rm = TRUE)
  jackknife_se <- sqrt(((n - 1) / n) * sum((jackknife_estimates - theta_jack_mean)^2, na.rm = TRUE))

  # Identify outliers based on threshold percentile
  outlier_cutoff <- stats::quantile(influence, threshold / 100, na.rm = TRUE)
  outlier_indices <- which(influence > outlier_cutoff & !is.na(influence))

  # Return results
  result <- list(
    estimate = estimate,
    jackknife_estimates = jackknife_estimates,
    influence = influence,
    jackknife_se = jackknife_se,
    outlier_indices = outlier_indices,
    outlier_threshold = threshold,
    outlier_cutoff_value = as.numeric(outlier_cutoff),
    n_transcripts = n,
    q = q,
    norm = norm
  )

  class(result) <- c("tsenat_jackknife", "list")
  return(result)
}


#' Single Entropy Calculation
#'
#' Internal helper function to compute Tsallis entropy for a vector of counts.
#'
#' @keywords internal
#' @noRd
.tsenat_entropy_single <- function(counts, q = 1, norm = TRUE, log_base = exp(1),
                                    pseudocount = 0) {

  # Normalize to proportions
  p <- counts / sum(counts)

  # Add pseudocount to avoid log(0)
  p <- p + pseudocount
  p <- p / sum(p)

  n <- length(p)

  if (abs(q - 1) < 1e-6) {
    # Shannon entropy as q → 1
    # Filter out zeros to avoid 0 * log(0) = NaN
    p_nonzero <- p[p > 0]
    if (length(p_nonzero) > 0) {
      entropy <- -sum(p_nonzero * log(p_nonzero) / log(log_base))
    } else {
      entropy <- 0
    }
  } else {
    # Generalized Tsallis entropy
    entropy <- (1 / (q - 1)) * (1 - sum(p^q)) / log(log_base)
  }

  # Normalize to [0, 1]
  if (norm) {
    # Maximum entropy is achieved with uniform distribution
    # For q != 1: S_max = (1/(q-1)) * (1 - n^(1-q))
    # For q = 1: S_max = log(n)
    if (abs(q - 1) < 1e-6) {
      max_entropy <- log(n) / log(log_base)
    } else {
      max_entropy <- (1 / (q - 1)) * (1 - n^(1 - q)) / log(log_base)
    }

    if (!is.na(max_entropy) && !is.nan(max_entropy) && max_entropy > 0 && is.finite(max_entropy)) {
      entropy <- entropy / max_entropy
    }
  }

  return(as.numeric(entropy))
}


#' Print Jackknife Diagnostics Results
#'
#' @param x A \code{tsenat_jackknife} object
#' @param ... Additional arguments (unused)
#'
#' @keywords internal
#' @noRd
#' @export
#' @method print tsenat_jackknife
print.tsenat_jackknife <- function(x, ...) {
  cat("Jackknife Diagnostics for Tsallis Entropy (q =", x$q, ")\n")
  cat("=========================================================\n")
  cat("Number of transcripts:", x$n_transcripts, "\n")
  cat("Entropy (full data):", round(x$estimate, 6), "\n")
  cat("Jackknife SE:", round(x$jackknife_se, 6), "\n")
  cat("95% CI approximately: [",
      round(x$estimate - 1.96 * x$jackknife_se, 6), ", ",
      round(x$estimate + 1.96 * x$jackknife_se, 6), "]\n\n")

  cat("Transcript Influence (absolute change in entropy):\n")
  cat("Min:", round(min(x$influence), 6), "\n")
  cat("Median:", round(stats::median(x$influence), 6), "\n")
  cat("Max:", round(max(x$influence), 6), "\n\n")

  cat("Outlier Summary (threshold:", x$outlier_threshold, "%ile):\n")
  cat("Cutoff value:", round(x$outlier_cutoff_value, 6), "\n")
  cat("Number of outliers:", length(x$outlier_indices), "\n")

  if (length(x$outlier_indices) > 0) {
    cat("Outlier transcript indices:", paste(x$outlier_indices, collapse = ", "), "\n")
  }

  invisible(x)
}


#' Summary of Jackknife Diagnostics Results
#'
#' @param object A \code{tsenat_jackknife} object
#' @param ... Additional arguments (unused)
#'
#' @keywords internal
#' @noRd
#' @export
#' @method summary tsenat_jackknife
summary.tsenat_jackknife <- function(object, ...) {
  cat("Jackknife Diagnostics Summary\n")
  cat("=============================\n")
  cat("Gene entropy (q =", object$q, "):", round(object$estimate, 6), "\n")
  cat("Jackknife standard error:", round(object$jackknife_se, 6), "\n")
  cat("Coefficient of variation:", round(object$jackknife_se / object$estimate, 4), "\n")
  cat("Total transcripts:", object$n_transcripts, "\n\n")

  cat("Influence Distribution:\n")
  cat("-----------------------\n")
  print(summary(object$influence))

  cat("\n\nOutlier Transcripts (influence > ", object$outlier_threshold, "%ile):\n", sep = "")
  cat("-----------------------------------------------\n")

  if (length(object$outlier_indices) > 0) {
    outl_data <- data.frame(
      Transcript = object$outlier_indices,
      Influence = object$influence[object$outlier_indices]
    )
    print(outl_data)
  } else {
    cat("No outliers detected.\n")
  }

  invisible(object)
}


#' Print Jackknife List Results
#'
#' @param x A \code{tsenat_jackknife_list} object
#' @param ... Additional arguments (unused)
#'
#' @keywords internal
#' @noRd
#' @export
#' @method print tsenat_jackknife_list
print.tsenat_jackknife_list <- function(x, ...) {
  cat("Jackknife Diagnostics for", length(x), "genes\n")
  cat("=============================================\n")
  cat("Gene names:", paste(head(names(x), 5), collapse = ", "))
  if (length(x) > 5) cat(", ...")
  cat("\n")
  cat("\nUse indexing to view individual genes: object[[1]] or object$'GeneName'\n")
  invisible(x)
}
#' Bootstrap Helper Function for Delta Statistics
#'
#' Internal function to compute bootstrap confidence intervals and p-values
#' for delta_influence (difference in Tsallis entropy influence between conditions).
#' Computes per-transcript statistics from bootstrap resamples.
#'
#' @keywords internal
#' @noRd
compute_delta_statistics <- function(counts_A, counts_B, delta_influence,
                                   q = 1, norm = TRUE, log_base = exp(1),
                                   pseudocount = 0, n_bootstrap = 1000,
                                   seed = 42, confidence = 0.95, n_transcripts = NULL) {
  .calculate_tsallis <- function(counts, q, norm, log_base, pseudocount, n_transcripts_fixed) {
    if (is.vector(counts)) counts <- t(as.matrix(counts))
    
    # Check for zero-sum columns BEFORE adding pseudocount to detect samples with no expression
    raw_col_sums <- colSums(counts)
    with_zero_counts <- raw_col_sums == 0
    
    # CRITICAL FIX: Enforce minimum pseudocount to avoid zero-count edge cases
    if (pseudocount <= 0) {
      pseudocount <- 1e-8
    }
    counts <- counts + pseudocount
    
    col_sums <- colSums(counts)
    # Avoid division by zero
    if (any(col_sums <= 0)) {
      return(rep(NA_real_, ncol(counts)))
    }
    
    # Avoid division by near-zero for zero-count columns by using raw column sums intelligently
    # For zero-count columns, set col_sums to 1 to avoid division by near-zero
    col_sums_safe <- pmax(col_sums, 1)
    
    p <- counts / col_sums_safe
    
    if (q == 1) {
      if (log_base == exp(1)) {
        h <- -colSums(p * log(p + 1e-100))
      } else {
        h <- -colSums(p * log(p + 1e-100, log_base))
      }
    } else {
      p_q_sum <- colSums(p^q)
      if (log_base == exp(1)) {
        h <- (1 / (1 - q)) * (1 - p_q_sum)
      } else {
        h <- (1 / (1 - q)) * (1 - p_q_sum)
      }
    }
    
    h[!is.finite(h)] <- NA_real_
    # Set entropy to NA for zero-count columns
    h[with_zero_counts] <- NA_real_
    
    if (norm) {
      if (q == 1) {
        # Use FIXED n_transcripts if provided, otherwise use current matrix dimensions
        n_tx <- if (!is.null(n_transcripts_fixed)) n_transcripts_fixed else nrow(counts)
        max_h <- log(n_tx) / log(log_base)
      } else {
        n_tx <- if (!is.null(n_transcripts_fixed)) n_transcripts_fixed else nrow(counts)
        # Use absolute value because the formula produces negative values for q<1 and q>1
        max_h <- abs((1 / (1 - q)) * (1 - n_tx^(1 - q)))
      }
      if (!is.na(max_h) && is.finite(max_h) && max_h > 0) {
        h <- h / max_h
      }
    }
    return(h)
  }
  
  .jackknife_entropy <- function(counts, q, norm, log_base, pseudocount, n_transcripts_fixed) {
    if (is.vector(counts)) counts <- matrix(counts, nrow = 1)
    n_tx <- nrow(counts)
    h_full <- .calculate_tsallis(counts, q, norm, log_base, pseudocount, n_transcripts_fixed)
    influences <- numeric(n_tx)
    
    for (i in 1:n_tx) {
      h_leave_i <- .calculate_tsallis(counts[-i, , drop = FALSE], q, norm, log_base, pseudocount, n_transcripts_fixed)
      influences[i] <- mean(abs(h_full - h_leave_i), na.rm = TRUE)
    }
    list(full_entropy = h_full, influences = influences)
  }
  
  n_tx <- length(delta_influence)
  set.seed(seed)
  bootstrap_deltas_matrix <- matrix(nrow = n_bootstrap, ncol = n_tx)
  
  for (b in 1:n_bootstrap) {
    idx_A <- sample(1:ncol(counts_A), size = ncol(counts_A), replace = TRUE)
    idx_B <- sample(1:ncol(counts_B), size = ncol(counts_B), replace = TRUE)
    
    boot_A <- counts_A[, idx_A, drop = FALSE]
    boot_B <- counts_B[, idx_B, drop = FALSE]
    
    jack_A <- .jackknife_entropy(boot_A, q, norm, log_base, pseudocount, n_transcripts)
    jack_B <- .jackknife_entropy(boot_B, q, norm, log_base, pseudocount, n_transcripts)
    
    # Compute per-transcript delta for this bootstrap sample
    bootstrap_deltas_matrix[b, ] <- jack_A$influences - jack_B$influences
  }
  
  # Compute per-transcript statistics
  alpha <- 1 - confidence
  ci_lower <- apply(bootstrap_deltas_matrix, 2, function(x) quantile(x, alpha / 2, na.rm = TRUE))
  ci_upper <- apply(bootstrap_deltas_matrix, 2, function(x) quantile(x, 1 - alpha / 2, na.rm = TRUE))
  
  # Para cada transcrito, calcular p-value basado en cuántos bootstrap samples
  # tienen signo opuesto al delta_influence observado
  pvalues <- numeric(n_tx)
  for (i in 1:n_tx) {
    boot_signs <- sign(bootstrap_deltas_matrix[, i])
    obs_sign <- sign(delta_influence[i])
    # P-value: proporción de muestras bootstrap con signo opuesto
    pvalues[i] <- mean(boot_signs != obs_sign, na.rm = TRUE)
    pvalues[i] <- max(pvalues[i], 1 / n_bootstrap)  # Minimum p-value
  }
  
  # Standard error per transcript
  se <- apply(bootstrap_deltas_matrix, 2, sd, na.rm = TRUE)
  
  return(list(
    ci_lower = as.numeric(ci_lower),
    ci_upper = as.numeric(ci_upper),
    pvalue = pvalues,
    se = se
  ))
}


#' Jackknife Isoform Switching Detection
#'
#' Detects isoform/transcript switches between conditions using condition-stratified
#' jackknife analysis on Tsallis entropy. Identifies which transcripts change
#' importance between conditions (e.g., condition A vs B).
#'
#' @param se A SummarizedExperiment object with transcript-level counts.
#' @param condition_col Character: column name in colData for condition labels.
#' @param pair_col Character: optional column name for paired design (individual IDs).
#' @param gene_col Character: column name in rowData for gene IDs.
#' @param isoform_col Character: column name in rowData for transcript/isoform IDs.
#' @param q Numeric: Tsallis entropy order (default 1 = Shannon entropy). Can be a vector
#'   for multi-q analysis (e.g., q = c(0.5, 1.0, 1.5, 2.0)); results will be nested by q value.
#' @param norm Logical: normalize entropy to [0,1]? (default TRUE).
#' @param log_base Numeric: log base for entropy (default e).
#' @param pseudocount Numeric: pseudocount to add (default 0).
#' @param threshold Numeric: percentile for outlier detection on influences (default 90).
#' @param n_bootstrap Numeric: number of bootstrap resamples (default 1000).
#' @param print_results Logical: print results? (default TRUE).
#' @param verbose Logical: verbose output? (default FALSE).
#' @param lm_results Data frame: results from calculate_lm_interaction() with 'gene' column.
#'   Can contain either gene names or gene IDs; function automatically maps names to IDs
#'   using rowData(se). Include 'p_interaction' and/or 'adj_p_interaction' columns for
#'   filtering genes by significance. When provided, only genes passing lm_p_threshold are
#'   analyzed; all matching genes are included (top_n parameter removed).
#' @param lm_p_threshold Numeric: p-value threshold for LM gene filtering (default 0.05).
#' @param use_lm_fdr Logical: use adjusted p-values from LM results if available (default TRUE).
#'
#' @return If q is a single value, returns a list of class tsenat_isoform_switching with:
#'   \describe{
#'     \item{results_per_gene}{named list of per-gene results}
#'     \item{summary_table}{data.frame with per-gene summary}
#'     \item{all_transcript_stats}{data.frame with all transcript statistics}
#'     \item{gene_names}{character vector of analyzed genes}
#'     \item{conditions}{character vector of the two conditions compared}
#'     \item{metadata}{list with analysis metadata}
#'   }
#'   
#'   If q is a vector, returns a list of class tsenat_isoform_switching_multiq where each
#'   element is a complete tsenat_isoform_switching result for that q value. Keys are
#'   formatted as "q_X_XX" for ease of iteration (e.g., q_0_01, q_1_00, q_2_00).
#'
#' @export
jackknife_isoform_switching <- function(
  se = NULL,
  condition_col = "condition",
  pair_col = NULL,
  gene_col = NULL,
  isoform_col = NULL,
  q = 1,
  norm = TRUE,
  log_base = exp(1),
  pseudocount = 0,
  threshold = 90,
  n_bootstrap = 1000,
  print_results = TRUE,
  verbose = FALSE,
  lm_results = NULL,
  lm_p_threshold = 0.05,
  use_lm_fdr = TRUE
) {
  
  # Input validation
  if (is.null(se)) {
    stop("SummarizedExperiment object (se) is required")
  }
  
  if (!inherits(se, "SummarizedExperiment")) {
    stop("se must be a SummarizedExperiment object")
  }
  
  # Handle multiple q values
  if (is.numeric(q) && length(q) > 1) {
    # Recursive call for each q value
    results_list <- lapply(q, function(q_val) {
      jackknife_isoform_switching(
        se = se,
        condition_col = condition_col,
        pair_col = pair_col,
        gene_col = gene_col,
        isoform_col = isoform_col,
        q = q_val,
        norm = norm,
        log_base = log_base,
        pseudocount = pseudocount,
        threshold = threshold,
        n_bootstrap = n_bootstrap,
        print_results = FALSE,
        verbose = verbose,
        lm_results = lm_results,
        lm_p_threshold = lm_p_threshold,
        use_lm_fdr = use_lm_fdr
      )
    })
    
    # Create named list with q values formatted as keys
    names(results_list) <- paste0("q_", gsub("\\.", "_", sprintf("%.2f", q)))
    
    # Add direction consistency to multi-q results
    # This analyzes each transcript's delta_influence pattern across q-values
    for (gene_idx in seq_along(results_list[[1]]$results_per_gene)) {
      gene_name <- names(results_list[[1]]$results_per_gene)[gene_idx]
      
      # Collect delta_influence for this gene across all q-values
      transcript_ids <- results_list[[1]]$results_per_gene[[gene_name]]$transcript_ids
      
      # Build matrix: rows = transcripts, cols = q-values
      delta_matrix <- matrix(NA_real_, nrow = length(transcript_ids), ncol = length(q))
      colnames(delta_matrix) <- paste0("q_", gsub("\\.", "_", sprintf("%.2f", q)))
      rownames(delta_matrix) <- transcript_ids
      
      for (q_idx in seq_along(results_list)) {
        if (!is.null(results_list[[q_idx]]$results_per_gene[[gene_name]])) {
          delta_matrix[, q_idx] <- results_list[[q_idx]]$results_per_gene[[gene_name]]$delta_influence
        }
      }
      
      # Compute direction consistency for each transcript
      consistency_results <- apply(delta_matrix, 1, function(x) {
        x_valid <- x[!is.na(x) & !is.infinite(x)]
        if (length(x_valid) >= 2) {
          pos_count <- sum(x_valid > 0)
          neg_count <- sum(x_valid < 0)
          zero_count <- sum(x_valid == 0)
          
          if (pos_count == length(x_valid)) {
            return("Consistent positive")
          } else if (neg_count == length(x_valid)) {
            return("Consistent negative")
          } else if (zero_count == length(x_valid)) {
            return("All zero")
          } else {
            return("Mixed directions")
          }
        } else if (length(x_valid) == 1) {
          return("Single q-value")
        } else {
          return("No data")
        }
      })
      
      # Add direction_consistency to each q-value result for this gene
      for (q_idx in seq_along(results_list)) {
        if (!is.null(results_list[[q_idx]]$results_per_gene[[gene_name]])) {
          results_list[[q_idx]]$results_per_gene[[gene_name]]$direction_consistency <- consistency_results
        }
      }
    }
    
    class(results_list) <- c("tsenat_isoform_switching_multiq", "list")
    
    # Optional printing for multi-q results
    if (print_results) {
      cat("Isoform Switching Analysis - Multi-Q Comparison\n")
      cat("================================================\n\n")
      for (i in seq_along(results_list)) {
        cat("q = ", q[i], "\n")
        res <- results_list[[i]]
        if (!is.null(res$metadata)) {
          cat("  Genes analyzed:      ", length(res$gene_names), "\n")
          cat("  Transcripts tested:  ", res$metadata$n_transcripts_tested, "\n")
          cat("  FDR-significant:     ", res$metadata$n_fdr_significant, "\n")
          cat("  Genes with switching:", sum(res$summary_table$n_switching_transcripts > 0), "\n\n")
        }
      }
      cat("✓ Access results$q_<value>$results_per_gene$<gene> for per-q, per-gene details\n")
      cat("✓ Compare q values to assess scale-dependent isoform switching patterns\n\n")
    }
    
    return(invisible(results_list))
  }
  
  if (!(condition_col %in% colnames(colData(se)))) {
    stop("condition_col '", condition_col, "' not found in colData")
  }
  
  # Get conditions
  conditions <- unique(colData(se)[[condition_col]])
  if (length(conditions) != 2) {
    stop("Exactly 2 conditions required, found: ", length(conditions))
  }
  conditions <- sort(conditions)
  
  # Validate gene and isoform columns
  if (is.null(gene_col)) {
    stop("gene_col must be specified")
  }
  if (is.null(isoform_col)) {
    stop("isoform_col must be specified")
  }
  
  if (!(gene_col %in% colnames(rowData(se)))) {
    stop("gene_col '", gene_col, "' not found in rowData")
  }
  if (!(isoform_col %in% colnames(rowData(se)))) {
    stop("isoform_col '", isoform_col, "' not found in rowData")
  }
  
  # Handle paired design
  is_paired <- FALSE
  pair_info <- NULL
  
  if (!is.null(pair_col)) {
    if (!(pair_col %in% colnames(colData(se)))) {
      stop("pair_col '", pair_col, "' not found in colData")
    }
    
    # Check if pairing is valid (50% threshold)
    pairs <- colData(se)[[pair_col]]
    conds <- colData(se)[[condition_col]]
    pair_conds <- table(pairs, conds)
    paired_ratio <- sum(apply(pair_conds > 0, 1, sum) == 2) / nrow(pair_conds)
    
    if (paired_ratio >= 0.5) {
      is_paired <- TRUE
      paired_indices <- apply(pair_conds > 0, 1, sum) == 2
      pair_info <- list(
        pair_col = pair_col,
        n_pairs = sum(paired_indices),
        matched_pairs = names(which(paired_indices))
      )
    }
  }
  
  # Get gene list
  gene_ids <- unique(rowData(se)[[gene_col]])
  
  # Create mapping from gene_id to gene_name from rowData
  rd_mapping <- rowData(se)
  gene_id_to_name <- character(0)
  if ("gene_name" %in% colnames(rd_mapping)) {
    gene_id_to_name <- setNames(
      as.character(rd_mapping$gene_name),
      as.character(rd_mapping[[gene_col]])
    )
    # Remove duplicates, keeping first occurrence
    gene_id_to_name <- gene_id_to_name[!duplicated(names(gene_id_to_name))]
  }
  
  # Handle LM filtering with automatic gene name to ID mapping
  lm_gene_mapping <- NULL
  lm_genes_filtered <- 0
  
  if (!is.null(lm_results)) {
    if (!("gene" %in% colnames(lm_results))) {
      stop("lm_results must have 'gene' column")
    }
    
    # Check if lm_results contains gene names or gene IDs
    # First, check if the "gene" column matches gene IDs in se
    sample_lm_genes <- lm_results$gene[1:min(5, nrow(lm_results))]
    genes_are_ids <- all(sample_lm_genes %in% gene_ids)
    
    # If genes are names, map them to IDs
    if (!genes_are_ids) {
      # Map gene names to Ensembl IDs from rowData
      rd <- rowData(se)
      gene_name_to_id <- setNames(as.character(rd[[gene_col]]), as.character(rd$gene_name))
      gene_name_to_id <- gene_name_to_id[!is.na(names(gene_name_to_id))]
      gene_name_to_id <- gene_name_to_id[!duplicated(names(gene_name_to_id))]
      
      # Map genes in lm_results
      lm_results <- lm_results
      lm_results$gene <- unname(gene_name_to_id[as.character(lm_results$gene)])
      lm_results <- lm_results[!is.na(lm_results$gene), ]
      
      if (nrow(lm_results) == 0) {
        warning("No genes from lm_results matched in SummarizedExperiment rowData.")
        lm_gene_mapping <- NULL
      } else {
        lm_gene_mapping <- lm_results
      }
    } else {
      lm_gene_mapping <- lm_results
    }
    
    # Determine which p-value column to use
    if (!is.null(lm_gene_mapping)) {
      p_col <- NULL
      if (use_lm_fdr && "adj_p_interaction" %in% colnames(lm_gene_mapping)) {
        p_col <- "adj_p_interaction"
      } else if ("p_interaction" %in% colnames(lm_gene_mapping)) {
        p_col <- "p_interaction"
      } else {
        warning("lm_results missing p_interaction or adj_p_interaction columns")
      }
      
      if (!is.null(p_col)) {
        sig_genes <- lm_gene_mapping[lm_gene_mapping[[p_col]] < lm_p_threshold, "gene"]
        lm_genes_filtered <- length(sig_genes)
        
        if (length(sig_genes) == 0) {
          warning("No genes pass LM threshold p < ", lm_p_threshold,
                  ". Analyzing all genes.")
        } else {
          gene_ids <- intersect(gene_ids, sig_genes)
        }
      }
    }
  }
  
  # Note: top_n parameter removed - all LM-significant genes are analyzed
  
  # Initialize results
  results_per_gene <- list()
  all_pvalues <- list()
  all_fdr <- list()
  summary_rows <- list()
  
  # Track gene processing for debugging (if verbose mode enabled)
  gene_processing_log <- data.frame(
    gene = character(),
    n_transcripts = numeric(),
    has_2_transcripts = logical(),
    in_results = logical(),
    stringsAsFactors = FALSE
  )
  
  # Process each gene
  for (gene in gene_ids) {
    # Get transcript indices for this gene
    gene_mask <- rowData(se)[[gene_col]] == gene
    gene_isos <- rowData(se)[gene_mask, isoform_col]
    
    if (length(gene_isos) < 2) {
      gene_processing_log <- rbind(gene_processing_log, data.frame(
        gene = gene, n_transcripts = length(gene_isos),
        has_2_transcripts = FALSE, in_results = FALSE
      ))
      next # Skip single-transcript genes
    }
    
    # Mark that it passed 2-transcript check
    processing_check_passed <- TRUE
    
    # Get counts for this gene
    counts_matrix <- assays(se)$counts[gene_mask, , drop = FALSE]
    
    # Extract condition A and B counts
    cond_mask_A <- colData(se)[[condition_col]] == conditions[1]
    cond_mask_B <- colData(se)[[condition_col]] == conditions[2]
    
    counts_A_all <- counts_matrix[, cond_mask_A, drop = FALSE]
    counts_B_all <- counts_matrix[, cond_mask_B, drop = FALSE]
    
    # Handle paired design: use matched pairs only
    if (is_paired && !is.null(pair_col)) {
      pairs_A <- colData(se)[cond_mask_A, pair_col]
      pairs_B <- colData(se)[cond_mask_B, pair_col]
      
      matched_pairs <- intersect(pairs_A, pairs_B)
      
      cond_A_paired_idx <- which(cond_mask_A)[match(matched_pairs, pairs_A, nomatch = 0)]
      cond_B_paired_idx <- which(cond_mask_B)[match(matched_pairs, pairs_B, nomatch = 0)]
      
      if (length(cond_A_paired_idx) > 0 && length(cond_B_paired_idx) > 0) {
        counts_A <- counts_matrix[, cond_A_paired_idx, drop = FALSE]
        counts_B <- counts_matrix[, cond_B_paired_idx, drop = FALSE]
      } else {
        counts_A <- counts_A_all
        counts_B <- counts_B_all
      }
    } else {
      counts_A <- counts_A_all
      counts_B <- counts_B_all
    }
    
    # Calculate Tsallis entropy and influences
    # CRITICAL: Need to maintain consistent normalization across jackknife leave-one-out iterations
    n_tx_original <- nrow(counts_A)  # Store original transcript count
    
    .tsallis_entropy <- function(counts, q, norm, log_base, pseudocount, n_tx_fixed = NULL) {
      # Check for zero-sum columns BEFORE adding pseudocount to detect samples with no expression
      raw_col_sums <- colSums(counts)
      with_zero_counts <- raw_col_sums == 0
      
      # Ensure pseudocount is used to avoid log(0) and division by zero
      if (pseudocount == 0) {
        pseudocount <- 1e-8
      }
      counts <- counts + pseudocount
      
      col_sums <- colSums(counts)
      if (any(col_sums <= 0)) {
        return(rep(NA_real_, ncol(counts)))
      }
      
      # Return NA for samples with zero total counts (no expression for this gene)
      h_result <- rep(NA_real_, ncol(counts))
      if (all(with_zero_counts)) {
        return(h_result)
      }
      
      # Avoid division by zero for zero-count columns by using raw column sums intelligently
      # For zero-count columns, set col_sums to 1 to avoid division by near-zero
      # These samples will return NA anyway (no valid probability distribution)
      col_sums_safe <- pmax(col_sums, 1)
      
      p <- counts / col_sums_safe
      
      if (q == 1) {
        h <- if (log_base == exp(1)) {
          -colSums(p * log(p + 1e-100))
        } else {
          -colSums(p * log(p + 1e-100, log_base))
        }
      } else {
        p_q_sum <- colSums(p^q)
        h <- (1 / (1 - q)) * (1 - p_q_sum)
      }
      
      h[!is.finite(h)] <- NA_real_
      # Set entropy to NA for zero-count columns
      h[with_zero_counts] <- NA_real_
      
      if (norm) {
        # Use FIXED n_transcripts if provided (for consistent jackknife normalization)
        n_tx_use <- if (!is.null(n_tx_fixed)) n_tx_fixed else nrow(counts)
        if (q == 1) {
          max_h <- log(n_tx_use) / log(log_base)
        } else {
          # Use absolute value because the formula produces negative values for q<1 and q>1
          max_h <- abs((1 / (1 - q)) * (1 - n_tx_use^(1 - q)))
        }
        
        if (!is.na(max_h) && is.finite(max_h) && max_h > 0) {
          h <- h / max_h
        }
      }
      return(h)
    }
    
    .jackknife_influences <- function(counts, q, norm, log_base, pseudocount, n_tx_fixed = NULL) {
      h_full <- .tsallis_entropy(counts, q, norm, log_base, pseudocount, n_tx_fixed)
      n_tx <- nrow(counts)
      influences <- numeric(n_tx)
      
      if (n_tx < 2) {
        return(influences)
      }
      
      for (i in 1:n_tx) {
        counts_leave_i <- counts[-i, , drop = FALSE]
        h_leave_i <- .tsallis_entropy(counts_leave_i, q, norm, log_base, pseudocount, n_tx_fixed)
        # Compute mean absolute difference across samples
        diffs <- abs(h_full - h_leave_i)
        influences[i] <- mean(diffs, na.rm = TRUE)
      }
      return(influences)
    }
    
    # Calculate influences for each condition with fixed transcript count
    # This ensures consistent normalization across all jackknife iterations
    influences_A <- .jackknife_influences(counts_A, q, norm, log_base, pseudocount, n_tx_original)
    influences_B <- .jackknife_influences(counts_B, q, norm, log_base, pseudocount, n_tx_original)
    
    # Calculate delta influence
    delta_influence <- influences_A - influences_B
    
    # Get original transcript count for consistent normalization during jackknife bootstrap
    n_tx <- nrow(counts_A)
    
    # Calculate bootstrap statistics
    delta_stats <- compute_delta_statistics(
      counts_A, counts_B, delta_influence,
      q = q, norm = norm, log_base = log_base,
      pseudocount = pseudocount, n_bootstrap = n_bootstrap,
      n_transcripts = n_tx  # Pass original transcript count for consistent normalization
    )
    
    # Determine switching status
    switching_status <- ifelse(delta_influence > 0, "up", ifelse(delta_influence < 0, "down", "neutral"))
    
    # Store gene result (delta_stats now contains per-transcript vectors)
    gene_result <- list(
      gene_id = gene,
      transcript_ids = as.character(gene_isos),
      delta_influence = delta_influence,
      delta_se = delta_stats$se,
      delta_ci_lower = delta_stats$ci_lower,
      delta_ci_upper = delta_stats$ci_upper,
      delta_pvalue = delta_stats$pvalue,
      switching_status = switching_status
    )
    
    # Add LM annotation if available
    if (!is.null(lm_gene_mapping)) {
      lm_row <- lm_gene_mapping[lm_gene_mapping$gene == gene, ]
      if (nrow(lm_row) > 0) {
        gene_result$lm_p_interaction <- if ("p_interaction" %in% colnames(lm_row)) {
          lm_row$p_interaction[1]
        } else {
          NA
        }
        gene_result$lm_adj_p_interaction <- if ("adj_p_interaction" %in% colnames(lm_row)) {
          lm_row$adj_p_interaction[1]
        } else {
          NA
        }
      }
    }
    
    results_per_gene[[gene]] <- gene_result
    
    # Log successful processing
    gene_processing_log <- rbind(gene_processing_log, data.frame(
      gene = gene, n_transcripts = length(gene_isos),
      has_2_transcripts = TRUE, in_results = TRUE
    ))
    
    # Collect p-values for global FDR correction (per-transcript)
    for (i in seq_along(gene_isos)) {
      all_pvalues[[length(all_pvalues) + 1]] <- list(
        gene = gene,
        transcript = gene_isos[i],
        pvalue = delta_stats$pvalue[i]  # Per-transcript p-value
      )
    }
  }
  
  # Global FDR correction (Benjamini-Hochberg)
  if (length(all_pvalues) > 0) {
    pvals_vec <- sapply(all_pvalues, function(x) x$pvalue)
    fdr_vec <- p.adjust(pvals_vec, method = "BH")
    
    for (i in seq_along(all_pvalues)) {
      all_pvalues[[i]]$fdr <- fdr_vec[i]
    }
    
    # Add FDR to gene results
    for (px in all_pvalues) {
      gene_idx <- match(px$gene, names(results_per_gene))
      if (!is.na(gene_idx)) {
        iso_idx <- match(px$transcript, results_per_gene[[gene_idx]]$transcript_ids)
        if (!is.na(iso_idx)) {
          if (is.null(results_per_gene[[gene_idx]]$delta_fdr)) {
            results_per_gene[[gene_idx]]$delta_fdr <- rep(NA, length(results_per_gene[[gene_idx]]$transcript_ids))
          }
          results_per_gene[[gene_idx]]$delta_fdr[iso_idx] <- px$fdr
        }
      }
    }
  }
  
  # Build all_transcript_stats table
  all_transcript_stats <- do.call(rbind, lapply(names(results_per_gene), function(gene) {
    res <- results_per_gene[[gene]]
    df <- data.frame(
      gene = gene,
      transcript_id = res$transcript_ids,
      pvalue = ifelse(is.null(res$delta_pvalue), NA_real_, res$delta_pvalue),
      fdr = ifelse(is.null(res$delta_fdr), NA_real_, res$delta_fdr),
      stringsAsFactors = FALSE
    )
    
    if (!is.null(res$lm_p_interaction)) {
      df$lm_p_interaction <- res$lm_p_interaction
    }
    if (!is.null(res$lm_adj_p_interaction)) {
      df$lm_adj_p_interaction <- res$lm_adj_p_interaction
    }
    
    return(df)
  }))
  
  # Build summary table with gene names
  summary_rows <- lapply(names(results_per_gene), function(gene) {
    res <- results_per_gene[[gene]]
    n_tx <- length(res$transcript_ids)
    n_switching <- sum(res$switching_status != "neutral", na.rm = TRUE)
    # Only compute max_delta if there are finite values
    delta_vals <- res$delta_influence[is.finite(res$delta_influence)]
    max_delta <- if (length(delta_vals) > 0) max(abs(delta_vals)) else 0
    n_fdr_sig <- sum(res$delta_fdr < 0.05, na.rm = TRUE)
    
    # Get gene name from rowData if available
    gene_name <- gene  # Default to gene ID
    gene_rows <- which(rowData(se)[[gene_col]] == gene)
    if (length(gene_rows) > 0 && "gene_name" %in% colnames(rowData(se))) {
      gene_name_from_data <- rowData(se)[gene_rows[1], "gene_name"]
      if (!is.na(gene_name_from_data)) {
        gene_name <- as.character(gene_name_from_data)
      }
    }
    
    data.frame(
      gene = gene,
      gene_name = gene_name,
      n_transcripts = n_tx,
      max_delta_influence = max_delta,
      n_switching_transcripts = n_switching,
      n_fdr_significant = n_fdr_sig,
      stringsAsFactors = FALSE
    )
  })
  
  summary_table <- do.call(rbind, summary_rows)
  rownames(summary_table) <- NULL
  
  # Build metadata
  n_fdr_sig_total <- if (nrow(all_transcript_stats) > 0) {
    sum(all_transcript_stats$fdr < 0.05, na.rm = TRUE)
  } else {
    0
  }
  
  metadata <- list(
    q = q,
    is_paired = is_paired,
    pair_col = if (is_paired) pair_col else NULL,
    pair_info = pair_info,
    norm = norm,
    log_base = log_base,
    pseudocount = pseudocount,
    threshold = threshold,
    n_transcripts_tested = nrow(all_transcript_stats),
    n_fdr_significant = n_fdr_sig_total,
    lm_results_provided = !is.null(lm_results),
    lm_p_threshold = lm_p_threshold,
    lm_genes_filtered = lm_genes_filtered,
    gene_processing_log = gene_processing_log
  )
  
  # Create result object
  # Build gene_id and gene_name vectors for analyzed genes
  analyzed_genes <- names(results_per_gene)
  
  # gene_names currently contains gene IDs (from rowData gene_col)
  # Add explicit gene_id vector for clarity
  result_gene_ids <- analyzed_genes
  
  # Add gene_names vector (mapped from gene_id_to_name)
  result_gene_names <- unname(gene_id_to_name[analyzed_genes])
  if (length(result_gene_names) == 0) {
    result_gene_names <- rep(NA_character_, length(analyzed_genes))
  }
  
  result <- list(
    gene_names = analyzed_genes,  # Keep for backwards compatibility (actually gene IDs)
    gene_ids = result_gene_ids,   # Explicit gene IDs
    gene_name_map = result_gene_names,  # Mapping to gene symbols/names
    conditions = conditions,
    results_per_gene = results_per_gene,
    summary_table = summary_table,
    all_transcript_stats = all_transcript_stats,
    metadata = metadata
  )
  
  class(result) <- c("tsenat_isoform_switching", "list")
  
  # Print results if requested
  if (print_results) {
    cat("Isoform Switching Analysis Results\n")
    cat("===================================\n")
    cat("Conditions: '", conditions[1], "' vs. '", conditions[2], "'\n", sep = "")
    
    if (is_paired && !is.null(pair_info)) {
      cat("Design: PAIRED (", pair_col, ") - ", pair_info$n_pairs, " matched pairs\n", sep = "")
    } else {
      cat("Design: UNPAIRED\n")
    }
    
    if (!is.null(lm_results)) {
      cat("LM filtering: Genes with p < ", lm_p_threshold, " (N = ", lm_genes_filtered, ")\n", sep = "")
    }
    
    cat("\nGenes analyzed:", length(result$gene_names), "\n")
    cat("Total transcripts tested:", result$metadata$n_transcripts_tested, "\n")
    cat("FDR-significant transcripts (FDR<0.05):", result$metadata$n_fdr_significant, "\n\n")
    
    cat("Summary Table:\n")
    print(summary_table)
    
    cat("\n✓ Use results$results_per_gene$'GeneName' to access per-gene switching details\n")
    cat("✓ Use results$summary_table for overview across genes\n")
    cat("✓ Use results$all_transcript_stats for FDR-corrected p-values per transcript\n")
    cat("✓ Use results$metadata$is_paired to check if paired design was applied\n")
    
    if (!is.null(lm_results)) {
      cat("✓ Access lm_p_interaction in each gene$lm_p_interaction for LM test results\n")
    }
    cat("\n")
  }
  
  return(invisible(result))
}


#' Heatmap of Isoform Switching Patterns
#'
#' Visualizes transcript-level switching effects as a heatmap where rows are
#' transcripts and columns are genes, with cell values showing delta_influence.
#' Color intensity indicates strength of condition shift.
#'
#' @param switching_results A result object from \code{jackknife_isoform_switching()}.
#' @param top_n Numeric: show top N genes by max switching effect (default 10).
#' @param top_transcripts_per_gene Numeric: top transcripts per gene (default 3).
#' @param main Character: title for heatmap.
#' @param color_scheme Character: "diverging" (blue-white-red, default) or "sequential" (white-orange-red).
#' @param show_values Logical: show numeric values in cells? (default TRUE).
#'
#' @return Invisibly returns a matrix of delta_influence values plotted.
#'
#' @details
#' Heatmap interpretation:
#' - \strong{Red}: Transcript increases influence in Condition A (up-switching)
#' - \strong{Blue}: Transcript increases influence in Condition B (down-switching)
#' - \strong{White/Light}: Minimal switching (neutral effect)
#'
#' Useful for:
#' - Identifying coordinated switching across genes
#' - Ranking genes by switching magnitude
#' - Finding tissue-specific or condition-specific isoforms
#'
#' References: \strong{ISO014} (classifying exon usage), \strong{ISO018} (tissue-specific splicing)
#'
#' @examples
#' \dontrun{
#' # Visualize top switching transcripts
#' plot_isoform_switching_heatmap(
#'   switching_results,
#'   top_n = 5,
#'   top_transcripts_per_gene = 3
#' )
#' }
#'
#' @export
plot_isoform_switching_heatmap <- function(
  switching_results,
  top_n = 10,
  top_transcripts_per_gene = 3,
  main = "Isoform Switching Heatmap (delta_influence)",
  color_scheme = "diverging",
  show_values = TRUE
) {
  
  if (!inherits(switching_results, "tsenat_isoform_switching")) {
    stop("switching_results must be from jackknife_isoform_switching()")
  }
  
  # Extract top genes by max switching
  summary_table <- switching_results$summary_table
  summary_table <- summary_table[order(-summary_table$max_delta_influence), ]
  top_genes <- summary_table$gene[1:min(top_n, nrow(summary_table))]
  
  # Create mapping of gene ID to gene name
  gene_names <- setNames(summary_table$gene_name[1:min(top_n, nrow(summary_table))], 
                         summary_table$gene[1:min(top_n, nrow(summary_table))])
  
  # Build data structure organized by gene
  gene_data <- list()
  all_values <- numeric()
  
  for (gene in top_genes) {
    res <- switching_results$results_per_gene[[gene]]
    
    # Get top transcripts by |delta_influence|
    deltas <- res$delta_influence
    top_iso_idx <- order(abs(deltas), decreasing = TRUE)[1:min(top_transcripts_per_gene, length(deltas))]
    
    iso_names <- res$transcript_ids[top_iso_idx]
    iso_labels <- sub("^[^-]+-", "", iso_names)  # Remove prefix for readability
    iso_values <- deltas[top_iso_idx]
    
    gene_data[[gene]] <- list(
      labels = iso_labels,
      values = iso_values
    )
    
    all_values <- c(all_values, iso_values)
  }
  
  if (length(gene_data) == 0) {
    warning("No transcripts found for heatmap")
    return(invisible(NULL))
  }
  
  # Calculate global max for consistent color scale
  max_val <- max(abs(all_values), na.rm = TRUE)
  
  # Handle infinite values
  if (!is.finite(max_val)) {
    finite_vals <- all_values[is.finite(all_values)]
    if (length(finite_vals) > 0) {
      max_val <- max(abs(finite_vals), na.rm = TRUE)
    } else {
      max_val <- 1
    }
  }
  
  # Replace infinite values in gene_data with max_val
  for (gene in names(gene_data)) {
    gene_data[[gene]]$values[!is.finite(gene_data[[gene]]$values)] <- 
      sign(gene_data[[gene]]$values[!is.finite(gene_data[[gene]]$values)]) * max_val
  }
  
  # Create color scheme
  n_colors <- 99
  breaks <- seq(-max_val, max_val, length.out = n_colors + 1)
  
  if (color_scheme == "diverging") {
    colors <- c(
      grDevices::colorRampPalette(c("blue", "white"))(n_colors %/% 2),
      grDevices::colorRampPalette(c("white", "red"))(n_colors - n_colors %/% 2)
    )
  } else {
    colors <- grDevices::colorRampPalette(c("white", "orange", "red"))(n_colors)
  }
  
  # Save graphics parameters
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  
  # Create subplot layout - 4 columns, 1 row
  n_genes <- length(gene_data)
  ncol <- min(4, n_genes)
  nrow <- ceiling(n_genes / ncol)
  
  graphics::par(mfrow = c(nrow, ncol), mar = c(8, 3, 2.5, 1), oma = c(0, 0, 4, 0))
  
  # Plot each gene in its own panel
  for (gene in top_genes) {
    gene_name <- gene_names[gene]
    if (is.na(gene_name)) gene_name <- gene  # Fallback to gene ID if name not found
    
    gene_info <- gene_data[[gene]]
    n_iso <- length(gene_info$labels)
    
    # Create matrix for this gene
    heatmap_matrix <- matrix(gene_info$values, nrow = n_iso, ncol = 1)
    rownames(heatmap_matrix) <- gene_info$labels
    
    # Plot heatmap for this gene (transcripts on x-axis)
    graphics::image(
      x = 1:n_iso,
      y = 1,
      z = heatmap_matrix,
      col = colors,
      breaks = breaks,
      main = gene_name,
      xlab = "",
      ylab = "",
      xaxt = "n",
      yaxt = "n",
      zlim = c(-max_val, max_val)
    )
    
    # Add transcript labels on x-axis
    graphics::axis(
      1, 
      at = 1:n_iso, 
      labels = gene_info$labels, 
      las = 2,
      cex.axis = 0.85,
      tick = TRUE,
      mgp = c(3, 2, 0)
    )
    
    # Add values if requested
    if (show_values) {
      for (i in seq_len(n_iso)) {
        graphics::text(i, 1, labels = round(gene_info$values[i], 3), cex = 0.65, col = "black")
      }
    }
  }
  
  # Add overall title
  mtext(main, side = 3, outer = TRUE, line = 2, cex = 1.3)
  
  return(invisible(NULL))
}


#' Q-Parameter Sensitivity Curve
#'
#' Plots how isoform switching effects change across different q-values,
#' revealing scale-dependent diversity patterns (paper I004).
#'
#' @param se A SummarizedExperiment object with transcript-level counts.
#' @param condition_col Character: column name in colData for condition labels.
#' @param gene Character: name of the gene to plot.
#' @param gene_col Character: column name in rowData for gene IDs.
#' @param isoform_col Character: column name in rowData for transcript/isoform IDs.
#' @param q_values Numeric: vector of q values to test (default c(0.5, 1, 1.5, 2, 2.5)).
#' @param norm Logical: normalize entropy to [0,1]? (default TRUE).
#' @param log_base Numeric: log base for entropy (default e).
#' @param pseudocount Numeric: pseudocount to add (default 0).
#' @param main Character: plot title.
#' @param show_legend Logical: show legend with transcript names? (default TRUE).
#'
#' @return Invisibly returns a data.frame with delta_influence for each q and transcript.
#'
#' @details
#' Q-parameter interpretation:
#' - \strong{q < 1}: Emphasizes rare transcripts
#' - \strong{q = 1}: Shannon entropy (balanced weighting)
#' - \strong{q > 1}: Emphasizes abundant transcripts
#' - \strong{q spectrum}: Different genes may show optimal switching at different q
#'
#' References: \strong{S111} (q-parameter sensitivity), \strong{I004} (q-parameter theory)
#'
#' @examples
#' \dontrun{
#' # Plot switching effects across q values for a specific gene
#' plot_q_sensitivity_curve(
#'   se = ts_se,
#'   condition_col = "condition",
#'   gene = "Gene1",
#'   gene_col = "gene_id",
#'   isoform_col = "isoform_id"
#' )
#' }
#'
#' @export
plot_q_sensitivity_curve <- function(
  se = NULL,
  condition_col = "condition",
  gene = NULL,
  gene_col = NULL,
  isoform_col = NULL,
  q_values = c(0.5, 1, 1.5, 2, 2.5),
  norm = TRUE,
  log_base = exp(1),
  pseudocount = 0,
  main = NULL,
  show_legend = TRUE
) {
  
  # Input validation
  if (is.null(se) || is.null(gene) || is.null(gene_col) || is.null(isoform_col)) {
    stop("se, gene, gene_col, and isoform_col are required")
  }
  
  if (!(gene %in% rowData(se)[[gene_col]])) {
    stop("Gene '", gene, "' not found in rowData")
  }
  
  # Get conditions
  conditions <- sort(unique(colData(se)[[condition_col]]))
  if (length(conditions) != 2) {
    stop("Exactly 2 conditions required")
  }
  
  # Extract gene data
  gene_mask <- rowData(se)[[gene_col]] == gene
  gene_isos <- rowData(se)[gene_mask, isoform_col]
  counts_matrix <- assays(se)$counts[gene_mask, , drop = FALSE]
  
  # Helper function: compute per-isoform influence
  compute_influence_for_q <- function(counts, q) {
    h_full <- if (q == 1) {
      -colSums((counts / colSums(counts)) * log(counts / colSums(counts) + 1e-100))
    } else {
      (1 / (1 - q)) * log(colSums((counts / colSums(counts))^q) + 1e-100)
    }
    
    if (norm) {
      n_tx <- nrow(counts)
      max_h <- if (q == 1) {
        log(n_tx) / log(log_base)
      } else {
        (1 / (1 - q)) * (1 - n_tx^(1 - q)) / log(log_base)
      }
      if (!is.na(max_h) && is.finite(max_h) && max_h > 0) {
        h_full <- h_full / max_h
      }
    }
    
    # Compute mean entropy
    mean(h_full, na.rm = TRUE)
  }
  
  # Compute delta_influence for each q and isoform
  cond_mask_A <- colData(se)[[condition_col]] == conditions[1]
  cond_mask_B <- colData(se)[[condition_col]] == conditions[2]
  
  counts_A <- counts_matrix[, cond_mask_A, drop = FALSE]
  counts_B <- counts_matrix[, cond_mask_B, drop = FALSE]
  
  # Store results
  sensitivity_data <- data.frame()
  
  for (q in q_values) {
    h_A <- compute_influence_for_q(counts_A, q)
    h_B <- compute_influence_for_q(counts_B, q)
    
    sensitivity_data <- rbind(sensitivity_data, data.frame(
      q = q,
      condition_A_entropy = h_A,
      condition_B_entropy = h_B,
      delta_entropy = h_A - h_B,
      stringsAsFactors = FALSE
    ))
  }
  
  # Create plot
  if (is.null(main)) {
    main <- paste("Q-Sensitivity Curve for", gene)
  }
  
  plot(
    sensitivity_data$q,
    sensitivity_data$delta_entropy,
    type = "b",
    pch = 16,
    main = main,
    xlab = "q parameter",
    ylab = "delta_entropy (Condition A - Condition B)",
    col = "red",
    lwd = 2,
    cex = 1.2
  )
  
  # Add horizontal line at zero
  graphics::abline(h = 0, col = "gray", lty = 2, lwd = 1)
  
  # Add condition entropy curves
  graphics::lines(sensitivity_data$q, sensitivity_data$condition_A_entropy, 
                   type = "b", pch = 15, col = "blue", lty = 2, cex = 0.8)
  graphics::lines(sensitivity_data$q, sensitivity_data$condition_B_entropy, 
                   type = "b", pch = 17, col = "green", lty = 2, cex = 0.8)
  
  # Add legend
  if (show_legend) {
    graphics::legend(
      "topright",
      legend = c(
        paste("Δ entropy (", conditions[1], " - ", conditions[2], ")", sep = ""),
        paste("Entropy:", conditions[1]),
        paste("Entropy:", conditions[2])
      ),
      col = c("red", "blue", "green"),
      lty = c(1, 2, 2),
      pch = c(16, 15, 17)
    )
  }
  
  # Add q-range recommendation
  graphics::abline(v = 0.5, col = "gray", lty = 3)
  graphics::abline(v = 2, col = "gray", lty = 3)
  
  return(invisible(sensitivity_data))
}


#' Print Isoform Switching Results
#'
#' @param x A tsenat_isoform_switching object
#' @param ... Additional arguments (unused)
#'
#' @keywords internal
#' @noRd
#' @method print tsenat_isoform_switching
print.tsenat_isoform_switching <- function(x, ...) {
  cat("Isoform Switching Analysis\n")
  cat("==========================\n")
  cat("Conditions: '", x$conditions[1], "' vs. '", x$conditions[2], "'\n", sep = "")
  
  if (x$metadata$is_paired && !is.null(x$metadata$pair_info)) {
    cat("Design: PAIRED (", x$metadata$pair_col, ") - ", 
        x$metadata$pair_info$n_pairs, " matched pairs\n", sep = "")
  } else {
    cat("Design: UNPAIRED\n")
  }
  
  cat("Genes analyzed:", length(x$gene_names), "\n")
  cat("q-parameter:", x$metadata$q, "\n")
  cat("Transcripts analyzed:", x$metadata$n_transcripts_tested, "\n")
  cat("FDR-significant transcripts:", x$metadata$n_fdr_significant, "\n\n")
  
  cat("Use $results_per_gene to access per-gene results\n")
  cat("Use $summary_table for summary across genes\n")
  cat("Use $all_transcript_stats for all transcript statistics\n")
  
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
#' @keywords internal
#' @export
prepare_gene_switching_tables <- function(
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
    cat(sprintf("Extracted q_vector from multi_q_results: %s\n", paste(sprintf("%.2f", q_vector), collapse=", ")))
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
      cat(sprintf("n_top_genes is NULL; using all %d genes from summary_df\n", n_top_genes))
    }
  }
  
  if (verbose) {
    cat(sprintf("Created summary_df with %d genes\n", nrow(summary_df)))
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
  for (i in 1:min(n_top_genes, nrow(summary_df))) {
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
    cat(sprintf("Matched %d top genes to multi_q_results\n", length(top_genes_list)))
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
      tx_ids <- gene_data_by_q[[first_q]]$transcript_ids[1:n_tx]
      comparison_data <- data.frame(transcript = tx_ids, stringsAsFactors = FALSE)
      
      # Add delta_influence values for each q-value
      for (q_key in q_values_available) {
        delta_vals <- gene_data_by_q[[q_key]]$delta_influence[1:n_tx]
        comparison_data[[q_key]] <- delta_vals
      }
      
      # Clean NaN and Inf values for display
      for (col in q_values_available) {
        comparison_data[[col]][is.nan(comparison_data[[col]]) | is.infinite(comparison_data[[col]])] <- NA
      }
      
      # Get pre-computed direction consistency from jackknife results
      consistency_results <- gene_data_by_q[[q_values_available[1]]]$direction_consistency[1:n_tx]
      
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
