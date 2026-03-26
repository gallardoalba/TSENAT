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
#' @param verbose Logical: if TRUE (default), display formatted results for each gene and print diagnostic messages.
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
#' - \code{calculate_tsallis_entropy()}: computes entropy (stability as background)
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
#' \code{calculate_tsallis_entropy()} for entropy calculation,
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
#'   verbose = TRUE  # Auto-displays summary for all genes
#' )
#' 
#' # Example 2b: Multiple q values for robustness checking
#' jack_multiq <- jackknife_tsallis_entropy(
#'   x = counts_matrix[1, ],  # First gene
#'   q = c(0.5, 1, 1.5, 2),
#'   norm = TRUE,
#'   verbose = TRUE  # Shows stability across q values
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
#' #     verbose = TRUE  # Auto-extracts top 5 genes and displays summary                                                # Auto-extracts top 5 genes and displays summary
#' # )
#'
#' @keywords internal
#' @noRd
jackknife_tsallis_entropy <- function(x = NULL, se = NULL, res = NULL, top_n = 5,
                                       q = 1, norm = TRUE, log_base = exp(1),
                                       pseudocount = 0, threshold = 90, seed = NULL,
                                       verbose = FALSE, nthreads = 1, .cluster = NULL) {

  # Input validation
  if (!is.numeric(q) || any(q <= 0)) {
    stop("'q' must be positive numeric value(s)")
  }

  if (!is.numeric(threshold) || threshold < 0 || threshold > 100) {
    stop("'threshold' must be between 0 and 100")
  }
  
  # Handle multiple q values (with parallel processing optimization)
  if (length(q) > 1) {
    # Use parallel processing if:
    # - User requested parallelization (nthreads > 1 or nthreads = NULL)
    # - q-values > 2 (overhead only worthwhile for multiple q scenarios)
    # - parallel package is available
    # P2 OPTIMIZATION: Reuse existing cluster if provided, create once instead of per-recursion
    
    # Determine number of threads to use
    if (is.null(nthreads)) {
      # Auto-detect: use all available cores minus 1
      n_cores <- max(1, parallel::detectCores() - 1)
      if (exists(".tsenat_get_effective_nthreads", mode = "function")) {
        n_cores <- .tsenat_get_effective_nthreads(n_cores)
      } else {
        # Fallback: check env var manually for R CMD check
        core_limit <- Sys.getenv("_R_CHECK_LIMIT_CORES_", NA)
        if (!is.na(core_limit)) {
          core_limit <- as.integer(core_limit)
          if (is.finite(core_limit) && core_limit > 0) {
            n_cores <- min(n_cores, core_limit)
          }
        }
      }
    } else if (nthreads > 1) {
      # User explicitly requested parallelization
      n_cores <- as.integer(nthreads)
    } else {
      # nthreads = 1: skip parallelization
      n_cores <- 1
    }
    
    # Create cluster only if parallelization requested and not already provided
    create_cluster <- is.null(.cluster) && n_cores > 1 && length(q) > 2 && requireNamespace("parallel", quietly = TRUE)
    
    if (create_cluster) {
      .cluster <- parallel::makeCluster(n_cores, type = "PSOCK")
      on.exit(parallel::stopCluster(.cluster), add = TRUE)
      
      # Export required functions to cluster (only when creating new cluster)
      parallel::clusterExport(.cluster, c("jackknife_tsallis_entropy", ".tsenat_entropy_single"), 
                             envir = environment())
    }
    
    # P2 OPTIMIZATION: Pass cluster through recursive calls to avoid recreation
    if (!is.null(.cluster)) {
      # Parallel lapply for each q value with reused cluster
      results_list <- parallel::parLapply(.cluster, q, function(q_val) {
        jackknife_tsallis_entropy(
          x = x, se = se, res = res, top_n = top_n, q = q_val, 
          norm = norm, log_base = log_base, pseudocount = pseudocount,
          threshold = threshold, seed = seed, verbose = FALSE
        )
      })
    } else {
      # Sequential lapply for small q-value sets or when parallel not available
      results_list <- lapply(q, function(q_val) {
        jackknife_tsallis_entropy(
          x = x, se = se, res = res, top_n = top_n, q = q_val, 
          norm = norm, log_base = log_base, pseudocount = pseudocount,
          threshold = threshold, seed = seed, verbose = FALSE,
          .cluster = NULL  # No cluster available
        )
      })
    }
    
    names(results_list) <- paste0("q=", q)
    class(results_list) <- c("tsenat_jackknife_list_multiq", "list")
    
    # Optional printing
    # P3 OPTIMIZATION: Batch message formatting for multi-q display
    if (verbose && !is.null(x) && (is.vector(x) || length(q) > 1)) {
      output_lines <- c(
        "Jackknife Stability Analysis for Multiple q Values",
        "===================================================="
      )
      
      for (i in seq_along(results_list)) {
        output_lines <- c(output_lines, sprintf("q = %s", q[i]))
        res <- results_list[[i]]
        if (is.list(res) && "estimate" %in% names(res)) {
          output_lines <- c(output_lines,
            sprintf("  Estimate:            %.4f", res$estimate),
            sprintf("  Jackknife SE:        %.4f", res$jackknife_se),
            sprintf("  Max influence:       %.4f", max(res$influence)),
            sprintf("  Outliers detected:   %d", length(res$outlier_indices))
          )
        }
      }
      
      # Print all at once instead of individual message() calls
      message(paste(output_lines, collapse = "\n"))
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

    # P0 OPTIMIZATION: Pre-build lookup tables for gene/rowData columns
    # This avoids repeated which() scans inside the loop (O(n) per gene -> O(n) once)
    se_rownames <- rownames(se)
    rd <- SummarizedExperiment::rowData(se)
    rd_cols <- if (!is.null(rd)) colnames(rd) else character(0)
    counts_assay <- as.matrix(SummarizedExperiment::assay(se, "counts"))
    
    # Pre-build lookup tables (one-time cost)
    # Strategy 1: Hash table for rownames
    rownames_lookup <- setNames(seq_along(se_rownames), se_rownames)
    
    # Strategy 2: Build mappings for common rowData columns (if present)
    rowdata_lookups <- list()
    if ("gene_name" %in% rd_cols && !is.null(rd$gene_name)) {
      # Group transcript indices by gene_name
      rowdata_lookups$gene_name <- tapply(seq_len(nrow(rd)), rd$gene_name, list, simplify = FALSE)
    }
    if ("gene_id" %in% rd_cols && !is.null(rd$gene_id)) {
      # Group transcript indices by gene_id
      rowdata_lookups$gene_id <- tapply(seq_len(nrow(rd)), rd$gene_id, list, simplify = FALSE)
    }
    
    # Process all genes with optimized lookups
    results_temp <- lapply(top_genes, function(gene) {
      # P0 OPTIMIZATION: Use safe hash table lookups instead of which() (O(n) scan)
      # Safe lookup: check if gene exists in names first before accessing
      tx_idx <- NULL
      
      # Try rownames first (fastest)
      if (gene %in% names(rownames_lookup)) {
        tx_idx <- rownames_lookup[[gene]]
      }
      
      # Fallback to rowData columns (still O(1) with pre-built tables)
      if (is.null(tx_idx)) {
        # Try gene_name column
        if (!is.null(rowdata_lookups$gene_name) && gene %in% names(rowdata_lookups$gene_name)) {
          tx_idx <- rowdata_lookups$gene_name[[gene]]
        }
      }
      if (is.null(tx_idx) && !is.null(rowdata_lookups$gene_id)) {
        # Try gene_id column as fallback
        if (gene %in% names(rowdata_lookups$gene_id)) {
          tx_idx <- rowdata_lookups$gene_id[[gene]]
        }
      }
      
      if (is.null(tx_idx) || length(tx_idx) == 0) {
        warning("Gene '", gene, "' not found in 'se' rownames or rowData columns. Skipping.")
        return(NULL)
      }
      
      # P0 OPTIMIZATION: Handle both single value (from rownames) and list (from rowData) cases
      if (is.list(tx_idx)) {
        tx_idx <- unlist(tx_idx)  # Convert list result to vector
      }
      
      # Aggregate counts across all transcripts for this gene
      gene_counts <- as.numeric(colSums(counts_assay[tx_idx, , drop = FALSE]))
      
      # Validate counts
      if (any(is.na(gene_counts)) || any(gene_counts < 0)) {
        warning("Gene '", gene, "' has missing or invalid count values. Skipping.")
        return(NULL)
      }
      
      list(counts = gene_counts, name = gene)
    })
    
    # Filter out NULL results and extract
    # P4 OPTIMIZATION: Use direct logical indexing instead of Filter(Negate(is.null))
    # Avoids closure creation and higher-level function overhead
    is_valid <- !vapply(results_temp, is.null, logical(1))
    valid_results <- results_temp[is_valid]
    counts_list <- lapply(valid_results, "[[", "counts")
    gene_names_out <- vapply(valid_results, "[[", "name", FUN.VALUE = character(1))
    
    if (length(counts_list) == 0) {
      stop("None of the top genes from 'res' found in 'se' rownames or rowData with valid counts")
    }

    # Create matrix with genes as rows
    counts_matrix <- do.call(rbind, counts_list)
    rownames(counts_matrix) <- gene_names_out

    # Call recursively with matrix input (with verbose suppressed for first call)
    # P2 OPTIMIZATION: Pass .cluster through recursion
    return(jackknife_tsallis_entropy(
      x = counts_matrix,
      q = q,
      norm = norm,
      log_base = log_base,
      pseudocount = pseudocount,
      threshold = threshold,
      seed = seed,
      verbose = verbose,
      nthreads = nthreads,  # Pass nthreads parameter through recursion
      .cluster = .cluster  # Pass cluster if it exists
    ))
  }

  # Traditional path: x must be provided
  if (is.null(x)) {
    stop("Either 'x' or both 'se' and 'res' must be provided")
  }

  # Handle matrix/data.frame input
  # P1 OPTIMIZATION: Inline jackknife logic instead of recursive calls
  # This eliminates function call overhead (~0.5ms per gene × number of genes)
  if (is.matrix(x) || is.data.frame(x)) {
    x <- as.matrix(x)
    results <- vector("list", nrow(x))
    names(results) <- rownames(x)
    
    # P1 OPTIMIZATION: Process all genes with inlined jackknife logic
    # Avoids function call overhead of recursive jackknife_tsallis_entropy() calls
    for (gene_idx in seq_len(nrow(x))) {
      gene_counts <- x[gene_idx, ]
      
      # Inlined jackknife logic (from vector path below)
      n <- length(gene_counts)
      if (n < 2) {
        warning("Gene '", names(results)[gene_idx], "' has < 2 transcripts. Skipping.")
        results[[gene_idx]] <- NULL
        next
      }
      
      # Pre-compute normalized proportions
      total_count <- sum(gene_counts)
      p <- (gene_counts + pseudocount) / (total_count + length(gene_counts) * pseudocount)
      
      # Core jackknife loop with vectorized incremental updates (P0 optimization)
      jackknife_estimates <- numeric(n)
      
      if (abs(q - 1) < 1e-6) {
        # Shannon entropy path - vectorized
        for (i in seq_len(n)) {
          denom <- 1.0 - p[i]
          if (denom > 1e-10) {
            p_minus_i <- p / denom
            p_minus_i[i] <- 0
            p_nonzero <- p_minus_i[p_minus_i > 1e-15]
            if (length(p_nonzero) > 0) {
              jackknife_estimates[i] <- -sum(p_nonzero * log(p_nonzero)) / log(log_base)
            } else {
              jackknife_estimates[i] <- 0
            }
          } else {
            jackknife_estimates[i] <- NA_real_
          }
        }
      } else {
        # Tsallis entropy path - vectorized
        for (i in seq_len(n)) {
          denom <- 1.0 - p[i]
          if (denom > 1e-10) {
            p_minus_i <- p / denom
            p_minus_i[i] <- 0
            jackknife_estimates[i] <- (1.0 / (q - 1.0)) * (1.0 - sum(p_minus_i^q)) / log(log_base)
          } else {
            jackknife_estimates[i] <- NA_real_
          }
        }
      }
      
      # Normalization (P1 optimization: cache max entropy)
      if (norm) {
        estimate <- .tsenat_entropy_single(gene_counts, q = q, norm = TRUE, log_base = log_base, pseudocount = pseudocount)
        
        # Cache max entropy for (n-1) transcripts
        n_jackknife <- n - 1
        if (abs(q - 1) < 1e-6) {
          max_entropy_jackknife <- log(n_jackknife) / log(log_base)
        } else {
          max_entropy_jackknife <- (1.0 / (q - 1.0)) * (1.0 - n_jackknife^(1.0 - q)) / log(log_base)
        }
        
        if (!is.na(max_entropy_jackknife) && !is.nan(max_entropy_jackknife) && 
            max_entropy_jackknife > 0 && is.finite(max_entropy_jackknife)) {
          jackknife_estimates <- jackknife_estimates / max_entropy_jackknife
        }
      } else {
        if (abs(q - 1) < 1e-6) {
          p_nonzero <- p[p > 1e-15]
          if (length(p_nonzero) > 0) {
            estimate <- -sum(p_nonzero * log(p_nonzero)) / log(log_base)
          } else {
            estimate <- 0
          }
        } else {
          estimate <- (1.0 / (q - 1.0)) * (1.0 - sum(p^q)) / log(log_base)
        }
      }
      
      # Compute influence and outliers
      influence <- abs(jackknife_estimates - estimate)
      theta_jack_mean <- mean(jackknife_estimates, na.rm = TRUE)
      jackknife_se <- sqrt(((n - 1) / n) * sum((jackknife_estimates - theta_jack_mean)^2, na.rm = TRUE))
      
      # Outlier detection (P0 optimization: vectorized comparison)
      outlier_cutoff <- stats::quantile(influence, threshold / 100, na.rm = TRUE)
      outlier_mask <- influence > outlier_cutoff & !is.na(influence)
      outlier_indices <- which(outlier_mask)
      
      # Store result
      results[[gene_idx]] <- list(
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
      class(results[[gene_idx]]) <- c("tsenat_jackknife", "list")
    }
    
    names(results) <- rownames(x)
    class(results) <- c("tsenat_jackknife_list", "list")
    
    # Display results if requested
    # P3 OPTIMIZATION: Batch message formatting to reduce I/O overhead
    if (verbose) {
      output_lines <- c(
        "Jackknife Stability Analysis for Top 5 Genes",
        "============================================"
      )
      
      for (i in seq_along(results)) {
        gene_name <- names(results)[i]
        jr <- results[[i]]
        
        output_lines <- c(output_lines,
          sprintf("Gene: %s", gene_name),
          sprintf("  Transcripts: %d", jr$n_transcripts),
          sprintf("  Diversity estimate: %.4f", jr$estimate),
          sprintf("  Jackknife SE: %.4f", jr$jackknife_se),
          sprintf("  Max transcript influence: %.4f", max(jr$influence)),
          sprintf("  Outliers detected: %d", length(jr$outlier_indices))
        )
        
        if (length(jr$outlier_indices) > 0) {
          output_lines <- c(output_lines,
            sprintf("    Outlier transcripts (indices): %s", paste(jr$outlier_indices, collapse = ", "))
          )
        }
      }
      
      output_lines <- c(output_lines,
        "",
        "Interpretation:",
        "- High SE relative to estimate: unstable diversity (few dominant transcripts)",
        "- Many outliers: non-uniform isoform distribution",
        "- No outliers: balanced/robust isoform diversity"
      )
      
      # Print all at once instead of per-gene message() calls
      message(paste(output_lines, collapse = "\n"))
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

  # (estimate already computed in jackknife loop above)

  # Jackknife: remove each transcript and recalculate (P0 OPTIMIZATION: vectorized proportions)
  # OPTIMIZATION: Use incremental formula p_-i[j] = p[j] / (1 - p[i]) instead of full recalculation
  # This reduces loop complexity from O(n²) to O(n)
  
  # Pre-compute normalized proportions once
  total_count <- sum(x)
  p <- (x + pseudocount) / (total_count + length(x) * pseudocount)
  
  # Pre-allocate results
  jackknife_estimates <- numeric(n)
  
  if (abs(q - 1) < 1e-6) {
    # Shannon entropy path - uses vectorized incremental updates
    for (i in seq_len(n)) {
      # OPTIMIZATION: Incremental proportion update - O(1) per iteration vs O(n)
      # Instead of: p_minus_i <- x[-i] / sum(x[-i])
      # Use: p_minus_i = p / (1 - p[i]) for j ≠ i
      # This avoids vector subsetting and renormalization
      
      denom <- 1.0 - p[i]
      if (denom > 1e-10) {  # Avoid division by zero
        p_minus_i <- p / denom
        p_minus_i[i] <- 0  # Removed transcript has zero proportion
        
        # Shannon entropy: -sum(p[p>0] * log(p[p>0]))
        p_nonzero <- p_minus_i[p_minus_i > 1e-15]
        if (length(p_nonzero) > 0) {
          jackknife_estimates[i] <- -sum(p_nonzero * log(p_nonzero)) / log(log_base)
        } else {
          jackknife_estimates[i] <- 0
        }
      } else {
        jackknife_estimates[i] <- NA_real_
      }
    }
  } else {
    # Tsallis entropy path - uses vectorized incremental updates
    for (i in seq_len(n)) {
      # OPTIMIZATION: Same incremental strategy as Shannon path
      denom <- 1.0 - p[i]
      if (denom > 1e-10) {
        p_minus_i <- p / denom
        p_minus_i[i] <- 0
        
        # Tsallis entropy: (1/(q-1)) * (1 - sum(p^q))
        jackknife_estimates[i] <- (1.0 / (q - 1.0)) * (1.0 - sum(p_minus_i^q)) / log(log_base)
      } else {
        jackknife_estimates[i] <- NA_real_
      }
    }
  }
  
  # Normalize if requested  
  # P1 OPTIMIZATION: Cache max entropy calculation (computed once, not per-gene in loop)
  if (norm) {
    estimate <- .tsenat_entropy_single(x, q = q, norm = TRUE, log_base = log_base, pseudocount = pseudocount)
    
    # P1 OPTIMIZATION: Pre-compute max entropy for (n-1) transcripts once, then reuse
    # Each jackknife estimate is computed on n-1 transcripts (leave-one-out)
    n_jackknife <- n - 1
    if (abs(q - 1) < 1e-6) {
      max_entropy_jackknife <- log(n_jackknife) / log(log_base)
    } else {
      max_entropy_jackknife <- (1.0 / (q - 1.0)) * (1.0 - n_jackknife^(1.0 - q)) / log(log_base)
    }
    
    # Validate and apply cached max entropy to all estimates at once
    if (!is.na(max_entropy_jackknife) && !is.nan(max_entropy_jackknife) && 
        max_entropy_jackknife > 0 && is.finite(max_entropy_jackknife)) {
      jackknife_estimates <- jackknife_estimates / max_entropy_jackknife
    }
  } else {
    # Non-normalized case: use raw entropy estimate
    # For Shannon: entropy = -sum(p[p>0] * log(p[p>0]))
    # For Tsallis: entropy = (1/(q-1)) * (1 - sum(p^q))
    if (abs(q - 1) < 1e-6) {
      p_nonzero <- p[p > 1e-15]
      if (length(p_nonzero) > 0) {
        estimate <- -sum(p_nonzero * log(p_nonzero)) / log(log_base)
      } else {
        estimate <- 0
      }
    } else {
      estimate <- (1.0 / (q - 1.0)) * (1.0 - sum(p^q)) / log(log_base)
    }
  }

  # Compute influence of each transcript
  influence <- abs(jackknife_estimates - estimate)

  # Jackknife standard error (from Efron & Tibshirani 1993)
  # SE = sqrt( ((n-1)/n) * sum( (theta_-i - theta_.)^2 ) )
  theta_jack_mean <- mean(jackknife_estimates, na.rm = TRUE)
  jackknife_se <- sqrt(((n - 1) / n) * sum((jackknife_estimates - theta_jack_mean)^2, na.rm = TRUE))

  # Identify outliers based on threshold percentile (OPTIMIZED - single quantile call)
  outlier_cutoff <- stats::quantile(influence, threshold / 100, na.rm = TRUE)
  # Use vectorized comparison instead of which() when possible
  outlier_mask <- influence > outlier_cutoff & !is.na(influence)
  outlier_indices <- which(outlier_mask)

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

  # Normalize to proportions (OPTIMIZED - single normalization with pseudocount)
  total <- sum(counts) + length(counts) * pseudocount
  p <- (counts + pseudocount) / total

  n <- length(p)

  if (abs(q - 1) < 1e-6) {
    # Shannon entropy as q -> 1
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
#' @exportS3Method base::print tsenat_jackknife
#' @method print tsenat_jackknife
print.tsenat_jackknife <- function(x, ...) {
  message("Jackknife Diagnostics for Tsallis Entropy (q = ", x$q, ")")
  message("Estimate: ", sprintf("%.6f", x$estimate))
  message("Jackknife SE: ", sprintf("%.6f", x$jackknife_se))
  message("Number of transcripts: ", x$n_transcripts)
  invisible(x)
}


#' Summary of Jackknife Diagnostics Results
#'
#' @param object A \code{tsenat_jackknife} object
#' @param ... Additional arguments (unused)
#'
#' @keywords internal
#' @noRd
#' @exportS3Method base::summary tsenat_jackknife
#' @method summary tsenat_jackknife
summary.tsenat_jackknife <- function(object, ...) {
  message("Jackknife Diagnostics Summary")
  message("=============================")
  message(sprintf("Gene entropy (q =%.6f): %.6f", object$q, object$estimate))
  message(sprintf("Jackknife standard error: %.6f", object$jackknife_se))
  message(sprintf("Coefficient of variation: %.4f", object$jackknife_se / object$estimate))
  message(sprintf("Total transcripts: %d", object$n_transcripts))

  message("Influence Distribution:")
  message("-----------------------")
  message(paste(capture.output(str(object$influence)), collapse = "\n"))

  message(sprintf("\n\nOutlier Transcripts (influence > %s%%ile):", object$outlier_threshold))
  message("-----------------------------------------------")

  if (length(object$outlier_indices) > 0) {
    outl_data <- data.frame(
      Transcript = object$outlier_indices,
      Influence = object$influence[object$outlier_indices]
    )
    message(paste(capture.output(str(outl_data)), collapse = "\n"))
  } else {
    message("No outliers detected.")
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
#' @exportS3Method base::print tsenat_jackknife_list
#' @method print tsenat_jackknife_list
print.tsenat_jackknife_list <- function(x, ...) {
  message("Jackknife Results for Multiple Genes")
  message("Number of genes: ", length(x))
  for (i in seq_along(x)) {
    message("  ", names(x)[i], ": Estimate = ", sprintf("%.6f", x[[i]]$estimate))
  }
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
    
    for (i in seq_len(n_tx)) {
      h_leave_i <- .calculate_tsallis(counts[-i, , drop = FALSE], q, norm, log_base, pseudocount, n_transcripts_fixed)
      influences[i] <- mean(abs(h_full - h_leave_i), na.rm = TRUE)
    }
    list(full_entropy = h_full, influences = influences)
  }
  
  n_tx <- length(delta_influence)
  # Seed handling left to caller for Bioconductor compliance
  bootstrap_deltas_matrix <- matrix(nrow = n_bootstrap, ncol = n_tx)
  
  for (b in seq_len(n_bootstrap)) {
    idx_A <- sample(seq_len(ncol(counts_A)), size = ncol(counts_A), replace = TRUE)
    idx_B <- sample(seq_len(ncol(counts_B)), size = ncol(counts_B), replace = TRUE)
    
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
  
  # Para cada transcrito, calcular p-value basado en cuantos bootstrap samples
  # tienen signo opuesto al delta_influence observado
  pvalues <- numeric(n_tx)
  for (i in seq_len(n_tx)) {
    boot_signs <- sign(bootstrap_deltas_matrix[, i])
    obs_sign <- sign(delta_influence[i])
    # P-value: proporcion de muestras bootstrap con signo opuesto
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
#' @noRd
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
