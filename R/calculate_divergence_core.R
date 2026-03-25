# =========================================================================
# MAIN FUNCTION: Calculate Bootstrap Divergence
# =========================================================================

#' Calculate Bootstrap Divergence Confidence Intervals Across Genes
#'
#' **NEW ARCHITECTURE: Transcript-level counts -> Gene-level aggregation -> Tsallis divergence**
#' 
#' Computes bootstrap confidence intervals for Tsallis divergence comparing
#' two groups across multiple genes. Automatically aggregates transcript-level counts
#' to gene-level (per Paper I033: gene-level analysis for information-theoretic diversity).
#' Supports both sequential and parallel computation, with optional support for paired sample designs.
#' 
#' Returns a SummarizedExperiment object containing:
#' - **assay**: genes * q matrix of divergence estimates (one per q value)
#' - **rowData**: gene metadata including per-q divergence estimates, CIs, pattern classification
#' - **colData**: one row per q value with q-specific metadata
#' - **metadata**: processing parameters and summary statistics
#'
#' **INPUT & OUTPUT ARCHITECTURE:**
#' ```
#' calculate_divergence(se, res=NULL, ...)  
#'   Input:  SummarizedExperiment (raw TRANSCRIPT-level counts)
#'           Each row is a transcript; rowData must have gene_names/gene_name column
#'   Step 1: Auto-aggregates transcripts -> genes via colSums
#'   Step 2: Computes divergence for each gene across q values
#'   Output: SummarizedExperiment with:
#'           - assay: genes * q_values matrix (divergence estimates)
#'           - rowData: gene_name, per_q_pattern, estimate_q*, lower_ci_q*, etc.
#'           - colData: one row per q value
#'           - metadata: parameters, timing, sample sizes
#' ```
#' Matches `calculate_diversity()` input/output pattern: transcript counts SE -> gene-level derivative SE
#'
#' **DESIGN PRINCIPLE - Transcript-to-Gene Aggregation:**
#' Following Paper I033 ("Application of information theoretical approaches to assess diversity 
#' in single-cell transcriptomics"), divergence analysis operates on GENE-LEVEL expression profiles.
#' When input is transcript-level data (typical RNA-seq output), this function automatically:
#'   1. Identifies all transcripts for each gene (via rowData gene_names column)
#'   2. Sums counts across transcripts for each gene
#'   3. Computes divergence on aggregated gene-level counts
#' This ensures statistical validity (one observation per gene per sample) and biological relevance.
#'
#' @param se SummarizedExperiment object with transcript-level counts
#'           (assay called "counts", rowData with gene identifier columns)
#' @param group_col Character; colData column for group membership (optional).
#'            If NULL, auto-detects in this order: "group", "condition", "treatment", 
#'            "sample_type". If no match found, an error is raised.
#'            (default: NULL, auto-detect)
#' @param control_group Character; reference group name (optional).
#'            If NULL, auto-detects by: (1) looking for "Normal", "Control", "WT", etc.,
#'            or (2) selecting the group with fewer samples (typical case-control),
#'            or (3) first alphabetically.
#'            (default: NULL, auto-detect)
#' @param q Tsallis parameter (scalar or vector) (default: 1)
#' @param paired Logical; if TRUE or if paired_samples column detected, uses paired sample design.
#'               With bootstrap=TRUE, automatically detects paired samples from metadata
#'               column names (searched in order: "paired_samples", "pair_id", "pair_samples",
#'               "subject_id", "patient_id") and applies pair-respecting bootstrap resampling
#'               to preserve within-pair correlations (Papers C016, S102-S109).
#'               (default: FALSE)
#' @param bootstrap Logical; if TRUE, computes bootstrap confidence intervals (~2-3 sec/gene).
#'                  If FALSE, computes point estimates only (~0.02-0.05 sec/gene).
#'                  (default: FALSE)
#' @param nboot Number of bootstrap replicates (default: 1000)
#'               Note: ignored if bootstrap=FALSE
#' @param ci Confidence level (default: 0.95)
#' @param method Bootstrap method: "percentile" or "bca" (default: "percentile")
#' @param log_base Logarithm base (default: exp(1), natural log)
#' @param norm Logical or character; normalization/standardization mode (default: TRUE).
#'        Backward compatible: TRUE = "range", FALSE = "none".
#'        Options:
#'        - "none": Raw divergence values, no standardization
#'        - "range": Range standardization [0,1] per q (classic approach)
#'        - "zscore": Z-score standardization per q: (D_q - mean) / sd
#'          Useful for cross-study comparison; results in mean=0, sd=1
#'        - "log_odds_ratio": Log ratio relative to theoretical maximum
#'          D_norm = log(D_q / D_max) where D_max depends on q-value
#'          Interpretation: 0 = theoretical max, <0 = below max
#'        - "relative_reference": Ratio to reference group (requires group_col)
#'          Interpretation: Reference = 1, >1 higher than reference
#' @param pseudocount Pseudocount for stability (default: 0.5)
#' @param nthreads Number of CPU threads for parallel processing (default: 1).
#'                 Set to > 1 to parallelize gene-level bootstrap computations.
#'                 nthreads=NULL auto-detects available cores minus 1.
#' @param progress Logical; show progress bar and timing (default: TRUE)
#' @param seed Random seed (optional; NULL for non-reproducible)
#'
#' @return SummarizedExperiment object with:
#'   **assays** (genes * q matrices):
#'     - divergence: Divergence point estimates for each gene
#'   
#'   **rowData** (data frame with one row per gene):
#'     - gene_name: Gene identifier
#'     - estimate: Point estimate of divergence
#'     - lower_ci: Lower CI bound (NA if bootstrap=FALSE)
#'     - upper_ci: Upper CI bound (NA if bootstrap=FALSE)
#'     - ci_width: Width of confidence interval (NA if bootstrap=FALSE)
#'     - q: Tsallis parameter(s) used
#'     - nboot: Number of bootstrap replicates
#'     - method: Bootstrap method ("percentile", "bca", or NA)
#'     - computation_time_sec: Wall-clock time per gene (seconds)
#'     - error: Error message if computation failed, NA_character_ otherwise
#'   
#'   **colData**:
#'     - Inherited from input `se` (sample grouping, pairing info, etc.)
#'   
#'   **metadata** (list):
#'     - summary_stats: list with counts (total, successful, failed)
#'     - elapsed_time_sec: Total computation time
#'     - avg_time_per_gene: Average time per gene
#'     - genes_per_minute: Processing rate
#'     - bootstrap_config: list with bootstrap parameters (nboot, ci, method)
#'     - computation_mode: "sequential" or "parallel"
#'
#' @details
#' **Paired Sample Auto-Detection (NEW FEATURE):**
#' When bootstrap=TRUE, the function automatically detects paired sample metadata from colData:
#' - Searches for columns: "paired_samples", "pair_id", "pair_samples", "subject_id", "patient_id"
#' - If found, uses **pair-respecting bootstrap resampling**:
#'   * Resamples pair indices (not individual samples) with replacement
#'   * Preserves within-pair correlations critical for matched designs
#'   * Maintains statistical validity in paired experimental designs
#' - If paired=TRUE but no pairing detected, falls back to independent bootstrap with warning
#' 
#' **Scientific Justification (Papers validating auto-detection):**
#' - Papers C016: Bootstrap for confidence intervals requires preserving data structure
#' - Papers S102-S109: Paired design standards and statistical methods
#' - Papers I002-I004: Tsallis divergence mathematical foundation
#' 
#' **Computation Mode Selection:**
#' The function automatically selects between sequential and parallel processing:
#' - nthreads=1 (default): Direct sequential loop, minimal overhead
#' - nthreads > 1 & num_genes >= 5: Parallel PSOCK cluster
#' - nthreads > 1 & num_genes < 5: Falls back to sequential (overhead not warranted)
#'
#' **Performance Characteristics:**
#' - With bootstrap=TRUE (default):
#'   - Sequential: ~2-3 seconds per gene (nboot=1000, percentile method)
#'   - Parallel overhead: ~1-2 seconds initial cluster setup
#'   - Break-even point: ~10-20 genes
#' - With bootstrap=FALSE (point estimates only):
#'   - Sequential: ~0.02-0.05 seconds per gene (50-100* faster)
#'
#' **Gene Filtering:**
#' - Always process all genes in se
#'
#' **Database Verification (tsenat_papers.db):**
#' - Tsallis divergence mathematical foundation: Papers I001-I004 validate
#'   divergence formula and q-parameter effects
#' - Bootstrap methodology: Papers C016, S018, S030 validate percentile and BCa
#'   bootstrap for entropy/divergence estimates with confidence level >= 0.95
#' - Transcript aggregation: Paper C105 validates gene-level aggregation
#' - Divergence normalization: Papers C112, S196, S201 validate normalization
#'   approaches for effect size comparability (S197 - DESeq2 independent filtering)
#' @keywords internal
#' @noRd

calculate_divergence <- function(
    se,
    group_col = NULL,
    control_group = NULL,
    q = 1,
    paired = FALSE,
    bootstrap = FALSE,
    nboot = "auto",
    ci = 0.95,
    method = "percentile",
    norm = TRUE,
    log_base = exp(1),
    pseudocount = 0.5,
    nthreads = 1,
    progress = FALSE,
    verbose = TRUE,
    seed = NULL) {

  # =========================================================================
  # INPUT VALIDATION
  # =========================================================================

  # Normalize norm parameter: coerce logical to character for backward compatibility
  norm <- .normalize_norm_parameter(norm)
  
  # Validate and sort q parameter
  q <- .validate_and_sort_q_values(q)
  
  # Validate SE input
  .validate_se_input(se)

  # =========================================================================
  # AUTO-DETECT GROUP COLUMN AND CONTROL GROUP
  # =========================================================================
  
  if (is.null(group_col) || is.null(control_group)) {
    auto_groups <- .auto_detect_groups(se)
    
    if (is.null(group_col)) {
      if (is.na(auto_groups$group_col)) {
        stop("Could not auto-detect group column in colData. ",
             "Available columns: ", 
             paste(colnames(SummarizedExperiment::colData(se)), collapse = ", "),
             ". Please specify 'group_col' explicitly.",
             call. = FALSE)
      }
      group_col <- auto_groups$group_col
      if (progress) {
        message("[calculate_divergence] Auto-detected group_col='", group_col, "'")
      }
    }
    
    if (is.null(control_group)) {
      if (is.na(auto_groups$control_group)) {
        stop("Could not auto-detect control_group. Found groups: ",
             paste(auto_groups$groups, collapse = ", "),
             ". Please specify 'control_group' explicitly.",
             call. = FALSE)
      }
      control_group <- auto_groups$control_group
      if (progress) {
        message("[calculate_divergence] Auto-detected control_group='", control_group, "'")
      }
    }
  }

  # Get list of UNIQUE genes to process (not transcripts)
  rd <- SummarizedExperiment::rowData(se)
  
  # Identify gene name column
  gene_col <- .identify_gene_column(se)
  
  # Extract unique gene list
  all_gene_names <- .extract_gene_list(se, gene_col)

  # Process unique genes only
  genes_to_process <- all_gene_names
  gene_indices <- seq_along(all_gene_names)

  num_genes <- length(gene_indices)
  
  if (num_genes == 0) {
    stop("No genes to process. ",
         "se gene names (first 3): ", paste(head(all_gene_names, 3), collapse=", "))
  }

  # =========================================================================
  # BOOTSTRAP CONFIGURATION
  # =========================================================================

  if (!is.logical(bootstrap) || length(bootstrap) != 1) {
    stop("bootstrap must be a logical (TRUE/FALSE)")
  }

  if (!bootstrap) {
    nboot <- 0
  }

  # AUTO-SELECT NBOOT WHEN "auto"
  if (bootstrap && identical(nboot, "auto")) {
    num_genes <- nrow(se)
    use_bca <- method == "bca"
    nboot <- suggest_nboot(num_genes, use_bca = use_bca, nthreads = nthreads)
    if (progress) {
      message("Auto-selected nboot =", nboot, "for", num_genes, "genes")
    }
  }

  # Configure parallel execution
  parallel_config <- .configure_parallel_execution(nthreads, num_genes)
  nthreads <- parallel_config$nthreads
  use_parallel <- parallel_config$use_parallel

  # =========================================================================
  # PAIRED SAMPLE DETECTION (NEW FEATURE)
  # =========================================================================
  
  pair_ids <- NULL
  pairing_info <- ""
  
  if (bootstrap) {
    # Auto-detect paired samples when bootstrap=TRUE (or when paired=TRUE)
    pair_detected <- .detect_pair_ids(se)
    
    if (pair_detected$num_pairs > 0) {
      pair_ids <- pair_detected$pair_ids
      pairing_info <- sprintf(" [paired: %d unique pairs from '%s' column]", 
                              pair_detected$num_pairs, 
                              pair_detected$column_name)
      
      if (paired == FALSE && progress) {
        message("NOTE: Paired sample structure detected in '", 
            pair_detected$column_name, "' column.\n",
            "      Using pair-respecting bootstrap resampling.")
      }
    } else {
      if (paired == TRUE && progress) {
        message("paired=TRUE but no pair ID column detected in colData.",
            " Using independent bootstrap resampling instead.")
      }
    }
  }

  if (progress) {
    mode_desc <- if (bootstrap) {
      paste0("bootstrap with ", nboot, " replicates (", method, ")", pairing_info)
    } else {
      "point estimates only"
    }

    mode_str <- if (use_parallel) "Parallel" else "Sequential"
    thread_desc <- if (use_parallel) paste0(" on ", nthreads, " threads") else ""
    message(mode_str, " mode: ", num_genes, " genes", thread_desc, " [", mode_desc, "]")
  }

  # =========================================================================
  # COMPUTATION (SEQUENTIAL OR PARALLEL)
  # =========================================================================

  start_time <- Sys.time()
  results_list <- list()

  if (!use_parallel) {
    # ====== SEQUENTIAL PROCESSING ======
    for (i in seq_along(gene_indices)) {
      gene_idx <- gene_indices[i]
      
      # Use consolidated gene processing helper
      results_list[[i]] <- .process_single_gene(
        gene_idx, all_gene_names, se, gene_col, rd,
        group_col, control_group, q, nboot, ci, method,
        log_base, pseudocount, seed, pair_ids
      )

      if (progress && (i %% 10 == 0)) {
        elapsed <- as.numeric(Sys.time() - start_time, units = "secs")
        rate <- (i / elapsed) * 60
        message("[", i, "/", num_genes, "] (", sprintf("%.1f genes/min", rate), ")")
      }
    }

  } else {
    # ====== PARALLEL PROCESSING ======

    cl <- parallel::makeCluster(nthreads, type = "PSOCK")
    on.exit(parallel::stopCluster(cl), add = TRUE)

    # Export helper functions and dependencies to cluster
    parallel::clusterExport(cl, 
      c(".process_single_gene", "calculate_divergence_bootstrap", ".tsallis_divergence_scalar",
        ".aggregate_counts_for_gene", ".extract_group_counts", ".compute_divergence_per_q",
        ".build_bootstrap_args", ".make_error_result",
        "rd", "gene_col", "all_gene_names", "se", "q", "nboot", "ci", "method",
        "log_base", "pseudocount", "seed", "group_col", "control_group", "pair_ids"),
      envir = environment())

    parallel::clusterCall(cl, function() {
      requireNamespace("SummarizedExperiment", quietly = TRUE)
    })

    results_list <- parallel::parLapply(cl, seq_along(gene_indices), function(i) {
      gene_idx <- gene_indices[i]
      
      # Use consolidated gene processing helper
      .process_single_gene(
        gene_idx, all_gene_names, se, gene_col, rd,
        group_col, control_group, q, nboot, ci, method,
        log_base, pseudocount, seed, pair_ids
      )
    })
  }

  elapsed <- as.numeric(Sys.time() - start_time, units = "secs")

  # =========================================================================
  # RESULTS COMPILATION - Build genes * q matrix
  # =========================================================================

  # Initialize result matrices
  matrices <- .initialize_result_matrices(num_genes, q)
  assay_matrix <- matrices$assay
  row_data_df <- matrices$rowData
  
  # Populate matrices from results_list
  populated <- .populate_result_matrices(results_list, assay_matrix, row_data_df, q)
  assay_matrix <- populated$assay
  row_data_df <- populated$rowData
  
  # Set row names
  rownames(row_data_df) <- row_data_df$gene_name
  rownames(assay_matrix) <- row_data_df$gene_name
  
  # Classify per-q patterns: RARE_DRIVEN, ABUNDANT_DRIVEN, or BALANCED
  row_data_df$per_q_pattern <- NA_character_
  if (length(q) > 1) {
    for (i in seq_len(nrow(row_data_df))) {
      # Get divergence values across q for this gene
      per_q_divs <- assay_matrix[i, ]
      names(per_q_divs) <- paste0("q_", q)
      
      # Only classify if we have valid values
      if (sum(!is.na(per_q_divs)) >= 2) {
        pattern <- classify_q_pattern(per_q_divs)
        row_data_df$per_q_pattern[i] <- if (is.na(pattern)) "UNCLASSIFIED" else pattern
      }
    }
  }
  
  # Populate generic estimate/lower_ci/upper_ci columns using reference q value (q=1)
  # This is needed for downstream functions like effect_sizes_divergence()
  q_ref <- 1.0
  q_idx <- which.min(abs(q - q_ref))  # Find closest q to 1.0
  if (length(q_idx) > 0 && q_idx <= length(q)) {
    ref_q <- q[q_idx]
    estimate_col <- paste0("estimate_q", ref_q)
    lower_ci_col <- paste0("lower_ci_q", ref_q)
    upper_ci_col <- paste0("upper_ci_q", ref_q)
    ci_width_col <- paste0("ci_width_q", ref_q)
    
    # Copy per-q columns to generic columns for downstream compatibility
    if (estimate_col %in% colnames(row_data_df)) {
      row_data_df$estimate <- row_data_df[[estimate_col]]
      row_data_df$lower_ci <- row_data_df[[lower_ci_col]]
      row_data_df$upper_ci <- row_data_df[[upper_ci_col]]
      row_data_df$ci_width <- row_data_df[[ci_width_col]]
    }
  }
  
  # Apply normalization if requested (BEFORE creating SE)
  if (norm != "none") {
    if (progress) {
      message(sprintf("Applying '%s' normalization to divergence estimates...", norm))
    }
    
    normalized <- .apply_divergence_normalization(
      assay_matrix = assay_matrix,
      row_data_df = row_data_df,
      q_vals = q,
      norm = norm
    )
    assay_matrix <- normalized$assay
    row_data_df <- normalized$rowData
  }

  # Summary statistics
  num_success <- sum(is.na(row_data_df$error))
  num_errors <- num_genes - num_success

  if (progress) {
    message("\nDIVERGENCE COMPUTATION COMPLETE")
    message("Summary:")
    message("  Genes processed:        ", num_genes)
    message("  Successful:             ", num_success)
    message("  Failed:                 ", num_errors)
    message("  Total elapsed time:     ", sprintf("%.1f seconds", elapsed))
    message("  Average per gene:       ", sprintf("%.2f seconds", elapsed / num_genes))
    message("  Genes per minute:       ", sprintf("%.1f", (num_genes / elapsed) * 60))

    if (num_errors > 0) {
      message("Failed genes:")
      failed <- row_data_df[!is.na(row_data_df$error), ]
      for (i in seq_len(min(10, nrow(failed)))) {
        message(sprintf("  [%d] %s: %s", i, failed$gene_name[i], failed$error[i]))
      }
      if (num_errors > 10) {
        message("  ... and", num_errors - 10, "more")
      }
    }
  }

  # =========================================================================
  # CREATE SUMMARIZED EXPERIMENT
  # =========================================================================
  
  result_se <- .construct_result_se(
    assay_matrix = assay_matrix,
    row_data_df = row_data_df,
    q_vals = q,
    elapsed = elapsed,
    nboot = nboot,
    ci = ci,
    method = method,
    norm = norm,
    use_parallel = use_parallel,
    num_genes = num_genes,
    num_errors = num_errors
  )

  return(result_se)
}

# ============================================================================
# HELPER FUNCTIONS FOR DIVERGENCE CALCULATION
# Extracted to meet Bioconductor ≤50 line requirement (March 2026)
# ============================================================================

# INPUT VALIDATION & CONFIGURATION HELPERS
# ============================================================================

#' Normalize norm parameter for backward compatibility
#' Coerces logical values to character strings
#' @noRd
.normalize_norm_parameter <- function(norm) {
    if (is.logical(norm)) {
        norm <- if (norm) "range" else "none"
    }
    match.arg(norm, choices = c("none", "range", "zscore", 
                                 "log_odds_ratio", "relative_reference"))
}

#' Validate and sort q-parameter values
#' Ensures q >= 0 and returns sorted vector
#' @noRd
.validate_and_sort_q_values <- function(q) {
    q <- sort(as.numeric(q))
    if (any(q < 0)) {
        stop("q parameter must be >= 0. ",
             "Note: q should be in range [0, 3] for typical use. ",
             "q=0 represents uniform divergence. ",
             "Got: ", paste(q, collapse = ", "))
    }
    q
}

#' Validate SummarizedExperiment input
#' @noRd
.validate_se_input <- function(se) {
    if (!methods::is(se, "SummarizedExperiment")) {
        stop("se must be a SummarizedExperiment object")
    }
    TRUE
}

#' Identify gene name/ID column in rowData
#' Preference: gene_name (human-readable) > gene_id (ensembl)
#' @noRd
.identify_gene_column <- function(se) {
    rd <- SummarizedExperiment::rowData(se)
    gene_col_candidates <- c("gene_name", "gene_id")
    
    if (is.null(rd)) return(NA_character_)
    
    for (col in gene_col_candidates) {
        if (col %in% colnames(rd)) {
            return(col)
        }
    }
    NA_character_
}

#' Extract unique gene list from SummarizedExperiment
#' @noRd
.extract_gene_list <- function(se, gene_col) {
    rd <- SummarizedExperiment::rowData(se)
    
    if (!is.na(gene_col) && !is.null(rd)) {
        all_genes <- unique(as.character(rd[[gene_col]]))
    } else {
        all_genes <- rownames(se)
    }
    
    if (length(all_genes) == 0) {
        stop("se must have gene identifiers in rowData or rownames")
    }
    
    all_genes
}

#' Configure parallel execution parameters
#' Auto-detects cores and decides between sequential/parallel
#' @noRd
.configure_parallel_execution <- function(nthreads, num_genes) {
    if (is.null(nthreads)) {
        nthreads <- parallel::detectCores() - 1
        nthreads <- max(1, nthreads)
    }
    
    if (!is.numeric(nthreads) || nthreads < 1) {
        stop("nthreads must be a positive integer")
    }
    nthreads <- as.integer(nthreads)
    
    use_parallel <- (num_genes >= 5 && nthreads > 1)
    
    list(nthreads = nthreads, use_parallel = use_parallel)
}

# GENE PROCESSING HELPERS
# ============================================================================

#' Aggregate transcript-level counts to gene-level
#' Sums counts across all transcripts for a given gene
#' @noRd
.aggregate_counts_for_gene <- function(se, target_gene, gene_col, rd) {
    # Find ALL transcripts for this gene
    if (!is.na(gene_col) && !is.null(rd)) {
        gene_transcript_indices <- which(as.character(rd[[gene_col]]) == target_gene)
    } else {
        gene_transcript_indices <- which(rownames(se) == target_gene)
    }
    
    if (length(gene_transcript_indices) == 0) {
        return(NULL)
    }
    
    # Aggregate counts across all transcripts for this gene
    counts_matrix <- as.matrix(SummarizedExperiment::assay(se, "counts")[gene_transcript_indices, , drop = FALSE])
    colSums(counts_matrix)
}

#' Extract group-specific counts
#' Splits gene-level counts by group membership
#' @noRd
.extract_group_counts <- function(counts_gene, groups, control_group) {
    if (length(counts_gene) != length(groups)) {
        stop("Length mismatch: counts_gene and groups must have same length")
    }
    
    x <- counts_gene[groups == control_group]
    y <- counts_gene[groups != control_group]
    
    list(control = x, treatment = y)
}

#' Construct error result structure
#' Standardized format for failed gene computations
#' @noRd
.make_error_result <- function(gene_name, q_vals, error_msg, elapsed_sec = NA_real_) {
    list(
        gene_name = gene_name,
        results_per_q = rep(list(list(
            estimate = NA_real_, lower_ci = NA_real_, 
            upper_ci = NA_real_, method = NA_character_
        )), length(q_vals)),
        computation_time_sec = elapsed_sec,
        error = error_msg
    )
}

#' Build bootstrap arguments for calculate_divergence_bootstrap
#' Conditionally includes pair_ids if detected
#' @noRd
.build_bootstrap_args <- function(x, y, q_val, nboot, ci, method, 
                                   log_base, pseudocount, gene_name, 
                                   seed, pair_ids = NULL) {
    args <- list(
        x = x, y = y, q = q_val, nboot = nboot, ci = ci, method = method,
        log_base = log_base, pseudocount = pseudocount,
        gene_name = gene_name, print_results = FALSE, seed = seed,
        paired = !is.null(pair_ids)
    )
    
    if (!is.null(pair_ids)) {
        args$pair_ids <- pair_ids
    }
    
    args
}

#' Compute divergence for all q values for a single gene
#' Returns list of results, one per q value
#' @noRd
.compute_divergence_per_q <- function(x, y, q_vals, nboot, ci, method,
                                       log_base, pseudocount, gene_name,
                                       seed, pair_ids = NULL) {
    gene_results <- list()
    
    for (j in seq_along(q_vals)) {
        q_val <- q_vals[j]
        
        bootstrap_args <- .build_bootstrap_args(
            x, y, q_val, nboot, ci, method, 
            log_base, pseudocount, gene_name, 
            seed, pair_ids
        )
        
        result <- do.call(calculate_divergence_bootstrap, bootstrap_args)
        gene_results[[j]] <- result
    }
    
    gene_results
}

# RESULTS COMPILATION HELPERS
# ============================================================================

#' Initialize result matrices for results compilation
#' Creates assay matrix and rowData structure
#' @noRd
.initialize_result_matrices <- function(num_genes, q_vals) {
    num_q_vals <- length(q_vals)
    
    assay_matrix <- matrix(NA_real_, nrow = num_genes, ncol = num_q_vals,
                           dimnames = list(NULL, paste0("q_", q_vals)))
    
    row_data_df <- data.frame(
        gene_name = character(num_genes),
        error = character(num_genes),
        computation_time_sec = numeric(num_genes),
        stringsAsFactors = FALSE
    )
    
    # Add columns for each q value's metadata
    for (j in seq_len(num_q_vals)) {
        row_data_df[[paste0("estimate_q", q_vals[j])]] <- NA_real_
        row_data_df[[paste0("lower_ci_q", q_vals[j])]] <- NA_real_
        row_data_df[[paste0("upper_ci_q", q_vals[j])]] <- NA_real_
        row_data_df[[paste0("ci_width_q", q_vals[j])]] <- NA_real_
        row_data_df[[paste0("method_q", q_vals[j])]] <- NA_character_
        row_data_df[[paste0("nboot_q", q_vals[j])]] <- NA_integer_
    }
    
    list(assay = assay_matrix, rowData = row_data_df)
}

#' Populate result matrices from results_list
#' Extracts individual gene results and fills matrices
#' @noRd
.populate_result_matrices <- function(results_list, assay_matrix, row_data_df, q_vals) {
    num_q_vals <- length(q_vals)
    
    for (i in seq_along(results_list)) {
        result <- results_list[[i]]
        
        row_data_df$gene_name[i] <- result$gene_name
        row_data_df$error[i] <- if (is.na(result$error)) NA_character_ else result$error
        row_data_df$computation_time_sec[i] <- result$computation_time_sec
        
        if (is.na(result$error)) {
            for (j in seq_len(num_q_vals)) {
                q_res <- result$results_per_q[[j]]
                assay_matrix[i, j] <- q_res$estimate
                
                row_data_df[[paste0("estimate_q", q_vals[j])]][i] <- q_res$estimate
                row_data_df[[paste0("lower_ci_q", q_vals[j])]][i] <- q_res$lower_ci
                row_data_df[[paste0("upper_ci_q", q_vals[j])]][i] <- q_res$upper_ci
                row_data_df[[paste0("method_q", q_vals[j])]][i] <- q_res$method %||% NA_character_
                row_data_df[[paste0("nboot_q", q_vals[j])]][i] <- as.integer(q_res$nboot %||% 0)
                
                # Compute CI width
                if (!is.na(q_res$lower_ci) && !is.na(q_res$upper_ci)) {
                    row_data_df[[paste0("ci_width_q", q_vals[j])]][i] <- 
                        q_res$upper_ci - q_res$lower_ci
                }
            }
        }
    }
    
    list(assay = assay_matrix, rowData = row_data_df)
}

# NORMALIZATION HELPERS
# ============================================================================

#' Apply range normalization [0,1]
#' @noRd
.normalize_range <- function(assay_matrix, row_data_df, q_vals) {
    for (j in seq_len(ncol(assay_matrix))) {
        col_vals <- assay_matrix[, j]
        valid_vals <- col_vals[!is.na(col_vals)]
        
        if (length(valid_vals) > 1) {
            min_val <- min(valid_vals)
            max_val <- max(valid_vals)
            range_val <- max_val - min_val
            
            if (range_val > 0) {
                assay_matrix[, j] <- (col_vals - min_val) / range_val
                
                # Apply same to estimate and CI bounds
                estimate_col <- paste0("estimate_q", q_vals[j])
                lower_col <- paste0("lower_ci_q", q_vals[j])
                upper_col <- paste0("upper_ci_q", q_vals[j])
                
                if (estimate_col %in% colnames(row_data_df)) {
                    row_data_df[[estimate_col]] <- (row_data_df[[estimate_col]] - min_val) / range_val
                    row_data_df[[lower_col]] <- (row_data_df[[lower_col]] - min_val) / range_val
                    row_data_df[[upper_col]] <- (row_data_df[[upper_col]] - min_val) / range_val
                }
            }
        }
    }
    
    list(assay = assay_matrix, rowData = row_data_df)
}

#' Apply z-score normalization
#' @noRd
.normalize_zscore <- function(assay_matrix, row_data_df, q_vals) {
    for (j in seq_len(ncol(assay_matrix))) {
        col_vals <- assay_matrix[, j]
        valid_vals <- col_vals[!is.na(col_vals)]
        
        if (length(valid_vals) > 1) {
            mean_val <- mean(valid_vals)
            sd_val <- sd(valid_vals)
            
            if (sd_val > 0) {
                assay_matrix[, j] <- (col_vals - mean_val) / sd_val
                
                # Apply same to estimate and CI bounds
                estimate_col <- paste0("estimate_q", q_vals[j])
                lower_col <- paste0("lower_ci_q", q_vals[j])
                upper_col <- paste0("upper_ci_q", q_vals[j])
                
                if (estimate_col %in% colnames(row_data_df)) {
                    row_data_df[[estimate_col]] <- (row_data_df[[estimate_col]] - mean_val) / sd_val
                    row_data_df[[lower_col]] <- (row_data_df[[lower_col]] - mean_val) / sd_val
                    row_data_df[[upper_col]] <- (row_data_df[[upper_col]] - mean_val) / sd_val
                }
            }
        }
    }
    
    list(assay = assay_matrix, rowData = row_data_df)
}

#' Apply log odds ratio normalization
#' D_norm = log(D_q / D_max) where D_max depends on q
#' @noRd
.normalize_log_odds_ratio <- function(assay_matrix, row_data_df, q_vals) {
    for (j in seq_len(ncol(assay_matrix))) {
        q_val <- q_vals[j]
        col_vals <- assay_matrix[, j]
        valid_vals <- col_vals[!is.na(col_vals)]
        
        # Maximum divergence depends on q
        d_max <- if (q_val < 1.0) log(2) else 1.0
        
        if (d_max > 0) {
            assay_matrix[, j] <- log(pmax(col_vals, 1e-10) / d_max)
            
            # Apply same to estimate and CI bounds
            estimate_col <- paste0("estimate_q", q_vals[j])
            lower_col <- paste0("lower_ci_q", q_vals[j])
            upper_col <- paste0("upper_ci_q", q_vals[j])
            
            if (estimate_col %in% colnames(row_data_df)) {
                row_data_df[[estimate_col]] <- log(pmax(row_data_df[[estimate_col]], 1e-10) / d_max)
                row_data_df[[lower_col]] <- log(pmax(row_data_df[[lower_col]], 1e-10) / d_max)
                row_data_df[[upper_col]] <- log(pmax(row_data_df[[upper_col]], 1e-10) / d_max)
            }
        }
    }
    
    list(assay = assay_matrix, rowData = row_data_df)
}

#' Apply relative reference normalization
#' Ratio to reference group mean per q value
#' @noRd
.normalize_relative_reference <- function(assay_matrix, row_data_df, q_vals) {
    for (j in seq_len(ncol(assay_matrix))) {
        col_vals <- assay_matrix[, j]
        valid_vals <- col_vals[!is.na(col_vals)]
        
        if (length(valid_vals) > 0) {
            reference_mean <- mean(valid_vals, na.rm = TRUE)
            
            if (reference_mean > 0) {
                assay_matrix[, j] <- col_vals / reference_mean
                
                # Apply same to estimate and CI bounds
                estimate_col <- paste0("estimate_q", q_vals[j])
                lower_col <- paste0("lower_ci_q", q_vals[j])
                upper_col <- paste0("upper_ci_q", q_vals[j])
                
                if (estimate_col %in% colnames(row_data_df)) {
                    row_data_df[[estimate_col]] <- row_data_df[[estimate_col]] / reference_mean
                    row_data_df[[lower_col]] <- row_data_df[[lower_col]] / reference_mean
                    row_data_df[[upper_col]] <- row_data_df[[upper_col]] / reference_mean
                }
            }
        }
    }
    
    list(assay = assay_matrix, rowData = row_data_df)
}

#' Apply normalization to results
#' Dispatcher that calls appropriate normalization method
#' @noRd
.apply_divergence_normalization <- function(assay_matrix, row_data_df, q_vals, norm) {
    if (norm == "none") {
        return(list(assay = assay_matrix, rowData = row_data_df))
    }
    
    if (norm == "range") {
        return(.normalize_range(assay_matrix, row_data_df, q_vals))
    }
    
    if (norm == "zscore") {
        return(.normalize_zscore(assay_matrix, row_data_df, q_vals))
    }
    
    if (norm == "log_odds_ratio") {
        return(.normalize_log_odds_ratio(assay_matrix, row_data_df, q_vals))
    }
    
    if (norm == "relative_reference") {
        return(.normalize_relative_reference(assay_matrix, row_data_df, q_vals))
    }
    
    list(assay = assay_matrix, rowData = row_data_df)
}

#' Generate computation summary statistics
#' Formats timing information for progress messages
#' @noRd
.generate_computation_summary <- function(elapsed, num_genes, num_errors, row_data_df) {
    num_success <- num_genes - num_errors
    
    list(
        total_elapsed = elapsed,
        avg_per_gene = elapsed / num_genes,
        genes_per_minute = (num_genes / elapsed) * 60,
        successful = num_success,
        failed = num_errors,
        failed_details = if (num_errors > 0) {
            failed <- row_data_df[!is.na(row_data_df$error), ]
            failed[seq_len(min(10, nrow(failed))), c("gene_name", "error")]
        } else NULL
    )
}

#' Construct final SummarizedExperiment output
#' Builds SE with assays, rowData, colData, and metadata
#' @noRd
.construct_result_se <- function(assay_matrix, row_data_df, q_vals, 
                                  elapsed, nboot, ci, method, norm, 
                                  use_parallel, num_genes, num_errors) {
    assays_list <- list(divergence = assay_matrix)
    
    col_data_output <- data.frame(
        q_value = q_vals,
        sample_type = rep("divergence_estimate", length(q_vals)),
        row.names = paste0("q_", q_vals)
    )
    
    num_success <- num_genes - num_errors
    
    SummarizedExperiment::SummarizedExperiment(
        assays = assays_list,
        rowData = row_data_df,
        colData = col_data_output,
        metadata = list(
            summary_stats = list(
                total_genes = num_genes,
                successful = num_success,
                failed = num_errors
            ),
            elapsed_time_sec = elapsed,
            avg_time_per_gene = elapsed / num_genes,
            genes_per_minute = (num_genes / elapsed) * 60,
            bootstrap_config = list(
                nboot = nboot,
                ci = ci,
                method = method
            ),
            normalization = norm,
            computation_mode = if (use_parallel) "parallel" else "sequential"
        )
    )
}

#' Process a single gene for divergence computation
#' Consolidates logic shared between sequential and parallel processing
#' @noRd
.process_single_gene <- function(gene_idx, all_gene_names, se, gene_col, rd, 
                                 group_col, control_group, q, nboot, ci, method,
                                 log_base, pseudocount, seed, pair_ids) {
    target_gene <- all_gene_names[gene_idx]
    gene_name <- target_gene
    gene_start <- Sys.time()
    
    tryCatch({
        # Get gene-level counts via transcript aggregation
        counts_gene <- .aggregate_counts_for_gene(se, target_gene, gene_col, rd)
        
        if (is.null(counts_gene)) {
            return(.make_error_result(gene_name, q, "No transcripts found for gene"))
        }
        
        # Extract group-specific counts
        groups <- se[[group_col]]
        group_counts <- .extract_group_counts(counts_gene, groups, control_group)
        x <- group_counts$control
        y <- group_counts$treatment
        
        if (length(x) == 0 || length(y) == 0) {
            return(.make_error_result(gene_name, q, "Insufficient group samples"))
        }
        
        # Compute divergence for each q value
        gene_results <- .compute_divergence_per_q(x, y, q, nboot, ci, method,
                                                   log_base, pseudocount, gene_name,
                                                   seed, pair_ids)
        
        gene_elapsed <- as.numeric(Sys.time() - gene_start, units = "secs")
        
        list(
            gene_name = gene_name,
            results_per_q = gene_results,
            computation_time_sec = gene_elapsed,
            error = NA_character_
        )
    }, error = function(e) {
        .make_error_result(gene_name, q, as.character(e$message),
                          as.numeric(Sys.time() - gene_start, units = "secs"))
    })
}
