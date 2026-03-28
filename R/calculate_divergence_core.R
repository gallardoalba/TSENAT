# =========================================================================
# NEW ORCHESTRATOR HELPERS (March 2026 refactoring)
# =========================================================================

#' Validate group columns and auto-detect missing ones
#' Consolidates duplicated if-pattern for group_col and control_group

#' @noRd
.tsenat_validate_and_auto_detect_groups <- function(se, group_col, control_group, progress) {
  if (is.null(group_col) || is.null(control_group)) {
    auto_groups <- .tsenat_auto_detect_groups(se)
    
    if (is.null(group_col)) {
      if (is.na(auto_groups$group_col)) {
        stop("Could not auto-detect group column in colData. ",
             "Available columns: ", 
             paste(colnames(SummarizedExperiment::colData(se)), collapse = ", "),
             ". Please specify 'group_col' explicitly.",
             call. = FALSE)
      }
      group_col <- auto_groups$group_col
      if (progress) message("[calculate_divergence] Auto-detected group_col='", group_col, "'")
    }
    
    if (is.null(control_group)) {
      if (is.na(auto_groups$control_group)) {
        stop("Could not auto-detect control_group. Found groups: ",
             paste(auto_groups$groups, collapse = ", "),
             ". Please specify 'control_group' explicitly.",
             call. = FALSE)
      }
      control_group <- auto_groups$control_group
      if (progress) message("[calculate_divergence] Auto-detected control_group='", control_group, "'")
    }
  }
  
  list(group_col = group_col, control_group = control_group)
}

#' Prepare genes for processing (identification + extraction)
#' Consolidates gene column identification and unique gene extraction

#' @noRd
.tsenat_prepare_genes_processing <- function(se) {
  rd <- SummarizedExperiment::rowData(se)
  gene_col <- .tsenat_identify_gene_column(se)
  all_gene_names <- .tsenat_extract_gene_list(se, gene_col)
  gene_indices <- seq_along(all_gene_names)
  
  if (length(gene_indices) == 0) {
    stop("No genes to process. ",
         "se gene names (first 3): ", paste(head(all_gene_names, 3), collapse=", "))
  }
  
  list(
    gene_col = gene_col,
    all_gene_names = all_gene_names,
    gene_indices = gene_indices,
    rd = rd,
    num_genes = length(gene_indices)
  )
}

#' Configure bootstrap and parallel execution parameters
#' Fixes nboot bug and consolidates configuration logic

#' @noRd
.tsenat_bootstrap_configure_parallel <- function(bootstrap, nboot, method, num_genes, nthreads, progress) {
  # Validate bootstrap flag - use isTRUE to safely handle NA
  if (!isTRUE(bootstrap) && !isFALSE(bootstrap)) {
    bootstrap <- FALSE  # Default to no bootstrap if invalid
  }
  
  if (!isTRUE(bootstrap)) {
    nboot <- 0
  }
  
  # AUTO-SELECT NBOOT WHEN "auto"
  if (isTRUE(bootstrap) && identical(nboot, "auto")) {
    use_bca <- !is.null(method) && identical(method, "bca")
    nboot <- .suggest_nboot(num_genes, use_bca = use_bca, nthreads = nthreads)
    if (isTRUE(progress)) {
      message("Auto-selected nboot =", nboot, "for", num_genes, "genes")
    }
  }
  
  # Configure parallel execution (fixes line 257 bug by using num_genes parameter)
  parallel_config <- .tsenat_configure_parallel(nthreads, num_genes)
  
  list(
    nboot = nboot,
    nthreads = parallel_config$nthreads,
    use_parallel = parallel_config$use_parallel
  )
}

#' Prepare paired sample and progress information
#' Consolidates paired detection and progress message assembly

#' @noRd
.tsenat_prepare_divergence_execution <- function(se, bootstrap, paired, nboot, method, nthreads, progress) {
  pair_ids <- NULL
  pairing_info <- ""
  
  if (isTRUE(bootstrap)) {  # Use isTRUE to safely handle NA
    # Auto-detect paired samples
    pair_detected <- .tsenat_detect_pair_ids(se)
    
    if (pair_detected$num_pairs > 0) {
      pair_ids <- pair_detected$pair_ids
      pairing_info <- sprintf(" [paired: %d unique pairs from '%s' column]", 
                              pair_detected$num_pairs, 
                              pair_detected$column_name)
      
      if (isFALSE(paired) && progress) {  # Use isFALSE to safely handle NA
        message("NOTE: Paired sample structure detected in '", 
            pair_detected$column_name, "' column.\n",
            "      Using pair-respecting bootstrap resampling.")
      }
    } else {
      if (isTRUE(paired) && progress) {  # Use isTRUE to safely handle NA
        message("paired=TRUE but no pair ID column detected in colData.",
            " Using independent bootstrap resampling instead.")
      }
    }
  }
  
  if (progress) {
    mode_desc <- if (isTRUE(bootstrap)) {  # Use isTRUE to safely handle NA
      paste0("bootstrap with ", nboot, " replicates (", method, ")", pairing_info)
    } else {
      "point estimates only"
    }
    
    mode_str <- if (!is.null(nthreads) && nthreads > 1) "Parallel" else "Sequential"
    thread_desc <- if (!is.null(nthreads) && nthreads > 1) paste0(" on ", nthreads, " threads") else ""
    message(mode_str, " mode: ", mode_desc, thread_desc)
  }
  
  list(pair_ids = pair_ids, pairing_info = pairing_info)
}

#' Execute divergence computation (abstracted seq vs parallel dispatch)
#' Consolidates nearly-identical sequential and parallel blocks

#' @noRd
.tsenat_compute_divergence_worker <- function(gene_indices, all_gene_names, se, gene_col, rd,
                                            group_col, control_group, q, nboot, ci, method,
                                            log_base, pseudocount, seed, pair_ids, 
                                            nthreads, use_parallel, progress) {
  start_time <- Sys.time()
  num_genes <- length(gene_indices)
  results_list <- list()
  
  if (!use_parallel) {
    # ====== SEQUENTIAL PROCESSING ======
    for (i in seq_along(gene_indices)) {
      gene_idx <- gene_indices[i]
      
      results_list[[i]] <- .tsenat_process_single_gene_div(
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
    
    parallel::clusterExport(cl, 
      c(".tsenat_process_single_gene_div", ".calculate_divergence_bootstrap", ".tsenat_tsallis_divergence_scalar",
        ".tsenat_compute_aggregate_counts", ".tsenat_extract_group_counts_gene", ".tsenat_compute_divergence_q",
        ".tsenat_bootstrap_build_args", ".tsenat_make_error_result",
        "rd", "gene_col", "all_gene_names", "se", "q", "nboot", "ci", "method",
        "log_base", "pseudocount", "seed", "group_col", "control_group", "pair_ids"),
      envir = environment())
    
    parallel::clusterCall(cl, function() {
      requireNamespace("SummarizedExperiment", quietly = TRUE)
    })
    
    results_list <- parallel::parLapply(cl, seq_along(gene_indices), function(i) {
      gene_idx <- gene_indices[i]
      .tsenat_process_single_gene_div(
        gene_idx, all_gene_names, se, gene_col, rd,
        group_col, control_group, q, nboot, ci, method,
        log_base, pseudocount, seed, pair_ids
      )
    })
  }
  
  elapsed <- as.numeric(Sys.time() - start_time, units = "secs")
  list(results = results_list, elapsed = elapsed)
}

#' Finalize result matrices from results list
#' Consolidates matrix initialization, population, and reference q handling

#' @noRd
.tsenat_finalize_divergence_matrices <- function(results_list, num_genes, q, norm, progress) {
  # Initialize result matrices
  matrices <- .tsenat_initialize_matrices(num_genes, q)
  assay_matrix <- matrices$assay
  row_data_df <- matrices$rowData
  
  # Populate matrices from results_list
  populated <- .tsenat_populate_matrices(results_list, assay_matrix, row_data_df, q)
  assay_matrix <- populated$assay
  row_data_df <- populated$rowData
  
  # Set row names
  rownames(row_data_df) <- row_data_df$gene_name
  rownames(assay_matrix) <- row_data_df$gene_name
  
  # Populate generic estimate/lower_ci/upper_ci columns using reference q value (q=1)
  q_ref <- 1.0
  q_idx <- which.min(abs(q - q_ref))
  if (length(q_idx) > 0 && q_idx <= length(q)) {
    ref_q <- q[q_idx]
    estimate_col <- paste0("estimate_q", ref_q)
    lower_ci_col <- paste0("lower_ci_q", ref_q)
    upper_ci_col <- paste0("upper_ci_q", ref_q)
    ci_width_col <- paste0("ci_width_q", ref_q)
    
    if (estimate_col %in% colnames(row_data_df)) {
      row_data_df$estimate <- row_data_df[[estimate_col]]
      row_data_df$lower_ci <- row_data_df[[lower_ci_col]]
      row_data_df$upper_ci <- row_data_df[[upper_ci_col]]
      row_data_df$ci_width <- row_data_df[[ci_width_col]]
    }
  }
  
  # Apply normalization if requested
  if (norm != "none") {
    if (progress) {
      message(sprintf("Applying '%s' normalization to divergence estimates...", norm))
    }
    
    normalized <- .tsenat_normalize_divergence_matrix(
      assay_matrix = assay_matrix,
      row_data_df = row_data_df,
      q_vals = q,
      norm = norm
    )
    assay_matrix <- normalized$assay
    row_data_df <- normalized$rowData
  }
  
  # Classify per-q patterns
  row_data_df$per_q_pattern <- NA_character_
  if (length(q) > 1) {
    for (i in seq_len(nrow(row_data_df))) {
      per_q_divs <- assay_matrix[i, ]
      names(per_q_divs) <- paste0("q_", q)
      
      if (sum(!is.na(per_q_divs)) >= 2) {
        pattern <- .classify_q_pattern(per_q_divs)
        row_data_df$per_q_pattern[i] <- if (is.na(pattern)) "UNCLASSIFIED" else pattern
      }
    }
  }
  
  list(assay = assay_matrix, rowData = row_data_df)
}

#' Print divergence computation summary
#' Consolidates logging and summary statistics reporting

#' @noRd
.tsenat_print_divergence_summary <- function(num_genes, num_errors, elapsed, row_data_df, progress) {
  num_success <- num_genes - num_errors
  
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
}


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
#' .calculate_divergence(se, res=NULL, ...)  
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
#' Matches `.calculate_diversity()` input/output pattern: transcript counts SE -> gene-level derivative SE
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

#' @noRd

.calculate_divergence <- function(
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
  # INPUT VALIDATION - Must be done BEFORE implementation
  # Following Bioconductor guidelines: fail fast with clear error messages
  # =========================================================================
  
  # Validate bootstrap parameter
  if (!is.logical(bootstrap) || length(bootstrap) != 1 || is.na(bootstrap)) {
    stop("bootstrap must be a logical", call. = FALSE)
  }
  
  # Validate paired parameter
  if (!is.logical(paired) || length(paired) != 1 || is.na(paired)) {
    stop("paired must be a logical", call. = FALSE)
  }
  
  # Validate method parameter
  if (!is.null(method) && (!is.character(method) || length(method) == 0 || is.na(method[1]))) {
    stop("method must be a character string", call. = FALSE)
  }
  
  # Validate nthreads parameter
  if (!is.null(nthreads)) {
    if (!is.numeric(nthreads) || length(nthreads) != 1 || is.na(nthreads)) {
      stop("nthreads must be a positive integer", call. = FALSE)
    }
  }
  
  # Validate progress parameter
  if (!is.logical(progress) || length(progress) != 1 || is.na(progress)) {
    stop("progress must be a logical", call. = FALSE)
  }
  
  # Call implementation directly - errors will propagate clearly
  .tsenat_calculate_divergence_impl(
    se, group_col, control_group, q, paired, bootstrap, nboot, ci, method,
    norm, log_base, pseudocount, nthreads, progress, verbose, seed
  )
}

#' Implementation of calculate_divergence with parameter validation

#' @noRd
.tsenat_calculate_divergence_impl <- function(
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
  # =========================================================================
  # INPUT VALIDATION & SETUP
  # =========================================================================

  # Parameters are validated in .calculate_divergence() 
  # Normalize/coerce for internal use
  if (!isTRUE(bootstrap)) {
    bootstrap <- FALSE
  }
  if (!isTRUE(paired)) {
    paired <- FALSE
  }
  
  method <- if (!is.null(method) && is.character(method) && length(method) > 0) {
    as.character(method[1])
  } else {
    "percentile"
  }
  
  nthreads <- if (is.numeric(nthreads) && length(nthreads) == 1 && !is.na(nthreads)) {
    as.integer(max(1, nthreads))
  } else if (is.null(nthreads)) {
    1L
  } else {
    1L
  }
  
  progress <- if (is.logical(progress) && length(progress) == 1) {
    progress
  } else {
    FALSE
  }
  
  norm <- .tsenat_validate_norm_parameter(norm)
  q <- .tsenat_validate_and_sort_q_values(q)
  .tsenat_validate_se_input(se)

  # =========================================================================
  # AUTO-DETECT GROUP COLUMN AND CONTROL GROUP
  # =========================================================================

  group_info <- .tsenat_validate_and_auto_detect_groups(se, group_col, control_group, progress)
  group_col <- group_info$group_col
  control_group <- group_info$control_group

  # =========================================================================
  # PREPARE GENES FOR PROCESSING
  # =========================================================================

  genes_info <- .tsenat_prepare_genes_processing(se)
  gene_col <- genes_info$gene_col
  all_gene_names <- genes_info$all_gene_names
  gene_indices <- genes_info$gene_indices
  rd <- genes_info$rd
  num_genes <- genes_info$num_genes

  # =========================================================================
  # BOOTSTRAP & PARALLEL CONFIGURATION
  # =========================================================================

  boot_config <- .tsenat_bootstrap_configure_parallel(
    bootstrap, nboot, method, num_genes, nthreads, progress
  )
  nboot <- boot_config$nboot
  nthreads <- boot_config$nthreads
  use_parallel <- boot_config$use_parallel

  # =========================================================================
  # PAIRED SAMPLE DETECTION & EXECUTION SETUP
  # =========================================================================

  exec_setup <- .tsenat_prepare_divergence_execution(
    se, bootstrap, paired, nboot, method, nthreads, progress
  )
  pair_ids <- exec_setup$pair_ids

  # =========================================================================
  # EXECUTE DIVERGENCE COMPUTATION (SEQUENTIAL OR PARALLEL)
  # =========================================================================

  comp_result <- .tsenat_compute_divergence_worker(
    gene_indices, all_gene_names, se, gene_col, rd,
    group_col, control_group, q, nboot, ci, method,
    log_base, pseudocount, seed, pair_ids,
    nthreads, use_parallel, progress
  )
  results_list <- comp_result$results
  elapsed <- comp_result$elapsed

  # =========================================================================
  # FINALIZE MATRICES & APPLY NORMALIZATION
  # =========================================================================

  matrices_final <- .tsenat_finalize_divergence_matrices(
    results_list, num_genes, q, norm, progress
  )
  assay_matrix <- matrices_final$assay
  row_data_df <- matrices_final$rowData

  # =========================================================================
  # SUMMARY STATISTICS & LOGGING
  # =========================================================================

  num_errors <- sum(!is.na(row_data_df$error))
  .tsenat_print_divergence_summary(num_genes, num_errors, elapsed, row_data_df, progress)

  # =========================================================================
  # CREATE & RETURN SUMMARIZED EXPERIMENT
  # =========================================================================

  result_se <- .tsenat_construct_result_se(
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
.tsenat_validate_norm_parameter <- function(norm) {
    if (is.logical(norm)) {
        norm <- if (norm) "range" else "none"
    }
    match.arg(norm, choices = c("none", "range", "zscore", 
                                 "log_odds_ratio", "relative_reference"))
}

#' Validate and sort q-parameter values
#' Ensures q >= 0 and returns sorted vector
#' @noRd
.tsenat_validate_and_sort_q_values <- function(q) {
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
.tsenat_validate_se_input <- function(se) {
    if (!methods::is(se, "SummarizedExperiment")) {
        stop("se must be a SummarizedExperiment object")
    }
    TRUE
}

#' Identify gene name/ID column in rowData
#' Preference: gene_name (human-readable) > gene_id (ensembl)
#' @noRd
.tsenat_identify_gene_column <- function(se) {
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
.tsenat_extract_gene_list <- function(se, gene_col) {
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
.tsenat_configure_parallel <- function(nthreads, num_genes) {
    if (is.null(nthreads) || is.na(nthreads)) {
        nthreads <- parallel::detectCores() - 1
        nthreads <- max(1, nthreads)
    }
    
    if (!is.numeric(nthreads) || nthreads < 1) {
        stop("'nthreads' must be a positive integer")
    }
    nthreads <- as.integer(nthreads)
    
    # Safe boolean check: only use parallel if num_genes is numeric and nthreads > 1
    use_parallel <- (!is.na(num_genes) && num_genes >= 5 && nthreads > 1)
    
    list(nthreads = nthreads, use_parallel = use_parallel)
}

# GENE PROCESSING HELPERS
# ============================================================================

#' Aggregate transcript-level counts to gene-level
#' Sums counts across all transcripts for a given gene
#' @noRd
.tsenat_compute_aggregate_counts <- function(se, target_gene, gene_col, rd) {
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
.tsenat_extract_group_counts_gene <- function(counts_gene, groups, control_group) {
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
.tsenat_make_error_result <- function(gene_name, q_vals, error_msg, elapsed_sec = NA_real_) {
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
.tsenat_bootstrap_build_args <- function(x, y, q_val, nboot, ci, method, 
                                   log_base, pseudocount, gene_name, 
                                   seed, pair_ids = NULL) {
    args <- list(
        x = x, y = y, q = q_val, nboot = nboot, ci = ci, method = method,
        log_base = log_base, pseudocount = pseudocount,
        gene_name = gene_name, verbose = FALSE, seed = seed,
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
.tsenat_compute_divergence_q <- function(x, y, q_vals, nboot, ci, method,
                                       log_base, pseudocount, gene_name,
                                       seed, pair_ids = NULL) {
    gene_results <- list()
    
    for (j in seq_along(q_vals)) {
        q_val <- q_vals[j]
        
        bootstrap_args <- .tsenat_bootstrap_build_args(
            x, y, q_val, nboot, ci, method, 
            log_base, pseudocount, gene_name, 
            seed, pair_ids
        )
        
        result <- do.call(.calculate_divergence_bootstrap, bootstrap_args)
        gene_results[[j]] <- result
    }
    
    gene_results
}

# RESULTS COMPILATION HELPERS
# ============================================================================

#' Initialize result matrices for results compilation
#' Creates assay matrix and rowData structure
#' @noRd
.tsenat_initialize_matrices <- function(num_genes, q_vals) {
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
.tsenat_populate_matrices <- function(results_list, assay_matrix, row_data_df, q_vals) {
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
.tsenat_normalize_range_matrix <- function(assay_matrix, row_data_df, q_vals) {
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
.tsenat_divergence_normalize_zscore <- function(assay_matrix, row_data_df, q_vals) {
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
.tsenat_divergence_normalize_log_odds_ratio <- function(assay_matrix, row_data_df, q_vals) {
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
.tsenat_normalize_reference <- function(assay_matrix, row_data_df, q_vals) {
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
.tsenat_normalize_divergence_matrix <- function(assay_matrix, row_data_df, q_vals, norm) {
    if (norm == "none") {
        return(list(assay = assay_matrix, rowData = row_data_df))
    }
    
    if (norm == "range") {
        return(.tsenat_normalize_range_matrix(assay_matrix, row_data_df, q_vals))
    }
    
    if (norm == "zscore") {
        return(.tsenat_divergence_normalize_zscore(assay_matrix, row_data_df, q_vals))
    }
    
    if (norm == "log_odds_ratio") {
        return(.tsenat_divergence_normalize_log_odds_ratio(assay_matrix, row_data_df, q_vals))
    }
    
    if (norm == "relative_reference") {
        return(.tsenat_normalize_reference(assay_matrix, row_data_df, q_vals))
    }
    
    list(assay = assay_matrix, rowData = row_data_df)
}

#' Generate computation summary statistics
#' Formats timing information for progress messages
#' @noRd
.tsenat_generate_summary <- function(elapsed, num_genes, num_errors, row_data_df) {
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
.tsenat_construct_result_se <- function(assay_matrix, row_data_df, q_vals, 
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
.tsenat_process_single_gene_div <- function(gene_idx, all_gene_names, se, gene_col, rd, 
                                 group_col, control_group, q, nboot, ci, method,
                                 log_base, pseudocount, seed, pair_ids) {
    target_gene <- all_gene_names[gene_idx]
    gene_name <- target_gene
    gene_start <- Sys.time()
    
    tryCatch({
        # Get gene-level counts via transcript aggregation
        counts_gene <- .tsenat_compute_aggregate_counts(se, target_gene, gene_col, rd)
        
        if (is.null(counts_gene)) {
            return(.tsenat_make_error_result(gene_name, q, "No transcripts found for gene"))
        }
        
        # Extract group-specific counts
        groups <- se[[group_col]]
        group_counts <- .tsenat_extract_group_counts_gene(counts_gene, groups, control_group)
        x <- group_counts$control
        y <- group_counts$treatment
        
        if (length(x) == 0 || length(y) == 0) {
            return(.tsenat_make_error_result(gene_name, q, "Insufficient group samples"))
        }
        
        # Compute divergence for each q value
        gene_results <- .tsenat_compute_divergence_q(x, y, q, nboot, ci, method,
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
        .tsenat_make_error_result(gene_name, q, as.character(e$message),
                          as.numeric(Sys.time() - gene_start, units = "secs"))
    })
}
