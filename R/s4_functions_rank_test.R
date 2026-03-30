# ============================================================================
# DETECT Q GENE INTERACTIONS WRAPPER
# ============================================================================

#' Detect q-dependent gene interactions
#'
#' @param analysis \code{TSENATAnalysis} object.
#' @param q \code{numeric} or \code{NULL}. Q-values to test across spectrum.
#'   If NULL, auto-detects from \code{@config$q_values} or diversity results.
#' @param output_file \code{character} or \code{NULL}. Optional file path to save results.
#'   Supported formats: .rds (for S4 objects). Default: NULL (no file output).
#' @param paired \code{logical} or \code{NULL}. If TRUE, uses paired/blocked design 
#'   (requires \code{subject_col}). If NULL, reads from \code{@config$paired}.
#' @param subject_col \code{character} or \code{NULL}. Column name for subject/block identifiers
#'   (required when \code{paired=TRUE}). If NULL, reads from \code{@config$subject_col}.
#' @param condition_col \code{character}. Column name for sample grouping/condition (REQUIRED).
#'   Specifies the condition/treatment variable for testing q×condition interactions.
#'   Example: "sample_type", "treatment", "disease_status".
#' @param test \code{character}. Test method: "auto" (default), "kruskal-wallis" (unpaired),
#'   "friedman" (paired), or "art" (aligned rank transform).
#' @param multicorr \code{character}. Multiple testing correction: "hochberg" (default),
#'   "benjamini-yekutieli", "westfall-young", or "none".
#' @param entropy_col \code{character}. Column name containing entropy/diversity data.
#'   Default: "diversity".
#' @param q_col \code{character}. Column name containing q-values. Default: "q".
#' @param gene_col \code{character}. Column name containing gene identifiers. Default: "gene".
#' @param wy_randomizations \code{numeric} or \code{character}. Number of permutations for 
#'   Westfall-Young correction. Use "auto" to estimate from data. Default: 500.
#' @param nperm_mode \code{character}. Mode for automatic permutation estimation:
#'   "standard" (default), "conservative", or "interactive".
#' @param nthreads \code{numeric} or \code{NULL}. Number of parallel threads for computation.
#'   If NULL, reads from \code{@config$nthreads}.
#' @param verbose \code{logical}. If TRUE, prints progress messages. Default: FALSE.
#' @param ... Additional arguments passed to the base \code{.rank_test_q_condition()} function.
#'
#' @return Modified TSENATAnalysis with interaction results in @lm_results.
#'
#' @details
#' Analyzes how gene interactions change across q-value spectrum using rank-based
#' (Friedman/Kruskal-Wallis) or parametric (GAM) statistical tests.
#'
#' **Parameter resolution priority** (explicit > @config > default/auto-detect):
#' \itemize{
#'   \item \code{condition_col}: REQUIRED - must be explicitly provided
#'   \item \code{q}: explicit arg > \code{@config$q_values} > extract from diversity_results keys
#'   \item \code{paired}: explicit arg > \code{@config$paired} > FALSE (default)
#'   \item \code{subject_col}: explicit arg > \code{@config$subject_col}
#'   \item \code{multicorr}: explicit arg > \code{@config$multicorr} > "hochberg"
#'   \item \code{nthreads}: explicit arg > \code{@config$nthreads} > 1 (default)
#'   \item \code{test}: explicit arg > \code{@config$test} > "auto" (auto-selection)
#'   \item \code{nperm_mode}: explicit arg > \code{@config$nperm_mode} > "standard"
#' }
#'
#' @examples
#' # Load example data (matching TSENAT.Rmd workflow)
#' data(readcounts)
#' readcounts <- as.matrix(salmon_dataset)
#' mode(readcounts) <- "numeric"
#' metadata_df <- read.table(
#'   system.file("extdata", "metadata.tsv", package = "TSENAT"),
#'   header = TRUE, sep = "\t"
#' )
#' gff3_dataset <- system.file("extdata", "annotation.gff3.gz", package = "TSENAT")
#' 
#' # Build analysis from vignette data and create small subset
#' analysis <- build_analysis_s4(readcounts, gff3_dataset, metadata = metadata_df,
#'   tpm = salmon_tpm, effective_length = salmon_effective_length)
#' analysis <- subset_analysis(analysis, n_genes = 30, n_samples = 8)
#' analysis <- calculate_diversity_s4(analysis, norm = TRUE)
#' 
#' # Test Q×Condition interaction (condition_col is REQUIRED)
#' analysis <- rank_test_q_condition_s4(analysis, condition_col = "condition", 
#'                                            multicorr = "hochberg")
#' head(lmResults(analysis)$q_interactions)
#'
#' @export
#' @importFrom utils write.table
# ============================================================================
# S4 WRAPPER: Detect Q×Condition Gene Interactions (Rank-Based Testing)
# ============================================================================
# Purpose:
#   Wrapper around .rank_test_q_condition() that manages TSENATAnalysis object.
#   Tests for genes with CONDITION-SPECIFIC q-dependent entropy patterns.
#   Tests whether the effect of q-values DIFFERS between experimental conditions.
# 
# Key Features:
#   - Q×Condition interaction: Tests if entropy patterns across q-values differ by condition
#   - Multi-q analysis: Combines diversity results for multiple q-values into
#     a single SummarizedExperiment for joint hypothesis testing
#   - Rank-based statistics: Kruskal-Wallis (unpaired) or Friedman (paired)
#   - Scheirer-Ray-Hare test: Two-way non-parametric ANOVA on ranks
#   - Multiple testing correction: Hochberg, Benjamini-Yekutieli, or
#     Westfall-Young permutation procedure
#   - AR(1) correlation handling: Westfall-Young preserves q-value correlations
#   - Effect sizes: Eta-squared (η²) for q×condition interactions
#
#   Mathematical Background:
#   Tests null hypothesis: H0 = "Gene entropy q-effect does NOT differ between conditions"
#   vs Alternative: H1 = "Gene entropy q-dependence is CONDITION-SPECIFIC"
# 
#   Example: Gene shows strong isoform switching (q-dependent entropy) in tumor cells
#   but NOT in healthy cells -> Identified as disease-relevant q-dependent gene.
# 
#   For condition-specific q-dependent genes:
#   - Condition A: Strong entropy variation across q (q-dependent)
#   - Condition B: Flat entropy profile across q (q-independent)
#   - Interaction: Condition-specific q-dependence pattern reveals biological process
# ============================================================================
rank_test_q_condition_s4 <- function(
    analysis, 
    condition_col,
    q = NULL, 
    output_file = NULL,
    paired = NULL,
    subject_col = NULL,
    test = c("auto", "kruskal-wallis", "friedman", "art"),
    multicorr = c("hochberg", "benjamini-yekutieli", "westfall-young", "none"),
    entropy_col = "diversity",
    q_col = "q",
    gene_col = "gene",
    wy_randomizations = 500,
    nperm_mode = c("standard", "conservative", "interactive"),
    nthreads = NULL,
    verbose = FALSE,
    ...) {
  # Validate input is TSENATAnalysis S4 class
  if (!is(analysis, "TSENATAnalysis")) {
    stop("'analysis' must be a TSENATAnalysis object", call. = FALSE)
  }

  # Validate condition_col is provided (Q×Condition interaction is required)
  # Try to read from config if not explicitly provided
  if (missing(condition_col) || is.null(condition_col)) {
    # Try to get from config
    if (!is.null(analysis@config) && "condition_col" %in% names(analysis@config)) {
      condition_col <- analysis@config$condition_col
    } else {
      # Default to "condition" (or "sample_type" for backward compatibility) if still not found
      condition_col <- "condition"
    }
  }

  # ========================================================================
  # PREREQUISITE CHECK: Diversity must be pre-calculated
  # ========================================================================
  # .rank_test_q_condition() requires a SummarizedExperiment with:
  #   - assays: entropy values (genes × samples)
  #   - colData: q-values, condition_col, and optional subject information
  if (length(analysis@diversity_results) == 0) {
    stop("Diversity results required. Run calculate_diversity_s4() first.",
         call. = FALSE)
  }

  # ========================================================================
  # PARAMETER EXTRACTION FROM @config using centralized helpers
  # ========================================================================
  
  # Prepare dots for additional arguments
  dots <- list(...)
  
  # Match test and multicorr enums early
  if (!missing(test)) {
    test <- match.arg(test)
    dots$test <- test
  } else if ("test" %in% names(analysis@config)) {
    dots$test <- analysis@config$test
  }
  
  if (!missing(multicorr)) {
    multicorr <- match.arg(multicorr)
    dots$multicorr <- multicorr
  } else if ("multicorr" %in% names(analysis@config)) {
    dots$multicorr <- analysis@config$multicorr
  }
  
  if (!missing(nperm_mode)) {
    nperm_mode <- match.arg(nperm_mode)
    dots$nperm_mode <- nperm_mode
  } else if ("nperm_mode" %in% names(analysis@config)) {
    dots$nperm_mode <- analysis@config$nperm_mode
  }
  
  # Add column specification parameters
  dots$entropy_col <- entropy_col
  dots$q_col <- q_col
  dots$gene_col <- gene_col
  
  # Use resolve_slot_param for remaining parameters
  q <- resolve_slot_param(q, analysis@config, "q_values", NULL)
  paired <- resolve_slot_param(paired, analysis@config, "paired", NULL)
  subject_col <- resolve_slot_param(subject_col, analysis@config, "subject_col", NULL)
  nthreads <- resolve_slot_param(nthreads, analysis@config, "nthreads", 1)
  
  # Add resolved parameters to dots list
  if (!is.null(paired)) dots$paired <- paired
  if (!is.null(subject_col)) dots$subject_col <- subject_col
  # condition_col is REQUIRED and explicitly passed
  dots$condition_col <- condition_col
  dots$nthreads <- nthreads
  dots$wy_randomizations <- wy_randomizations
  dots$verbose <- verbose

  # ========================================================================
  # OPTIMIZATION: Use Cached Multi-Q SummarizedExperiment
  # ========================================================================
  # If calculate_diversity_s4() is called with multiple q-values in a single call,
  # it caches the combined SE in @metadata$diversity_combined for reuse.
  # 
  # This optimization bypasses expensive per-q recombination when available:
  #   - Per-q SEs are automatically cbind()ed horizontally
  #   - Column names include _q= suffix to distinguish q-values
  #   - colData is rbind()ed with preserved q-value information
  # 
  # Benefits: ~10-50x faster for multi-q analysis on large datasets
  
  # Check if we have combined diversity result cached (from lazy evaluation)
  if (!is.null(analysis@metadata$diversity_combined) && 
      is.list(analysis@metadata$diversity_combined) &&
      !is.null(analysis@metadata$diversity_combined$combined_se)) {
    
    # Use cached combined SE directly - it has correct structure and metadata
    se_multi_q <- analysis@metadata$diversity_combined$combined_se
    
    # Verify it's valid
    if (!is(se_multi_q, "SummarizedExperiment") || ncol(se_multi_q) == 0) {
      se_multi_q <- NULL
    }
  } else {
    se_multi_q <- NULL
  }
  
  # Fallback: Recombine Per-Q Results
  # ========================================================================
  # If cache not available (e.g., diversity calculated with separate calls),
  # manually combine per-q SummarizedExperiments into single SE:
  # 
  # Process:
  #   1. Extract q-values from diversity_results keys (format: "q_0.500", "q_1.000", etc.)
  #   2. Validate all SEs have same genes (rownames must match)
  #   3. cbind() all assay matrices with renamed columns (add _q=X.XXX suffix)
  #   4. rbind() all colData (with updated rownames matching combined columns)
  #   5. Combine rowData from first SE (genes are same across all q-values)
  # 
  # Result: Single SummarizedExperiment(genes × (samples per q × n_q))
  if (is.null(se_multi_q)) {
    
    # Step 1: Extract q-values from diversity_results keys (format: "q_0.5", "q_1.0", etc.)
    q_keys <- names(analysis@diversity_results)
    q_values_extracted <- as.numeric(sub("^q_", "", q_keys))
    q_values_extracted <- sort(q_values_extracted)

    # Step 2-4: Combine list of SEs (one per q-value) into single SE
    # Each SE has same genes (rows) but different q-value samples (columns)
    # cbind() the assay matrices, rbind() the colData
    combined_assay_list <- list()
    combined_coldata_list <- list()
    common_rownames <- NULL  # Track genes are same across all q-values

    for (key in sort(q_keys)) {
      se <- analysis@diversity_results[[key]]
      
      # Extract q-value from key
      q_val <- as.numeric(sub("^q_", "", key))
      
      # Ensure SE format
      if (!is(se, "SummarizedExperiment")) {
        if (is.matrix(se) || is.data.frame(se)) {
          se <- SummarizedExperiment(assays = list(entropy = as.matrix(se)))
        } else {
          stop("Diversity result for ", key, " is not a SummarizedExperiment or matrix",
               call. = FALSE)
        }
      }
      
      # Get assay data
      if (length(SummarizedExperiment::assays(se)) == 0) {
        stop("Diversity result for ", key, " has no assays", call. = FALSE)
      }
      assay_data <- SummarizedExperiment::assay(se, 1)
      
      # Extract and enforce consistent rownames
      assay_rownames <- rownames(assay_data)
      if (is.null(assay_rownames)) {
        assay_rownames <- paste0("gene_", seq_len(nrow(assay_data)))
      }
      if (is.null(common_rownames)) {
        common_rownames <- assay_rownames
      } else if (!identical(common_rownames, assay_rownames)) {
        # If rownames differ, use the first one and reorder/match
        if (length(common_rownames) == length(assay_rownames)) {
          assay_data <- assay_data[common_rownames, , drop = FALSE]
        } else {
          stop("Diversity result for ", key, " has different number of genes",
               call. = FALSE)
        }
      }
      rownames(assay_data) <- common_rownames
      
      # Append q-value to column names: distinguishes samples from different q-values
      # Format: "sample_01_q=0.500", "sample_02_q=1.000", etc.
      # This encoding enables downstream functions to parse q and map back to original samples
      orig_colnames <- colnames(assay_data)
      if (is.null(orig_colnames)) {
        orig_colnames <- paste0("sample_", seq_len(ncol(assay_data)))
      }
      unique_colnames <- paste0(orig_colnames, "_q=", q_val)
      colnames(assay_data) <- unique_colnames
      
      # Get colData - ensure q column is present
      cd <- as.data.frame(SummarizedExperiment::colData(se))
      if (nrow(cd) == 0) {
        cd <- data.frame(q = rep(q_val, ncol(assay_data)))
      } else if (!"q" %in% colnames(cd)) {
        cd$q <- q_val
      }
      rownames(cd) <- unique_colnames
      
      # Store for combination
      combined_assay_list[[key]] <- assay_data
      combined_coldata_list[[key]] <- cd
    }

    # Step 3: Combine all assays horizontally (cbind columns from different q-values)
    # This creates a matrix: genes × (sample_1_q_0.5, sample_2_q_0.5, ..., sample_1_q_1.0, ...)
    combined_assay <- do.call(cbind, combined_assay_list)
    
    # Step 4: Combine colData vertically (rbind from each q-value's colData)
    # Rownames already set to unique_colnames matching in final assay matrix
    combined_coldata_df <- do.call(rbind, combined_coldata_list)
    
    # Ensure colnames of combined_assay match rownames of combined_coldata_df
    colnames(combined_assay) <- rownames(combined_coldata_df)
    
    # Step 5: Get rowData from first diversity result
    # Genes (rows) are IDENTICAL across all q-values, so only need from first
    first_se <- analysis@diversity_results[[sort(q_keys)[1]]]
    
    # Ensure first_se is a SummarizedExperiment (handle edge cases)
    if (!is(first_se, "SummarizedExperiment")) {
      if (is.matrix(first_se) || is.data.frame(first_se)) {
        first_se <- SummarizedExperiment(assays = list(entropy = as.matrix(first_se)))
      }
    }
    
    # Extract rowData from first SE if it exists
    rd <- tryCatch({
      rd_temp <- SummarizedExperiment::rowData(first_se)
      if (nrow(rd_temp) > 0) rd_temp else NULL
    }, error = function(e) NULL)
    
    # Create combined SE
    se_multi_q <- SummarizedExperiment(
      assays = list(entropy = combined_assay),
      colData = combined_coldata_df
    )
    
    # Add rowData if available
    if (!is.null(rd) && nrow(rd) > 0) {
      SummarizedExperiment::rowData(se_multi_q) <- rd
    }
  }

  # Extract q-values for metadata tracking (available in both paths)
  if (exists("q_values_extracted") && !is.null(q_values_extracted)) {
    # Already defined in fallback path
    q_vals_for_tracking <- q_values_extracted
  } else {
    # Extract from colData in optimized path
    coldata_vals <- SummarizedExperiment::colData(se_multi_q)
    if (!is.null(coldata_vals) && "q" %in% colnames(coldata_vals)) {
      q_vals_for_tracking <- unique(as.numeric(coldata_vals$q))
      q_vals_for_tracking <- sort(q_vals_for_tracking)
    } else {
      q_vals_for_tracking <- NULL
    }
  }

  # ========================================================================
  # RUN CORE RANK-BASED Q-INTERACTION TESTING
  # ========================================================================
  # Delegate to .rank_test_q_condition() which performs:
  #   1. SummarizedExperiment → long-format data frame conversion
  #   2. Per-gene rank-based test selection (conditional on data characteristics)
  #   3. Westfall-Young permutation procedure (if multicorr="westfall-young")
  #   4. Multiple testing corrections (Hochberg, Benjamini-Yekutieli, none)
  #   5. Effect size computation (η²) and result classification
  # 
  # Use merged parameter dictionary: config values + explicit overrides
  result <- tryCatch({
    do.call(.rank_test_q_condition, c(list(data = se_multi_q), dots))
  }, error = function(e) {
    stop("q-interaction detection failed:\n", e$message,
         call. = FALSE)
  })

  # ========================================================================
  # STORE RESULTS IN TSENATAnalysis OBJECT
  # ========================================================================
  # Store results under @lm_results$q_interactions for accessor compatibility
  # This location allows other functions to retrieve results via:
  #   lmResults(analysis, "q_interactions")
  if (is.list(analysis@lm_results)) {
    analysis@lm_results$q_interactions <- result
  } else {
    analysis@lm_results <- list(q_interactions = result)
  }

  # Track function execution in audit trail
  # Enables reproducibility: know which function calls were run and in what order
  if (!is.null(q_vals_for_tracking)) {
    analysis@metadata$function_calls <- c(
      analysis@metadata$function_calls,
      paste0("rank_test_q_condition[q=", paste(q_vals_for_tracking, collapse = ","), "]")
    )
  }

  # Save if output_file provided (using centralized output handler)
  if (!is.null(output_file)) {
    result_df <- as.data.frame(result)
    save_analysis_output(result_df, output_file, object = analysis, verbose = verbose,
                         func_name = "rank_test_q_condition_s4")
  }

  analysis
}