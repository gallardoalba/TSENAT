# ============================================================================
# DETECT Q GENE INTERACTIONS WRAPPER
# ============================================================================

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
#' # Build analysis from vignette data and create manageable subset
#' analysis <- build_analysis_s4(readcounts, gff3_dataset, metadata = metadata_df,
#'   tpm = salmon_tpm, effective_length = salmon_effective_length)
#' analysis <- filter_analysis_s4(analysis, min_samples = 1, subset_n_genes = 200)
#' analysis <- calculate_diversity_s4(analysis, q = c(0.5, 1.0, 1.5), verbose = FALSE)
#' 
#' # Test Q×Condition interaction (condition_col is REQUIRED)
#' analysis <- rank_test_q_condition_s4(analysis, condition_col = "condition", 
#'                                            multicorr = "hochberg")
#' # View results
#' results <- lmResults(analysis)
#' if (!is.null(results)) head(results$q_interactions)
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
  
  # PHASE 1: Validate input and prerequisites
  condition_col <- .validate_rank_test_input(analysis, condition_col)
  
  # PHASE 2: Resolve parameters from config + explicit args
  param_result <- .resolve_rank_test_params(analysis, test, multicorr, nperm_mode,
                                            q, paired, subject_col, nthreads, 
                                            wy_randomizations, entropy_col, q_col, gene_col)
  dots <- param_result$dots
  dots$condition_col <- condition_col
  dots$verbose <- verbose
  
  # PHASE 3: Prepare multi-Q SummarizedExperiment
  se_multi_q <- .prepare_multi_q_se(analysis)
  
  # PHASE 4: Run core rank-based testing
  result <- tryCatch({
    do.call(.rank_test_q_condition, c(list(data = se_multi_q), dots))
  }, error = function(e) {
    stop("q-interaction detection failed:\n", e$message, call. = FALSE)
  })
  
  # PHASE 5: Store results and save if requested
  analysis <- .store_rank_test_results(analysis, result, output_file, verbose)
  
  analysis
}

#' Internal: Validate rank test input and prerequisites
#'
#' @noRd
.validate_rank_test_input <- function(analysis, condition_col) {
  if (!is(analysis, "TSENATAnalysis")) {
    stop("'analysis' must be a TSENATAnalysis object", call. = FALSE)
  }
  
  if (missing(condition_col) || is.null(condition_col)) {
    if (!is.null(analysis@config) && "condition_col" %in% names(analysis@config)) {
      condition_col <- analysis@config$condition_col
    } else {
      condition_col <- "condition"
    }
  }
  
  if (length(analysis@diversity_results) == 0) {
    stop("Diversity results required. Run calculate_diversity_s4() first.",
         call. = FALSE)
  }
  
  condition_col
}

#' Internal: Resolve rank test parameters from config
#'
#' @noRd
.resolve_rank_test_params <- function(analysis, test, multicorr, nperm_mode, 
                                     q, paired, subject_col, nthreads, wy_randomizations,
                                     entropy_col, q_col, gene_col) {
  dots <- list()
  
  # Match enums early
  if (!missing(test)) {
    test <- match.arg(test, c("auto", "kruskal-wallis", "friedman", "art"))
    dots$test <- test
  } else if ("test" %in% names(analysis@config)) {
    dots$test <- analysis@config$test
  }
  
  if (!missing(multicorr)) {
    multicorr <- match.arg(multicorr, c("hochberg", "benjamini-yekutieli", "westfall-young", "none"))
    dots$multicorr <- multicorr
  } else if ("multicorr" %in% names(analysis@config)) {
    dots$multicorr <- analysis@config$multicorr
  }
  
  if (!missing(nperm_mode)) {
    nperm_mode <- match.arg(nperm_mode, c("standard", "conservative", "interactive"))
    dots$nperm_mode <- nperm_mode
  } else if ("nperm_mode" %in% names(analysis@config)) {
    dots$nperm_mode <- analysis@config$nperm_mode
  }
  
  # Add column parameters
  dots$entropy_col <- entropy_col
  dots$q_col <- q_col
  dots$gene_col <- gene_col
  
  # Use resolve_slot_param for remaining parameters
  q <- resolve_slot_param(q, analysis@config, "q_values", NULL)
  paired <- resolve_slot_param(paired, analysis@config, "paired", NULL)
  subject_col <- resolve_slot_param(subject_col, analysis@config, "subject_col", NULL)
  nthreads <- resolve_slot_param(nthreads, analysis@config, "nthreads", 1)
  
  if (!is.null(paired)) dots$paired <- paired
  if (!is.null(subject_col)) dots$subject_col <- subject_col
  dots$nthreads <- nthreads
  dots$wy_randomizations <- wy_randomizations
  
  list(dots = dots, q_extracted = q)
}

#' Internal: Prepare multi-Q SummarizedExperiment for testing
#'
#' @noRd
.prepare_multi_q_se <- function(analysis) {
  # Check cache first
  if (!is.null(analysis@metadata$diversity_combined) && 
      is.list(analysis@metadata$diversity_combined) &&
      !is.null(analysis@metadata$diversity_combined$combined_se)) {
    
    se_multi_q <- analysis@metadata$diversity_combined$combined_se
    if (is(se_multi_q, "SummarizedExperiment") && ncol(se_multi_q) > 0) {
      return(se_multi_q)
    }
  }
  
  # Fallback: combine per-Q results
  q_keys <- names(analysis@diversity_results)
  combined_assay_list <- list()
  combined_coldata_list <- list()
  common_rownames <- NULL
  
  for (key in sort(q_keys)) {
    se <- analysis@diversity_results[[key]]
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
    
    assay_data <- SummarizedExperiment::assay(se, 1)
    assay_rownames <- rownames(assay_data)
    
    if (is.null(assay_rownames)) {
      assay_rownames <- paste0("gene_", seq_len(nrow(assay_data)))
    }
    if (is.null(common_rownames)) {
      common_rownames <- assay_rownames
    } else if (!identical(common_rownames, assay_rownames)) {
      if (length(common_rownames) == length(assay_rownames)) {
        assay_data <- assay_data[common_rownames, , drop = FALSE]
      } else {
        stop("Diversity result for ", key, " has different number of genes", call. = FALSE)
      }
    }
    rownames(assay_data) <- common_rownames
    
    # Rename columns with q-value suffix
    orig_colnames <- colnames(assay_data)
    if (is.null(orig_colnames)) orig_colnames <- paste0("sample_", seq_len(ncol(assay_data)))
    unique_colnames <- paste0(orig_colnames, "_q=", q_val)
    colnames(assay_data) <- unique_colnames
    
    # Get and update colData
    cd <- as.data.frame(SummarizedExperiment::colData(se))
    if (nrow(cd) == 0) {
      cd <- data.frame(q = rep(q_val, ncol(assay_data)))
    } else if (!"q" %in% colnames(cd)) {
      cd$q <- q_val
    }
    rownames(cd) <- unique_colnames
    
    combined_assay_list[[key]] <- assay_data
    combined_coldata_list[[key]] <- cd
  }
  
  # Combine horizontally
  combined_assay <- do.call(cbind, combined_assay_list)
  combined_coldata_df <- do.call(rbind, combined_coldata_list)
  colnames(combined_assay) <- rownames(combined_coldata_df)
  
  # Get rowData from first SE
  first_se <- analysis@diversity_results[[sort(q_keys)[1]]]
  if (!is(first_se, "SummarizedExperiment")) {
    first_se <- SummarizedExperiment(assays = list(entropy = as.matrix(first_se)))
  }
  rd <- tryCatch(SummarizedExperiment::rowData(first_se), error = function(e) NULL)
  
  se_multi_q <- SummarizedExperiment(assays = list(entropy = combined_assay),
                                      colData = combined_coldata_df)
  if (!is.null(rd) && nrow(rd) > 0) {
    SummarizedExperiment::rowData(se_multi_q) <- rd
  }
  
  se_multi_q
}

#' Internal: Store rank test results in analysis object
#'
#' @noRd
.store_rank_test_results <- function(analysis, result, output_file, verbose) {
  if (is.list(analysis@lm_results)) {
    analysis@lm_results$q_interactions <- result
  } else {
    analysis@lm_results <- list(q_interactions = result)
  }
  
  if (!is.null(output_file)) {
    result_df <- as.data.frame(result)
    save_analysis_output(result_df, output_file, object = analysis, verbose = verbose,
                         func_name = "rank_test_q_condition_s4")
  }
  
  analysis
}
