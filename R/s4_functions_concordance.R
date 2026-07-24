#' Compare method concordance for differential analysis results
#'
#' Compares statistical results from two different methods (typically SAIT/GAM for
#' continuous data and Conover-Iman Rank Transform tests) to assess agreement and identify
#' genes detected by one method but not the other.
#'
#' @aliases calculate_concordance,TSENATAnalysis-method
#'
#' @param analysis_sait \code{TSENATAnalysis} object containing SAIT/GAM analysis results
#'   (from \code{calculate_sait()}).
#' @param analysis_rank \code{TSENATAnalysis} object or NULL. If NULL, uses legacy 
#'   single-object API with analysis_sait containing both results. If provided, 
#'   compares SAIT results from analysis_sait with rank-test results from analysis_rank.
#' @param verbose \code{logical}. Print progress messages (default: FALSE).
#' @param output_file \code{character} or NULL. Optional file path to save results.
#'   Supported formats: .rds (for S4 objects). Default: NULL (no file output).
#' @param ... Additional arguments for future extensibility.
#'
#' @return Modified TSENATAnalysis object with concordance results stored in:
#'   \code{@metadata$method_concordance}:
#'   \describe{
#'     \item{comparison_df}{Data frame comparing results from both methods}
#'     \item{spearman_rho}{Spearman correlation between adjusted p-values}
#'     \item{high_confidence}{Genes with strong agreement}
#'     \item{agreement_table}{Contingency table of significant/non-significant calls}
#'     \item{sait_method}{Method name used for SAIT/GAM analysis}
#'     \item{rank_method}{Method name used for rank-based analysis}
#'     \item{timestamp}{When concordance was computed}
#'   }
#'
#' @details
#' Compares results from two different statistical methods (typically GAM
#' for continuous and Conover-Iman Rank Transform for rank-based analysis) on the same data.
#' Identifies:
#' - Genes significant in both methods (high confidence)
#' - Genes detected by one method only (potential false positives or method-specific signal)
#' - Spearman correlation of p-values (overall agreement trends)
#'
#' @usage
#' calculate_concordance(analysis_sait, analysis_rank = NULL, ...)
#'
#' \S4method{calculate_concordance}{TSENATAnalysis}(
#'   analysis_sait,
#'   analysis_rank = NULL,
#'   verbose = FALSE,
#'   output_file = NULL,
#'   ...
#' )
#'
#' @examples
#' # Compare results from SAIT and rank-based testing
#' # (Requires pre-computed analysis objects from calculate_sait and calculate_rank_transform)
#' # results_df <- results(calculate_concordance(analysis_sait, analysis_rank))
#'
#' @aliases calculate_concordance
#' @export
setGeneric("calculate_concordance", function(analysis_sait, analysis_rank = NULL, ...) {
    standardGeneric("calculate_concordance")
})

#' @rdname calculate_concordance
#' @keywords internal
#' @param analysis_sait TSENATAnalysis object
#' @param analysis_rank TSENATAnalysis object or NULL
#' @param ... Additional arguments
#' @return NULL (stops on error)
#' @noRd
.validate_concordance_inputs <- function(analysis_sait, analysis_rank, ...) {
    # Check for unexpected arguments
    extra_args <- list(...)
    if (length(extra_args) > 0) {
        arg_names <- paste(names(extra_args), collapse = ", ")
        stop("The following argument(s) are not recognized and cannot be used: ", 
             arg_names, call. = FALSE)
    }
    
    # Validate analysis_sait
    if (!is(analysis_sait, "TSENATAnalysis")) {
        stop("'analysis_sait' must be a TSENATAnalysis object", call. = FALSE)
    }
    
    # Validate analysis_rank if provided
    if (!is.null(analysis_rank) && !is(analysis_rank, "TSENATAnalysis")) {
        stop("'analysis_rank' must be a TSENATAnalysis object", call. = FALSE)
    }
}

#' Helper: Handle two-object concordance API
#' @param analysis_sait TSENATAnalysis object
#' @param analysis_rank TSENATAnalysis object
#' @param verbose Logical; print progress
#' @return List with concordance_result, sait_method, rank_method
#' @noRd
.concordance_two_objects <- function(analysis_sait, analysis_rank, verbose) {
    if (verbose) {
        message("[calculate_concordance] Using two TSENATAnalysis objects")
    }
    
    concordance_result <- tryCatch({
        .calculate_concordance(analysis_sait = analysis_sait, analysis_rank = analysis_rank)
    }, error = function(e) {
        stop("[calculate_concordance] ", conditionMessage(e), call. = FALSE)
    })
    
    list(
        concordance_result = concordance_result,
        sait_method = concordance_result$sait_method,
        rank_method = concordance_result$rank_method
    )
}

#' Helper: Handle legacy single-object concordance API
#' @param analysis_sait TSENATAnalysis object
#' @param verbose Logical; print progress
#' @return List with concordance_result, sait_method, rank_method
#' @noRd
.concordance_legacy_api <- function(analysis_sait, verbose) {
    if (verbose) {
        message("[calculate_concordance] Using legacy single-object API")
    }
    
    # Validate SAIT results
    if (is.null(analysis_sait@sait_results) || length(analysis_sait@sait_results) == 0) {
        stop("No SAIT results found in analysis_sait@sait_results. Run calculate_sait() first.",
            call. = FALSE)
    }
    
    # Auto-detect SAIT method
    default_sait_method <- names(analysis_sait@sait_results)[1]
    
    # Auto-detect rank method
    rank_method <- "rank_test"
    if (is.null(analysis_sait@rank_test_results) || 
        length(analysis_sait@rank_test_results) == 0 ||
        !("rank_test" %in% names(analysis_sait@rank_test_results))) {
        if (is.null(analysis_sait@rank_test_results) || 
            length(analysis_sait@rank_test_results) == 0) {
            stop("No rank test results found. Run calculate_rank_transform() first.",
                call. = FALSE)
        }
        rank_method <- names(analysis_sait@rank_test_results)[1]
    }
    
    # Extract and validate results
    sait_results_final <- analysis_sait@sait_results[[default_sait_method]]
    rank_test_results <- analysis_sait@rank_test_results[[rank_method]]
    
    if (!is.data.frame(sait_results_final)) {
        stop("SAIT results ('", default_sait_method, "') must be a data.frame", 
             call. = FALSE)
    }
    
    if (!is.data.frame(rank_test_results)) {
        stop("Rank test results ('", rank_method, "') must be a data.frame", 
             call. = FALSE)
    }
    
    if (verbose) {
        message("[calculate_concordance] Computing concordance between ", 
                default_sait_method, " and ", rank_method)
    }
    
    # Create temporary analysis objects for the refactored function
    temp_sait <- analysis_sait
    temp_sait@sait_results <- list(temp = sait_results_final)
    temp_rank <- analysis_sait
    temp_rank@rank_test_results <- list(temp = rank_test_results)
    
    concordance_result <- tryCatch({
        .calculate_concordance(analysis_sait = temp_sait, analysis_rank = temp_rank)
    }, error = function(e) {
        stop("[calculate_concordance] ", conditionMessage(e), call. = FALSE)
    })
    
    list(
        concordance_result = concordance_result,
        sait_method = default_sait_method,
        rank_method = rank_method
    )
}

#' Helper: Store concordance results in metadata
#' @param analysis TSENATAnalysis object
#' @param concordance_result List from .calculate_concordance()
#' @param sait_method Character; SAIT method name
#' @param rank_method Character; rank method name
#' @param verbose Logical; print progress
#' @return TSENATAnalysis object with updated metadata
#' @noRd
.store_concordance_metadata <- function(analysis, concordance_result, 
                                        sait_method, rank_method, verbose) {
    # Store results in metadata
    analysis@metadata$method_concordance <- list(
        comparison_df = concordance_result$comparison_df,
        spearman_rho = concordance_result$spearman_rho,
        high_confidence = concordance_result$high_conf,
        agreement_table = concordance_result$agreement_table,
        sait_method = sait_method,
        rank_method = rank_method,
        timestamp = Sys.time()
    )
    
    # Track function call
    analysis@metadata$function_calls <- c(
        analysis@metadata$function_calls,
        sprintf("calculate_concordance[%s vs %s]", sait_method, rank_method)
    )
    
    if (verbose) {
        message("[calculate_concordance] Concordance computed successfully")
        if (!is.na(concordance_result$spearman_rho)) {
            message("[calculate_concordance] Spearman correlation = ", 
                    round(concordance_result$spearman_rho, 3))
        }
    }
    
    analysis
}

#' Helper: Write concordance results to file
#' @param analysis TSENATAnalysis object
#' @param output_file Character; file path
#' @param verbose Logical; print progress
#' @return TSENATAnalysis object with updated metadata
#' @noRd
.write_concordance_file <- function(analysis, output_file, verbose) {
    if (is.null(output_file)) {
        return(analysis)
    }
    
    # Ensure .txt extension
    if (!grepl("\\.txt$", output_file, ignore.case = TRUE)) {
        output_file <- paste0(output_file, ".txt")
    }
    
    # Generate and write results
    concordance_text <- results(analysis, type = "concordance")
    writeLines(concordance_text, con = output_file)
    
    if (verbose) {
        message("[calculate_concordance] Results written to: ", output_file)
    }
    
    # Store path in metadata
    analysis@metadata$concordance_results_file <- output_file
    
    analysis
}

setMethod("calculate_concordance", "TSENATAnalysis", function(analysis_sait, 
    analysis_rank = NULL, verbose = FALSE, output_file = NULL, ...) {
    
    # Validate inputs
    .validate_concordance_inputs(analysis_sait, analysis_rank, ...)
    
    # Route to appropriate API
    if (!is.null(analysis_rank)) {
        result_list <- .concordance_two_objects(analysis_sait, analysis_rank, verbose)
    } else {
        result_list <- .concordance_legacy_api(analysis_sait, verbose)
    }
    
    # Store metadata
    analysis_sait <- .store_concordance_metadata(
        analysis_sait, 
        result_list$concordance_result,
        result_list$sait_method, 
        result_list$rank_method,
        verbose
    )
    
    # Write output file if specified
    analysis_sait <- .write_concordance_file(analysis_sait, output_file, verbose)
    
    analysis_sait
})

