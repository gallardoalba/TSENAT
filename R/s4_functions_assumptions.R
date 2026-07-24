#' Test statistical assumptions on diversity data in TSENATAnalysis
#'
#' @param analysis \code{TSENATAnalysis} object with diversity results stored
#'   in \code{@diversity_results}.
#' @param q \code{numeric}. Q-value(s) to extract from diversity results.
#'   If NULL, uses the first available diversity result or q=1.0.
#' @param checks \code{character}. Which assumptions to test (default: 'rank').
#'   Presets:
#'   - 'rank': core assumption checks (exchangeability, monotonicity, consistency)
#'   - 'all': all checks including GAM diagnostics
#'   Explicit: character vector like \code{c('exchangeability', 'monotonicity')}.
#' @param alpha \code{numeric}. Significance level for tests (default: 0.05).
#' @param format \code{character}. Output format when used with \code{results()}.
#'   'text' (default): formatted text output for display
#'   'list': returns structured list for programmatic access.
#' @param ... Additional arguments (output_file, verbose for file output).
#'
#' @return Modified TSENATAnalysis object with assumption test results stored
#'   in \code{@metadata$rankbased_assumptions}.
#'
#' @details
#' This wrapper calls \code{.calculate_assumptions()} on diversity data
#' extracted from the analysis object. Evaluates data stability and 
#' consistency across dimensions. Results include:
#'
#' \describe{
#'   \item{exchangeability}{Permutation test for independence and temporal/spatial structure}
#'   \item{monotonicity}{Spearman correlation consistency across rows}
#'   \item{consistency}{Kendall's W concordance and ICC across samples}
#' }
#'
#' **Data Extraction Priority:**
#' 1. If q specified: uses diversity result for that q-value
#' 2. If q NULL: uses first available diversity result
#' 3.  If no diversity results:
#'  extracts from cached combined result (\code{@metadata$diversity_combined})
#'
#' @examples
#' # Load example data (matching TSENAT.Rmd workflow)
#' data(readcounts)
#' readcounts <- as.matrix(readcounts)
#' mode(readcounts) <- 'numeric'
#' metadata_df <- read.table(
#'   system.file('extdata', 'metadata.tsv', package = 'TSENAT'),
#'   header = TRUE, sep = '\t'
#' )
#' gff3_dataset <- system.file('extdata', 'annotation.gff3.gz', package =
#' 'TSENAT')
#' 
#' # Build analysis from vignette data and create small subset
#' config <- TSENAT_config(
#'   sample_col = 'sample',
#'   condition_col = 'condition',
#'   q_values = seq(0, 2, by = 0.05),
#'   paired = FALSE
#' )
#' analysis <- build_analysis(
#'   readcounts = readcounts,
#'   metadata = metadata_df,
#'   tx2gene = gff3_dataset,
#'   config = config,
#'   tpm = tpm,
#'   effective_length = effective_length
#' )
#' analysis <- filter_analysis(analysis, min_samples = 1, subset_n_genes = 200)
#' analysis <- calculate_diversity(analysis, q = c(0.5, 1.0, 1.5))
#' analysis <- calculate_assumptions(analysis, q = 1.0)
#' # Check results using rank_test accessor
#' results_df <- results(analysis, type = 'rank_test')
#'
#' @export
#' @rdname calculate_assumptions
setGeneric("calculate_assumptions", function(analysis, q = NULL, checks = "rank",
    alpha = 0.05, format = "text", ...) {
    standardGeneric("calculate_assumptions")
})

#' @rdname calculate_assumptions
setMethod("calculate_assumptions", signature(analysis = "TSENATAnalysis"), function(analysis,
    q = NULL, checks = "rank", alpha = 0.05, format = "text", ...) {

    # Validate inputs
    .validate_object_class(analysis, "TSENATAnalysis", "analysis")

    # Extract diversity data using consolidated helper
    diversity_result <- .extract_diversity_data(analysis, q = q)
    diversity_data <- diversity_result$diversity_data
    q_used <- diversity_result$q_used

    if (is.null(diversity_data)) {
        stop("No diversity results found in analysis object. ",
            "Run calculate_diversity() first.", call. = FALSE)
    }

    # Extract q-values from diversity_results names for GAM metrics
    q_values <- NULL
    if (q_used == "all" && length(analysis@diversity_results) > 1) {
        q_names <- names(analysis@diversity_results)
        q_values <- as.numeric(gsub("q_", "", q_names))
    }

    # Run assumptions test
    result <- tryCatch({
        .calculate_assumptions(data = diversity_data, checks = checks, alpha = alpha,
            q_values = q_values)
    }, error = function(e) {
        stop("Assumptions test failed:\n", e$message, call. = FALSE)
    })

    # Store results
    analysis@metadata$rankbased_assumptions <- list(result = result, q_value_tested = q_used,
        checks_performed = checks, alpha_used = alpha, timestamp = Sys.time())

    # Track function call
    analysis@metadata$function_calls <- c(analysis@metadata$function_calls, paste0("calculate_assumptions[q=",
        q_used, "]"))

    # Handle optional parameters and output file
    dots <- list(...)
    params <- .extract_dot_params(dots, c("output_file", "verbose"), 
                                   list(output_file = NULL, verbose = FALSE))

    if (!is.null(params$output_file)) {
        assumptions_df <- .format_assumptions_for_output(result)
        .handle_optional_output(analysis, params$output_file, assumptions_df,
            "calculate_assumptions", params$verbose)
    }

    analysis
})

# Helper function to extract q-value from key

.extract_q_from_key <- function(key) {
    # Extract numeric part from 'q_X.X' format
    as.numeric(sub("^q_", "", key))
}

# Helper: Format assumptions results for output
#' @noRd
.format_assumptions_for_output <- function(result) {
    # Use the new .process_assumptions_results formatter from
    # orchestration_results Convert structured assumptions to formatted data
    # frame suitable for TSV output

    if (!is.list(result)) {
        return(data.frame(check = "assumptions", status = "error", details = "Invalid result format"))
    }

    # Process using the orchestration_results formatter with format='list'
    processed <- .process_assumptions_results(result, format = "list")

    if (is.null(processed) || is.null(processed$assumptions_table)) {
        return(data.frame(check = "assumptions", status = "error", details = "No checks found"))
    }

    # Get the assumptions table
    assumptions_df <- processed$assumptions_table

    # Ensure it has the expected columns
    if (!all(c("Test", "Result", "Interpretation") %in% colnames(assumptions_df))) {
        return(data.frame(check = "assumptions", status = "error", details = "Invalid table structure"))
    }

    # Format for TSV output - keep original column names but convert to row
    # format for readability
    output_df <- data.frame(check = "assumption", test = assumptions_df$Test, 
        result = assumptions_df$Result, interpretation = assumptions_df$Interpretation,
        stringsAsFactors = FALSE)

    output_df
}

# ============================================================================
# OPTIMIZATION: Consolidated Helper Functions for Refactoring
# ============================================================================
# These helpers consolidate repeated code patterns to reduce duplication,
# improve testability, and increase coverage.

#' Extract parameters from dots list with defaults
#'
#' @param dots List from ... arguments
#' @param param_names Character vector of parameter names to extract
#' @param defaults List of default values (optional)
#' @return Named list with extracted parameters
#' @noRd
.extract_dot_params <- function(dots, param_names, defaults = NULL) {
    result <- defaults %||% list()
    
    for (param in param_names) {
        if (param %in% names(dots)) {
            result[[param]] <- dots[[param]]
        }
    }
    
    result
}

#' Handle optional output file writing
#'
#' @param analysis TSENATAnalysis object
#' @param output_file Character path or NULL
#' @param results_df Data frame to write
#' @param func_name Character; function name for logging
#' @param verbose Logical; print progress
#' @return Invisible NULL (writes file as side effect)
#' @noRd
.handle_optional_output <- function(analysis, output_file, results_df, func_name, verbose) {
    if (is.null(output_file)) {
        return(invisible(NULL))
    }
    
    tryCatch({
        if (verbose) {
            message("[", func_name, "] Writing results to: ", output_file)
        }
        save_analysis_output(results_df, output_file, object = analysis,
            verbose = verbose, func_name = func_name)
    }, error = function(e) {
        warning("[", func_name, "] Could not write results to file: ",
            conditionMessage(e), call. = FALSE)
    })
    
    invisible(NULL)
}

#' Extract diversity data with multi-path fallback
#'
#' Handles diversity data extraction with smart fallback:
#' 1. If q specified: uses diversity result for that q-value
#' 2. If q NULL and multiple results: combines all q-values
#' 3. If q NULL and single result: uses it
#' 4. Fallback: uses first available
#'
#' @param analysis TSENATAnalysis object
#' @param q Numeric; specific q-value or NULL
#' @return List(diversity_data = matrix, q_used = character/numeric)
#' @noRd
.extract_diversity_data <- function(analysis, q = NULL) {
    diversity_data <- NULL
    q_used <- q
    
    # If q is specified, try to get that specific q-value
    if (!is.null(q)) {
        q_key <- paste0("q_", q)
        if (q_key %in% names(analysis@diversity_results)) {
            div_se <- analysis@diversity_results[[q_key]]
            diversity_data <- assay(div_se, "diversity")
        }
    }
    
    # If q is NULL and multiple diversity results exist, combine all q-values
    if (is.null(diversity_data) && is.null(q) && length(analysis@diversity_results) > 1) {
        entropy_list <- lapply(analysis@diversity_results, function(se) {
            mat <- assay(se, "diversity")
            if (!is.matrix(mat)) mat <- as.matrix(mat)
            return(mat)
        })
        
        all_genes <- lapply(entropy_list, rownames)
        common_genes <- Reduce(intersect, all_genes)
        entropy_list <- lapply(entropy_list, function(mat) {
            mat[common_genes, , drop = FALSE]
        })
        
        diversity_data <- do.call(cbind, entropy_list)
        q_used <- "all"
    }
    
    # If q is NULL and only one result, use it
    if (is.null(diversity_data) && is.null(q) && length(analysis@diversity_results) == 1) {
        div_se <- analysis@diversity_results[[1]]
        diversity_data <- assay(div_se, "diversity")
        q_used <- .extract_q_from_key(names(analysis@diversity_results)[1])
    }
    
    # Fallback: use first diversity result
    if (is.null(diversity_data) && length(analysis@diversity_results) > 0) {
        div_se <- analysis@diversity_results[[1]]
        diversity_data <- assay(div_se, "diversity")
        if (is.null(q_used) || is.na(q_used)) {
            q_used <- .extract_q_from_key(names(analysis@diversity_results)[1])
        }
    }
    
    if (is.null(diversity_data)) {
        return(list(diversity_data = NULL, q_used = NULL))
    }
    
    if (!is.matrix(diversity_data)) {
        diversity_data <- as.matrix(diversity_data)
    }
    
    list(diversity_data = diversity_data, q_used = q_used)
}

#' Validate object class and provide clear error messages
#'
#' @param obj Object to validate
#' @param expected_class Character; expected class name
#' @param param_name Character; parameter name for error message
#' @return Invisible(TRUE) on success, stops on failure
#' @noRd
.validate_object_class <- function(obj, expected_class, param_name = "object") {
    if (!is(obj, expected_class)) {
        stop("'", param_name, "' must be a ", expected_class, " object",
            call. = FALSE)
    }
    invisible(TRUE)
}


