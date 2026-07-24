#' M-Estimation for Sample Quality (S4 Wrapper)
#'
#' S4 wrapper for \code{m_estimate} that performs robust M-estimation
#' on diversity results stored in a TSENATAnalysis object and stores results
#' back into the object.
#'
#' @param analysis \code{TSENATAnalysis} object with diversity results
#'   (typically via \code{\link{calculate_diversity}}).
#' @param condition_col \code{character}.
#'  Column name in sample metadata indicating
#'   condition/sample grouping. Auto-detected from \code{@config$condition_col}
#'   if available.
#' @param loss_type \code{character}. Type of loss function: 'huber' (default),
#'   'tukey', or 'lsq'. Determines robustness vs efficiency trade-off.
#' @param scale \code{numeric}.  Manual scale parameter.  If NULL,
#'  estimated from data.
#' @param max_iter \code{integer}.  Maximum iterations for  M-estimation.
#'  Default:  50.
#' @param tol \code{numeric}. Convergence tolerance. Default: 1e-6.
#' @param paired \code{logical}. If TRUE, adjusts degrees of freedom for paired
#'   designs. Auto-detected from \code{@config$paired} if available.
#'   Default: FALSE.
#' @param pcorr \code{character}.  P-value correction method.  Default:
#'  'BH' (Benjamini-Hochberg).
#' @param q_combine_method \code{character}. How to collapse multi-q results:
#'   'mean' (default) or 'median'.
#' @param influence_threshold \code{numeric}. Threshold for classifying samples
#'   as high-influence. Default: 0.75.
#' @param scale_method \code{character}.  Scale estimation method:
#'  'mad' (default),
#'   'proposal2', or 's-estimator'.
#' @param output_file \code{character} or  \code{NULL}.
#'  Optional file path to save results.
#' Supported formats: .rds (for S4 objects), .tsv, .csv, .txt (for tables).
#' Default: NULL (no file output).
#' @param verbose \code{logical}. Print status messages. Default: TRUE.
#'
#' @return
#' Modified TSENATAnalysis object with M-estimation results stored in
#' \code{analysis@metadata$m_estimate_results}. Contains data frame with
#' influence scores, robustness weights, entropy statistics, and QC classifications.
#' Returns visibly to support method chaining and piping.
#'
#' @details
#' This wrapper extracts diversity results from \code{analysis@diversity_results},
#' performs robust M-estimation on diversity values (entropy), and stores
#' results
#' in the analysis object metadata.
#'
#' **M-Estimation:** Robust regression technique that down-weights outliers
#' based
#' on their residuals. Useful for detecting low-quality samples that show
#' unusual diversity patterns.
#'
#' **M-Estimation Results include:**
#' \itemize{
#'   \item \code{sample_influence}: How much each sample affects the overall fit
#'   \item \code{robustness_weight}: Down-weighting factor
#'   \item \code{entropy_mean}: Average entropy for the sample
#'   \item \code{entropy_sd}: Entropy variability within the sample
#'   \item \code{Status}:  QC Classification 
#' }
#'
#' **Parameter resolution priority** (explicit > @config > error):
#' \itemize{
#'   \item \code{samples}: Uses explicit arg, else \code{@config$condition_col},
#' else error. Note: despite parameter name 'samples', maps to condition
#' grouping column
#' }
#'
#' **Data Requirements:**
#' \itemize{
#'   \item Diversity results must be computed via \code{calculate_diversity()}
#'   \item Sample grouping column required in colData (auto-detected from @config$condition_col
#'     or via 'samples' parameter)
#' }
#'
#' @examples
#' # Create test analysis and compute M-estimation
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
#' # Create config (metadata passed as explicit parameter to build_analysis)
#' config <- TSENAT_config(
#'   sample_col = 'sample',
#'   condition_col = 'condition',
#'   q_values = seq(0, 2, by = 0.05),
#'   paired = FALSE
#' )
#' 
#' # Build analysis from vignette data
#' analysis <- build_analysis(
#'   readcounts = readcounts,
#'   tx2gene = gff3_dataset,
#'   metadata = metadata_df,
#'   config = config,
#'   tpm = tpm,
#'   effective_length = effective_length
#' )
#' analysis <- filter_analysis(analysis, min_samples = 1, subset_n_genes = 200)
#' analysis <- calculate_diversity(analysis, q = c(0.5, 1.0, 1.5))
#' analysis <- calculate_m_estimator(
#'   analysis,
#'   condition_col = 'condition',
#'   loss_type = 'huber'
#' )
#'
#' @seealso
#' \code{\link{calculate_diversity}} for computing diversity
#'
#' @export
#' @importFrom utils write.table
calculate_m_estimator <- function(analysis, condition_col = NULL, loss_type = "huber",
    scale = NULL, max_iter = 50, tol = 1e-06, paired = NULL, pcorr = "BH", q_combine_method = "mean",
    influence_threshold = 0.75, scale_method = "mad", output_file = NULL, verbose = FALSE) {

    # Validate input
    if (!is(analysis, "TSENATAnalysis")) {
        stop("'analysis' must be a TSENATAnalysis object", call. = FALSE)
    }

    # Check for diversity results
    if (is.null(analysis@diversity_results) || length(analysis@diversity_results) ==
        0) {
        stop("Diversity results not found. Run calculate_diversity() first.", call. = FALSE)
    }

    # Validate that diversity_results is a properly structured named list
    if (!is.list(analysis@diversity_results) || is.null(names(analysis@diversity_results))) {
        stop("Diversity results must be a named list of SummarizedExperiment objects",
            call. = FALSE)
    }

    # Auto-detect condition_col if not provided
    if (is.null(condition_col)) {
        if ("condition_col" %in% names(analysis@config)) {
            condition_col <- analysis@config$condition_col
            if (is.null(condition_col) || !is.character(condition_col) || condition_col ==
                "") {
                stop("@config$condition_col must be a non-empty character value",
                  call. = FALSE)
            }
            if (verbose) {
                message(sprintf("Auto-detected 'condition_col' from config: %s",
                  condition_col))
            }
        } else {
            cd_cols <- colnames(SummarizedExperiment::colData(analysis@diversity_results[[1]]))
            stop("Sample grouping column not specified:\n", "  Available colData columns: ",
                paste(cd_cols, collapse = ", "), "\n\n", "SOLUTION: Set @config$condition_col or pass 'condition_col' parameter\n",
                "  Example: analysis@config$condition_col <- 'sample_type'\n", "  Or:      calculate_m_estimator(analysis, condition_col = 'sample_type')\n",
                call. = FALSE)
        }
    } else if (!is.character(condition_col) || length(condition_col) != 1) {
        stop("'condition_col' must be a single character value", call. = FALSE)
    }

    # Auto-detect paired parameter from @config if not explicitly provided
    paired <- resolve_slot_param(paired, analysis@config, "paired", FALSE)

    if (!is.logical(paired) || length(paired) != 1) {
        stop("'paired' must be a single logical value (TRUE or FALSE)", call. = FALSE)
    }

    # Extract diversity results - get first SE to access sample metadata
    diversity_se <- analysis@diversity_results[[1]]

    if (is.null(diversity_se) || nrow(diversity_se) == 0) {
        stop("Diversity SummarizedExperiment is empty", call. = FALSE)
    }

    # Verify condition_col exists
    sample_info <- SummarizedExperiment::colData(diversity_se)
    if (!(condition_col %in% colnames(sample_info))) {
        stop("Column '", condition_col, "' not found in sample metadata.\n", "Available columns: ",
            paste(colnames(sample_info), collapse = ", "), "\n\n", "SOLUTION: Use a valid column name\n",
            "  Example: calculate_m_estimator(analysis, condition_col = 'sample_type')\n",
            call. = FALSE)
    }

    # Combine all q-value diversity results into a single matrix (m_estimate
    # needs all diversity data in one SE)
    if (verbose) {
        message(sprintf("Combining %d q-value diversity results...", length(analysis@diversity_results)))
    }

    first_se <- analysis@diversity_results[[1]]
    combined_assay <- SummarizedExperiment::assay(first_se)
    combined_colnames <- colnames(first_se)

    # Add other q-values
    for (q_name in names(analysis@diversity_results)[-1]) {
        se_q <- analysis@diversity_results[[q_name]]
        combined_assay <- cbind(combined_assay, SummarizedExperiment::assay(se_q))
        combined_colnames <- c(combined_colnames, colnames(se_q))
    }

    # Update column names to reflect combined data
    colnames(combined_assay) <- combined_colnames
    single_colData <- SummarizedExperiment::colData(first_se)

    # Replicate colData for each q-value
    n_q_values <- length(analysis@diversity_results)
    combined_colData <- do.call(rbind, replicate(n_q_values, single_colData, simplify = FALSE))
    rownames(combined_colData) <- combined_colnames

    # Create combined SummarizedExperiment
    combined_se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = combined_assay),
        colData = combined_colData)

    # Run m_estimate
    if (verbose) {
        message("Running M-estimation on combined diversity...")
    }

    m_est_results <- tryCatch({
        result <- .calculate_m_estimator(x = combined_se, samples = condition_col,
            loss_type = loss_type, scale = scale, max_iter = max_iter, tol = tol,
            paired = paired, pcorr = pcorr, q_combine_method = q_combine_method,
            influence_threshold = influence_threshold, scale_method = scale_method)
        result
    }, error = function(e) {
        stop("M-estimation failed:\n", e$message, call. = FALSE)
    })

    # Store results in metadata
    analysis@metadata$m_estimate_results <- m_est_results

    # Save to output file if provided
    if (!is.null(output_file)) {
        if (!is.character(output_file) || length(output_file) != 1) {
            stop("'output_file' must be a character string (file path)", call. = FALSE)
        }
        # Use save_analysis_output for consistent handling (data frame for TSV,
        # RDS for S4)
        save_analysis_output(m_est_results, output_file, object = analysis, verbose = verbose,
            func_name = "calculate_m_estimator")
    }

    # Track function call
    analysis@metadata$function_calls <- c(analysis@metadata$function_calls, paste0("calculate_m_estimator[condition_col=",
        condition_col, ",loss_type=", loss_type, "]"))

    if (verbose) {
        message("M-estimation complete. Results stored in @metadata$m_estimate_results")
    }

    invisible(analysis)
}

