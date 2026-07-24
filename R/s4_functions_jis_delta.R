#' Plot Multi-Q Delta Influence Heatmaps from TSENATAnalysis Object
#'
#' S4 wrapper for  \code{. plot_jis_delta()} that 
#' extracts results
#' directly from a TSENATAnalysis object. Automatically retrieves jackknife
#' switching
#' results from the analysis object slots.
#'
#' @param analysis \code{TSENATAnalysis}. An S4 object containing completed
#'   jackknife isoform switching analysis across multiple q-values.
#'
#' @param n_genes \code{numeric}. Number of top genes to display in heatmaps
#'   (default: 4). Genes are ranked by LM p-values if available, otherwise
#'   by order of appearance in results.
#'
#' @param sait_results \code{data.frame} or \code{NULL}. Optional SAIT interaction
#' results for ranking genes (default: NULL). If NULL, attempts to extract
#' from
#'   \code{analysis@sait_results$sait_interaction}.
#'
#' @param verbose \code{logical}. If \code{TRUE}, print diagnostic messages
#'   during plot generation (default: FALSE).
#'
#' @param output_file \code{character} or  \code{NULL}.
#'  Optional file path to save the plot.
#'   Default: NULL (no file output).
#'
#' @param ... Additional arguments passed to the base function.
#'
#' @return A file path (character) to the saved heatmap PNG file, invisibly.
#'
#' @details
#' This function extracts the following from \code{analysis}:
#' \describe{
#'   \item{Jackknife results}{From \code{analysis@jackknife_results},  which 
#' should
#'         contain multi-q switching results keyed by q-value (e.g., 'q_1.00')}
#'   \item{SAIT results}{From \code{analysis@sait_results$sait_interaction} if not
#'         explicitly provided, for ranking genes by significance}
#' }
#'
#' The wrapper automatically handles parameter extraction and provides a
#' simplified
#' interface compared to the base function.
#'
#' @examples
#' # Plot 5: Multi-q delta influence (isoform switching) heatmaps
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
#' # Configure analysis parameters first
#' config <- TSENAT_config(
#'   sample_col = 'sample',
#'   condition_col = 'condition',
#'   subject_col = 'paired_samples',
#'   paired = TRUE,
#'   control = 'normal',
#'   q = seq(0, 2, by = 0.1)
#' )
#'
#' # Build analysis with configured parameters
#' analysis <- build_analysis(
#'   readcounts = readcounts,
#'   tx2gene = gff3_dataset,
#'   metadata = metadata_df,
#'   config = config,
#'   tpm = tpm,
#'   effective_length = effective_length
#' )
#'
#' analysis <- filter_analysis(analysis, stringency = 'severe')
#' analysis <- calculate_diversity(
#'   analysis,
#'   q = seq(0.2, 2, by = 0.4),
#'   verbose = FALSE
#' )
#' analysis <- calculate_divergence(
#'   analysis,
#'   q = seq(0.2, 2, by = 0.4)
#' )
#' analysis <- suppressWarnings(calculate_sait(analysis, method = 'lmm'))
#' analysis <- calculate_jis(
#'   analysis,
#'   q = seq(0.2, 2, by = 0.4),
#'   n_bootstrap = 20
#' )
#' heatmap_file <- plot_jis_delta(analysis, n_genes = 2)
#'
#' @seealso
#' \code{\link{calculate_jis}} for computing switching results
#'
#' @export
#' @importFrom methods is
plot_jis_delta <- function(analysis, n_genes = 4, sait_results = NULL, verbose = FALSE,
    output_file = NULL, ...) {

    # Load visualization dependencies (ggplot2, cowplot, pheatmap, etc.)
    .load_visualization_deps()

    # Extract verbose parameter if not provided
    verbose <- resolve_slot_param(verbose, analysis@config, "verbose", FALSE)

    # Validation
    if (!is(analysis, "TSENATAnalysis")) {
        stop("analysis must be a TSENATAnalysis object", call. = FALSE)
    }

    if (verbose)
        message("Extracting jackknife switching results from analysis object...")

    # Extract jackknife/switching results
    jackknife_results_list <- analysis@jackknife_results
    if (is.null(jackknife_results_list) || length(jackknife_results_list) == 0) {
        stop("No jackknife results found in analysis@jackknife_results. ", "Run calculate_jis() first.",
            call. = FALSE)
    }

    # Check for multi-q result (stored under 'multi_q' key when multiple
    # q-values provided)
    if ("multi_q" %in% names(jackknife_results_list)) {
        switching_results <- jackknife_results_list$multi_q
        if (verbose) {
            message("  [OK] Found multi-q result with class: ", class(switching_results)[1])
        }
    } else {
        # Fallback: use all results as list (for single or multiple q-values)
        switching_results <- jackknife_results_list
        if (verbose) {
            message("  [OK] Using individual q-value results (", length(switching_results),
                " q-values)")
        }
    }

    if (verbose) {
        message("  Q-values: ", paste(names(switching_results), collapse = ", "))
    }

    # Extract SAIT results if not provided
    if (is.null(sait_results)) {
        if (verbose)
            message("Extracting SAIT results from analysis@sait_results...")

        sait_results_list <- analysis@sait_results
        if (!is.null(sait_results_list)) {
            if (!is.null(sait_results_list$sait_interaction)) {
                if (is.data.frame(sait_results_list$sait_interaction$results)) {
                  sait_results <- sait_results_list$sait_interaction$results
                } else if (is.data.frame(sait_results_list$sait_interaction)) {
                  sait_results <- sait_results_list$sait_interaction
                }
            }

            if (!is.null(sait_results)) {
                if (verbose)
                  message("  [OK] Extracted SAIT results with ", nrow(sait_results),
                    " genes")
            } else if (verbose) {
                message("  SAIT results not found; genes will be ranked by appearance")
            }
        }
    }

    if (verbose)
        message("Calling .plot_jis_delta()...")

    # Call base function with extracted parameters Note: output_file parameter
    # can be used to save heatmap as PNG file
    result <- .plot_jis_delta(switching_results = switching_results, n_genes = n_genes,
        sait_results = sait_results, verbose = verbose, output_file = output_file, ...)

    if (verbose) {
        message("[OK] Heatmap plot generated successfully")
    }

    invisible(result)
}

