#' Plot Global Divergence q-Curve Across All Genes (S4 Wrapper)
#'
#' S4 wrapper that extracts divergence results from a TSENATAnalysis object
#' and visualizes the average Tsallis divergence across all genes (or specified
#' genes) as a function of q-value.
#'
#' @param analysis \code{TSENATAnalysis} object with divergence results
#'   (typically via \code{\link{calculate_divergence}}).
#' @param gene \code{character}. Optional specific gene name to plot.
#'   If NULL, plots global divergence curve (aggregated across all genes).
#' @param n_genes \code{integer}. Number of top genes to plot when showing
#'   multi-gene spectra. Default is 4. Genes are sorted by p-value significance.
#' @param ncol \code{integer}. Number of columns in grid layout for multi-gene
#'   plots. Default is 2. Number of rows is automatically calculated.
#' @param metric \code{character}. Summary statistic for global curve:
#'   'median' (default) or 'mean'. Only used when gene = NULL.
#' @param variability_metric \code{character}. Error bar type for global curve:
#'   'iqr' (default) or 'sd'. Only used when gene = NULL.
#' @param use_pvalue_ranking \code{logical}.  If TRUE,
#'  uses SAIT results to rank and
#' display top n_genes by p-value significance. If FALSE (default), plots
#' global
#'   divergence curve when gene = NULL. Default is FALSE.
#' @param output_file \code{character}. Optional file path to save the plot.
#'   If NULL, plot is returned but not saved.
#' @param width \code{numeric}. Plot width in inches. Default is 10.
#' @param height \code{numeric}. Plot height in inches. Default is 6.
#' @param verbose \code{logical}. Print status messages. Default is TRUE.
#' @param ... Additional arguments passed to the underlying plotting function.
#'
#' @return
#' Invisibly returns the file path if saved, otherwise the ggplot object.
#' If the plot cannot be created (missing data, ggplot2 not available),
#' returns NULL invisibly with an informative message.
#'
#' @details
#' This wrapper extracts the divergence SummarizedExperiment from
#' \code{analysis@divergence_results} and optionally the SAIT results from
#' \code{analysis@sait_results$sait_interaction} to pass to the base function.
#'
#' **Data Requirements:**
#' \itemize{
#'   \item Divergence must be computed via \code{calculate_divergence()}
#'   \item \code{@divergence_results$divergence_se} or direct divergence SE
#' }
#'
#' **Modes:**
#' \itemize{
#'   \item \strong{Global mode} (gene = NULL): Shows median/mean divergence
#'         across all genes with variability bands
#'   \item \strong{Gene-specific mode} (gene specified): Shows divergence
#'         spectrum for a single named gene
#'   \item \strong{Top genes mode} (gene = NULL, sait_res provided): Shows
#'         top n_genes by significance
#' }
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
#' config <- TSENAT_config(sample_col = 'sample', condition_col = 'condition')
#' analysis <- build_analysis(readcounts = readcounts, tx2gene =
#' gff3_dataset, metadata = metadata_df, config = config,
#'   tpm = tpm, effective_length = effective_length)
#' analysis <- filter_analysis(
#'   analysis,
#'   min_samples = 1,
#'   subset_n_genes = 200
#' )
#' analysis <- calculate_diversity(
#'   analysis,
#'   q = c(0.5, 1, 1.5),
#'   verbose = FALSE
#' )
#' analysis <- calculate_divergence(
#'   analysis,
#'   q = c(0.5, 1, 1.5),
#'   verbose = FALSE
#' )
#' p_global <- plot_divergence_spectrum(analysis)
#' # print(p_global)
#'
#' @seealso
#' \code{\link{calculate_divergence}} for computing divergence.
#'
#' @export
plot_divergence_spectrum <- function(analysis, gene = NULL, n_genes = 4, ncol = 2,
    metric = c("median", "mean"), variability_metric = c("iqr", "sd"), use_pvalue_ranking = FALSE,
    output_file = NULL, width = 12, height = NULL, verbose = FALSE, ...) {

    # Load visualization dependencies (ggplot2, cowplot, etc.)
    .load_visualization_deps()

    # Validate input FIRST (before accessing analysis@config)
    if (!is(analysis, "TSENATAnalysis")) {
        stop("'analysis' must be a TSENATAnalysis object", call. = FALSE)
    }

    # Extract verbose parameter if not provided
    verbose <- resolve_slot_param(verbose, analysis@config, "verbose", TRUE)

    # Match metric and variability_metric arguments
    metric <- match.arg(metric)
    variability_metric <- match.arg(variability_metric)

    # Extract divergence SE from analysis object
    if (is.null(analysis@divergence_results)) {
        stop("Divergence results not found in analysis@divergence_results. ", "Run calculate_divergence() first.",
            call. = FALSE)
    }

    # Handle both direct SE and wrapped 'divergence_se' key
    divergence_results_se <- if (is.list(analysis@divergence_results) && "divergence_se" %in%
        names(analysis@divergence_results)) {
        analysis@divergence_results$divergence_se
    } else if (is(analysis@divergence_results, "SummarizedExperiment")) {
        analysis@divergence_results
    } else {
        stop("Invalid divergence_results structure. Expected SummarizedExperiment or list with 'divergence_se' key",
            call. = FALSE)
    }

    if (nrow(divergence_results_se) == 0 || ncol(divergence_results_se) == 0) {
        stop("Divergence SummarizedExperiment is empty", call. = FALSE)
    }

    # Extract SAIT results for sait_res parameter (optional) Only use for ranking
    # if use_pvalue_ranking = TRUE
    sait_res <- NULL
    if (use_pvalue_ranking && !is.null(analysis@sait_results) && is.list(analysis@sait_results)) {
        if ("sait_interaction" %in% names(analysis@sait_results)) {
            sait_res <- analysis@sait_results$sait_interaction
        } else if (length(analysis@sait_results) > 0) {
            sait_res <- analysis@sait_results[[1]]
        }
    }

    # Validate SAIT results if using multi-gene mode
    if (is.null(gene) && use_pvalue_ranking && !is.null(sait_res)) {
        if (!is.data.frame(sait_res) || nrow(sait_res) == 0) {
            if (verbose) {
                message("Note: Invalid SAIT results. Plotting global curve without gene ranking.")
            }
            sait_res <- NULL
        }
    }

    # Calculate height if not provided (based on grid layout)
    if (is.null(height)) {
        n_rows <- ceiling(n_genes/ncol)
        height <- 3 + (2.5 * n_rows)  # 3" base + 2.5" per row
    }

    # Create the plot using base function
    p <- tryCatch({
        .plot_divergence_spectrum(divergence_results_se = divergence_results_se,
            gene = gene, sait_res = sait_res, n_genes = n_genes, ncol = ncol, metric = metric,
            variability_metric = variability_metric, analysis = analysis, ...)
    }, error = function(e) {
        if (verbose) {
            message("plot_divergence_spectrum failed: ", e$message)
        }
        return(NULL)
    })

    # If plot creation failed, return NULL invisibly
    if (is.null(p)) {
        return(invisible(NULL))
    }

    # Save plot to file if requested
    if (!is.null(output_file)) {
        save_analysis_output(p, output_file, object = analysis, verbose = verbose,
            func_name = "plot_divergence_spectrum", width = width, height = height)
        # Documented contract: return the file path invisibly when saved
        return(invisible(output_file))
    }

    # Return plot object when no file is requested
    p
}

# ============================================================================
# PLOT WRAPPER - Plot Method Concordance Comparison
# ============================================================================

