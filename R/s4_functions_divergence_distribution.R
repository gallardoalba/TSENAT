#' Plot Tsallis Divergence Effect Size Distribution (S4 Wrapper)
#'
#' S4 wrapper that extracts effect size results from a TSENATAnalysis object
#' and generates a histogram visualization of Tsallis divergence effect sizes
#' across genes.
#'
#' @param analysis \code{TSENATAnalysis} object with effect sizes computed
#'   (typically via \code{\link{calculate_effect_sizes}}).
#' @param threshold \code{numeric}. Effect size threshold for visual marking
#'   in the plot. Default is 0.1 (information-theoretic significance level).
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
#' This wrapper extracts the interaction results (with effect size columns)
#' from \code{analysis@metadata$effect_sizes_divergence$interaction_results}
#' and passes them to the base \code{.plot_divergence_distribution()} function.
#'
#' The function visualizes the distribution of effect sizes using the median
#' q-value's divergence (typically around q=1.0, close to Shannon entropy).
#' A red dashed line marks the information-theoretic significance threshold.
#'
#' **Data Requirements:**
#' \itemize{
#'   \item Effect sizes must be computed via \code{calculate_effect_sizes()}
#'   \item \code{@metadata$effect_sizes_divergence$interaction_results} must
#'         contain columns matching pattern \code{effect_size_D_q*}
#' }
#'
#' @examples
#' # Plot 2: Distribution of effect sizes across genes
#' data(readcounts)
#' readcounts <- as.matrix(readcounts)
#' mode(readcounts) <- 'numeric'
#' metadata_df <- read.table(
#'   system.file('extdata', 'metadata.tsv', package = 'TSENAT'),
#'   header = TRUE, sep = '\t'
#' )
#' gff3_dataset <- system.file('extdata', 'annotation.gff3.gz', package =
#' 'TSENAT')
#' # Configure analysis parameters first
#' config <- TSENAT_config(
#'   sample_col = 'sample',
#'   condition_col = 'condition',
#'   control = 'normal'
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
#'   q = seq(0.2, 2, by = 0.4),
#'   verbose = FALSE
#' )
#' analysis <- suppressWarnings(calculate_sait(analysis, method = 'lmm'))
#' analysis <- calculate_effect_sizes(analysis)
#' p_dist <- plot_divergence_distribution(analysis)
#' # print(p_dist)
#'
#' @seealso
#' \code{\link{calculate_effect_sizes}} for computing effect sizes.
#'
#' @export
plot_divergence_distribution <- function(analysis, threshold = 0.1, output_file = NULL,
    width = 12, height = 6, verbose = FALSE, ...) {

    # Load visualization dependencies (ggplot2, cowplot, etc.)
    .load_visualization_deps()

    # Extract verbose parameter if not provided
    verbose <- resolve_slot_param(verbose, analysis@config, "verbose", TRUE)

    # Validate input
    if (!is(analysis, "TSENATAnalysis")) {
        stop("'analysis' must be a TSENATAnalysis object", call. = FALSE)
    }

    # Extract effect sizes from metadata
    if (is.null(analysis@metadata$effect_sizes_divergence)) {
        stop("Effect sizes not found in analysis@metadata$effect_sizes_divergence. ",
            "Run calculate_effect_sizes() first.", call. = FALSE)
    }

    effect_sizes <- analysis@metadata$effect_sizes_divergence

    # Extract interaction results (with effect size columns)
    if (!is.list(effect_sizes) || is.null(effect_sizes$interaction_results)) {
        stop("Invalid effect size structure. Expected @metadata$effect_sizes_divergence$interaction_results",
            call. = FALSE)
    }

    interaction_results <- effect_sizes$interaction_results

    if (!is.data.frame(interaction_results) || nrow(interaction_results) == 0) {
        stop("interaction_results must be a non-empty data frame", call. = FALSE)
    }

    # Create the plot using base function
    p <- tryCatch({
        .plot_divergence_distribution(interaction_results = interaction_results,
            threshold = threshold, ...)
    }, error = function(e) {
        if (verbose) {
            message("plot_divergence_distribution failed: ", e$message)
        }
        return(NULL)
    })

    # If plot creation failed, return NULL invisibly
    if (is.null(p)) {
        return(invisible(NULL))
    }

    # Save to file if requested
    if (!is.null(output_file)) {
        save_analysis_output(p, output_file, object = analysis, verbose = verbose,
            func_name = "plot_divergence_distribution", width = width, height = height)
    }

    # Always return the plot object (not file path) Knitr will auto-manage
    # figure rendering File is saved separately if output_file provided
    invisible(p)
}

