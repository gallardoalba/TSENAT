#' Plot GAM q-curves from TSENATAnalysis object
#'
#' S4 wrapper that accepts a TSENATAnalysis object and generates GAM q-curve
#' plots
#' @param analysis \code{TSENATAnalysis} object with  diversity and 
#' SAIT interaction results.
#' @param n_top \code{integer}.
#'  Number of top genes (by adjusted p-value) to plot 
#'   (default: 6). Only used if genes = NULL.
#'
#' @param genes \code{character} vector. Optional specific gene names to plot. 
#'   If provided, these genes are plotted directly regardless of significance.
#'   If NULL (default), top n_top significant genes are selected.
#'
#' @param condition_col \code{character}. Column name in colData(se) specifying
#'   group assignments for samples. If NULL, attempts to auto-detect from
#'   \code{@config} (looks for 'condition_col' or 'condition'). If still NULL,
#'   defaults to 'sample_type'.
#'
#' @param sig_alpha \code{numeric}.  Significance threshold for 
#' adjusted p-values 
#'   (default: 0.05). Only used if genes = NULL; filters sait_res to significant 
#'   genes before selecting top n.
#'
#' @param assay_name \code{character}. Name of the assay in se to extract 
#'   (default: 'diversity').
#'
#' @param output_file \code{character} or  \code{NULL}.
#'  Optional file path to save the plot.
#'   Default: NULL (no file output).
#'
#' @param width \code{numeric}. Width of the plot in inches (default: 12).
#'   Only used if output_file is not NULL.
#'
#' @param height \code{numeric} or \code{NULL}. Height of the plot in inches.
#'   Default: NULL (automatically calculated based on width and aspect ratio).
#'   Only used if output_file is not NULL.
#'
#' @param verbose \code{logical}. If TRUE, print diagnostic messages during processing.
#'   (default: FALSE).
#'
#' @param ... Additional arguments passed to the base function.
#'
#' @return A single \code{ggplot} object with all selected genes arranged in a 
#'   grid layout. Can be saved with \code{ggplot2::ggsave()}.
#'
#' @details
#' This wrapper automatically:
#' 1. Extracts SummarizedExperiment from \code{@se} slot
#' 2. Extracts SAIT results from \code{@sait_results$sait_interaction} slot
#' 3. Detects condition_col from \code{@config} or uses default
#' 4. Calls \code{.plot_sait()} with extracted parameters
#'
#' **Parameter Resolution (condition_col):**
#' \enumerate{
#'   \item Explicit \code{condition_col} parameter (highest priority)
#'   \item \code{@config$condition_col} if available
#'   \item \code{@config$condition} if available  
#'   \item Default: 'sample_type'
#' }
#'
#' @seealso
#' \code{\link{calculate_sait}} for 
#' running LM analysis on TSENATAnalysis.
#'
#' @examples
#' \dontrun{
#' # Plot 3: GAM q-curves for genes with q-by-condition interactions
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
#'   q = seq(0.2, 2, by = 0.4)  # 5 unique q-values: 0.2, 0.6, 1.0, 1.4, 1.8
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
#' analysis <- calculate_diversity(analysis, q = seq(0.2, 2, by = 0.4))
#' analysis <- suppressWarnings(calculate_sait(analysis, method = 'gam'))
#' 
#' p_gam <- plot_sait(analysis, n_top = 2, sig_alpha = 0.15)
#' # print(p_gam)
#' }
#'
#' @export
plot_sait <- function(analysis, n_top = 6, genes = NULL, condition_col = NULL, sig_alpha = 0.05,
    assay_name = "diversity", output_file = NULL, width = 12, height = NULL, verbose = FALSE,
    ...) {
    # Load visualization dependencies (ggplot2, cowplot, mgcv, etc.)
    .load_visualization_deps()

    # =========================================================================
    # INPUT VALIDATION
    # =========================================================================
    if (!is(analysis, "TSENATAnalysis")) {
        stop("'analysis' must be a TSENATAnalysis object", call. = FALSE)
    }

    # Check that SAIT results exist
    if (is.null(analysis@sait_results) || is.null(analysis@sait_results$sait_interaction)) {
        stop("[plot_sait] No SAIT interaction results found in @sait_results$sait_interaction. ",
            "Run calculate_sait() first.", call. = FALSE)
    }

    sait_res <- analysis@sait_results$sait_interaction

    if (!is.data.frame(sait_res)) {
        stop("[plot_sait] @sait_results$sait_interaction must be a data.frame", call. = FALSE)
    }

    # Check that diversity results exist (needed for SE reconstruction)
    if (length(analysis@diversity_results) == 0) {
        stop("[plot_sait] No diversity results found in @diversity_results. ", "Run calculate_diversity() first.",
            call. = FALSE)
    }

    # =========================================================================
    # AUTO-DETECT condition_col FROM @config IF NOT PROVIDED
    # =========================================================================
    if (is.null(condition_col)) {
        cd_cols <- colnames(colData(analysis@se))
        condition_col <- auto_detect_column(cd_cols, analysis@config, "condition_col",
            c("condition", "sample_type", "group", "treatment"), verbose = verbose,
            param_name = "condition_col")
    }

    if (is.null(condition_col)) {
        stop("[plot_sait] Cannot auto-detect condition_col. ",
            "Available colData columns: ", paste(colnames(colData(analysis@se)), collapse = ", "),
            ".\n\nSOLUTION: Set @config$condition_col or pass explicit condition_col= parameter.",
            call. = FALSE)
    }

    # Validate that condition_col exists in colData
    if (!(condition_col %in% colnames(colData(analysis@se)))) {
        stop("[plot_sait] Specified condition_col='", condition_col, "' not found in colData. Available columns: ",
            paste(colnames(colData(analysis@se)), collapse = ", "), call. = FALSE)
    }

    # =========================================================================
    # RECONSTRUCT COMBINED DIVERSITY SE FOR PLOTTING (Same approach as in
    # calculate_sait)
    # =========================================================================
    # Extract q-values from diversity_results keys
    q_keys <- names(analysis@diversity_results)
    q_computed <- as.numeric(sub("^q_", "", q_keys))

    # Reconstruct combined diversity SE with all q-values
    diversity_combined <- tryCatch({
        .calculate_diversity(x = analysis@se, q = sort(q_computed), norm = TRUE,
            verbose = verbose, bootstrap = FALSE)
    }, error = function(e) {
        stop("[plot_sait] Failed to reconstruct diversity SE:\n", conditionMessage(e),
            call. = FALSE)
    })

    # =========================================================================
    # EXTRACT model_data FROM STORED RESULTS
    # =========================================================================
    model_data <- NULL
    if ("sait_interaction_model_data" %in% names(analysis@sait_results)) {
        model_data <- analysis@sait_results$sait_interaction_model_data
    }

    # =========================================================================
    # CALCULATE HEIGHT IF NOT PROVIDED
    # =========================================================================
    if (is.null(height)) {
        # Estimate number of genes to be plotted
        if (!is.null(genes)) {
            n_genes_plot <- length(genes)
        } else {
            # Count significant genes
            if ("adj_p_interaction" %in% colnames(sait_res)) {
                sig_genes <- sait_res$adj_p_interaction <= sig_alpha
            } else if ("p_interaction" %in% colnames(sait_res)) {
                sig_genes <- sait_res$p_interaction <= sig_alpha
            } else {
                sig_genes <- rep(TRUE, nrow(sait_res))
            }
            n_genes_plot <- min(sum(sig_genes), n_top)
        }
        # Calculate height: 2 rows per 3-gene group, ~3.5 inches per row
        n_rows <- ceiling(n_genes_plot/2)
        height <- 2 + (3.5 * n_rows)
    }

    # =========================================================================
    # CALL plot_sait_interaction_gam WITH RECONSTRUCTED DIVERSITY SE
    # =========================================================================
    result <- tryCatch({
        .plot_sait(se = diversity_combined, sait_res = sait_res, condition_col = condition_col,
            n_top = n_top, genes = genes, sig_alpha = sig_alpha, assay_name = assay_name,
            model_data = model_data, output_file = output_file, width = width, height = height,
            ...)
    }, error = function(e) {
        stop("[plot_sait]", conditionMessage(e), call. = FALSE)
    })

    # =========================================================================
    # RETURN PLOT
    # =========================================================================
    # Track that plotting occurred
    if (is.list(analysis@metadata)) {
        analysis@metadata$function_calls <- c(analysis@metadata$function_calls, paste0("plot_sait[n_top=",
            n_top, ", condition_col=", condition_col, "]"))
    }

    # Save plot to file if requested (only if result is a valid ggplot)
    if (!is.null(output_file) && inherits(result, "ggplot")) {
        save_analysis_output(result, output_file, object = analysis, verbose = verbose,
            func_name = "plot_sait", width = width, height = height)
    }

    # Return the plot object directly (not the analysis object)
    result
}


