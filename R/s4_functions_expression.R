#' Plot Top Transcripts from TSENATAnalysis Object
#'
#' S4 wrapper for  \code{. plot_expression()} that 
#' extracts data directly from
#' a TSENATAnalysis object. Automatically retrieves the SummarizedExperiment and
#' SAIT results for visualizing transcript abundance across conditions.
#'
#' @param analysis \code{TSENATAnalysis}. An S4 object containing a processed
#'   SummarizedExperiment and optional SAIT interaction results.
#'
#' @param gene \code{character} or  \code{NULL}.  Gene identifier(s) to plot.
#'  If a vector 
#' of multiple genes is provided, plots all of them. If NULL, automatically
#' selects
#'   the top genes from SAIT results based on \code{top_n} parameter (genes with 
#' lowest p-values).
#'   Default: NULL (auto-extract from sait_results).
#'
#' @param condition_col \code{character}. Column name in colData(se) specifying
#'   group assignments (default: 'sample_type').
#'
#' @param top_n \code{numeric}. Number of top transcripts to display for each
#'   condition (default: 3).
#'
#' @param output_file \code{character} or \code{NULL}. Path to save the plot as
#'   a PNG file. If NULL, saves to a temporary location (default: NULL).
#'
#' @param metric \code{character}. Method for ranking transcripts within genes.
#'   One of 'median', 'mean', 'variance', or 'iqr' (default: 'median').
#'
#' @param width \code{numeric} or \code{NULL}. Output image width in inches. 
#' If NULL, automatically calculated based on number of genes (default: ~13
#' inches per column).
#'
#' @param height \code{numeric} or \code{NULL}. Output image height in inches.
#' If NULL, automatically calculated based on number of genes (default: ~10
#' inches per row + headers).
#'
#' @param fontsize \code{numeric}. Base font size for heatmap titles and labels 
#'   (default: 16pt). Automatically scaled for readability.
#'
#' @param cellwidth \code{numeric}. Width of individual heatmap cells in pixels.
#'   If 0 (default), uses adaptive sizing based on data dimensions and layout.
#'   Set > 0 to override dynamic sizing.
#'
#' @param cellheight \code{numeric}.
#'  Height of individual heatmap cells in pixels.
#'   If 0 (default), uses adaptive sizing based on data dimensions and layout.
#'   Set > 0 to override dynamic sizing.
#'
#' @param layout_ncol \code{numeric}. Number of heatmaps per row in fixed layout
#'   (default: 2). If NULL, uses adaptive layout based on transcript counts.
#'
#' @param use_tpm \code{logical}.  If \code{TRUE},
#'  uses TPM (Transcripts Per Million) 
#' from metadata instead of raw counts (default: FALSE). TPM is normalized
#' for sequencing
#' depth and is recommended for comparing expression across samples.
#' Requires TPM data
#' in metadata from `build_analysis()` or `.build_se()` with `tpm`
#' parameter.
#'   Raises error if TPM not available and `use_tpm = TRUE`.
#'
#' @param quantity \code{character}. Quantity to visualize:
#'   `"abundance"` plots raw counts (or TPM when `use_tpm = TRUE`),
#'   `"usage"` plots each transcript's proportion of its gene total per
#'   sample (default: `"abundance"`).
#'
#' @param verbose \code{logical}. If \code{TRUE}, print diagnostic messages
#'   during plotting (default: FALSE).
#'
#' @param output_file \code{character} or  \code{NULL}.
#'  Optional file path to save the plot.
#'   Supported formats: .pdf, .png, .jpg. Default: NULL (no file output).
#' @param ... Additional arguments passed to the base plotting function.
#'
#' @return Invisibly returns the output file path (if `output_file`
#' provided), or invisible(NULL)
#' if rendering to active graphics device. Graphics are rendered to the
#' active grid device
#'   for capture during vignette compilation.
#'
#' @details
#' This wrapper extracts the following from \code{analysis}:
#' \describe{
#'   \item{SummarizedExperiment}{From \code{analysis@se} containing transcript counts}
#'   \item{SAIT results}{From \code{analysis@sait_results$sait_interaction} for 
#' gene selection}
#' }
#'
#' If no gene is specified, the function automatically selects the top gene from
#' the SAIT results (lowest p-value). This simplifies visualization of genes with
#' significant q x condition interaction effects.
#'
#' @examples
#' # Plot 6: Top transcripts across groups
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
#'   subject_col = 'paired_samples',
#'   paired = TRUE,
#'   control = 'normal',
#'   q_values = seq(0, 2, by = 0.1)
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
#' analysis <- suppressWarnings(calculate_sait(
#'   analysis,
#'   method = 'lmm',
#'   verbose = FALSE
#' ))
#' plot_file <- plot_expression(analysis, top_n = 2)
#' # print(plot_file)
#'
#' @seealso
#' \code{\link{TSENATAnalysis}} for object structure
#'
#' @export
#' @importFrom methods is
plot_expression <- function(analysis, gene = NULL, condition_col = NULL, top_n = 4,
    output_file = NULL, metric = c("median", "mean", "variance", "iqr"), use_tpm = TRUE,
    quantity = c("abundance", "usage"), width = NULL, height = NULL, fontsize = 16,
    cellwidth = 0, cellheight = 0, layout_ncol = 2, verbose = FALSE, ...) {

    # Load visualization dependencies (ggplot2, cowplot, pheatmap, etc.)
    .load_visualization_deps()

    # Validate input FIRST (before accessing analysis@config)
    if (!is(analysis, "TSENATAnalysis")) {
        stop("analysis must be a TSENATAnalysis object", call. = FALSE)
    }

    # Validate parameters per Bioconductor code syntax standards
    metric <- match.arg(metric)
    quantity <- match.arg(quantity)

    # Extract verbose parameter if not provided
    verbose <- resolve_slot_param(verbose, analysis@config, "verbose", FALSE)

    se <- analysis@se
    if (!inherits(se, "SummarizedExperiment")) {
        stop("analysis@se must be a SummarizedExperiment object", call. = FALSE)
    }

    # =========================================================================
    # AUTO-DETECT condition_col
    # =========================================================================
    if (is.null(condition_col)) {
        cd_cols <- colnames(colData(se))
        # Note: Always pass verbose=TRUE for condition_col to ensure users are
        # aware of auto-detection
        condition_col <- auto_detect_column(cd_cols, analysis@config, "condition_col",
            c("condition", "sample_type", "group", "treatment"), verbose = FALSE,
            param_name = "condition_col")
    }

    if (is.null(condition_col)) {
        stop("[plot_expression] Cannot auto-detect condition_col. ",
            "Available colData columns: ", paste(colnames(colData(se)), collapse = ", "),
            ".\n\nSOLUTION: Set @config$condition_col or pass explicit condition_col= parameter.",
            call. = FALSE)
    }

    # =========================================================================
    # EXTRACT LM RESULTS (for gene ranking if not specified)
    # =========================================================================
    sait_results_df <- NULL
    if (is.null(gene) && !is.null(analysis@sait_results)) {
        if ("sait_interaction" %in% names(analysis@sait_results)) {
            # Extract results data.frame from list structure
            if (is.data.frame(analysis@sait_results$sait_interaction)) {
                sait_results_df <- analysis@sait_results$sait_interaction
            } else if (is.list(analysis@sait_results$sait_interaction) && "results" %in%
                names(analysis@sait_results$sait_interaction)) {
                sait_results_df <- analysis@sait_results$sait_interaction$results
            }

            # Auto-select top gene from SAIT results
            if (!is.null(sait_results_df) && nrow(sait_results_df) > 0) {
                # Find p-value and gene columns
                p_col <- auto_detect_column(colnames(sait_results_df), analysis@config,
                  "p_col", c("p_interaction", "padj", "pvalue", "p.value", "p_value"),
                  verbose = FALSE, param_name = "p_col")

                gene_col <- auto_detect_column(colnames(sait_results_df), analysis@config,
                  "gene_col", c("gene", "gene_name", "gene_id"), verbose = FALSE,
                  param_name = "gene_col")

                if (!is.null(p_col) && !is.null(gene_col) && p_col %in% colnames(sait_results_df) &&
                  gene_col %in% colnames(sait_results_df)) {
                  # Get top genes (sorted by p-value, select top_n)
                  top_indices <- order(sait_results_df[[p_col]])[seq_len(min(top_n,
                    nrow(sait_results_df)))]
                  gene <- as.character(sait_results_df[top_indices, gene_col])

                  if (verbose) {
                    message("[plot_expression] Auto-selected top ", length(gene),
                      " genes from SAIT results")
                    message("  ", paste(gene, collapse = ", "))
                  }
                }
            }
        }
    }

    if (is.null(gene)) {
        stop("[plot_expression] No gene specified and cannot auto-detect from SAIT results. ",
            "Provide gene explicitly.", call. = FALSE)
    }

    # =========================================================================
    # CALL BASE FUNCTION
    # =========================================================================
    if (verbose) {
        if (length(gene) > 1) {
            message("[plot_expression] Calling .plot_expression() for genes: ", paste(gene,
                collapse = ", "))
        } else {
            message("[plot_expression] Calling .plot_expression() for gene: ", gene)
        }
    }

    plot_file <- tryCatch({
        .plot_expression(se = se, gene = gene, condition_col = condition_col, res = sait_results_df,
            top_n = top_n, output_file = output_file, metric = metric[1], use_tpm = use_tpm,
            quantity = quantity, width = width, height = height, fontsize = fontsize,
            cellwidth = cellwidth, cellheight = cellheight, layout_ncol = layout_ncol,
            ...)
    }, error = function(e) {
        stop("[plot_expression]", conditionMessage(e), call. = FALSE)
    })

    # Optionally save to file if output_file provided
    if (!is.null(output_file)) {
        if (verbose) {
            message("[plot_expression] Plot saved to: ", output_file)
        }
    }

    # Always return the plot object (ggplot) Knitr will auto-manage figure
    # rendering
    invisible(plot_file)
}

