#' Plot Top Transcripts from TSENATAnalysis Object
#'
#' S4 wrapper for \code{plot_top_transcripts()} that extracts data directly from
#' a TSENATAnalysis object. Automatically retrieves the SummarizedExperiment and
#' LM results for visualizing transcript abundance across conditions.
#'
#' @param analysis \code{TSENATAnalysis}. An S4 object containing a processed
#'   SummarizedExperiment and optional LM interaction results.
#'
#' @param gene \code{character} or \code{NULL}. Gene identifier(s) to plot. If a vector 
#'   of multiple genes is provided, plots all of them. If NULL, automatically selects 
#'   the top genes from LM results based on \code{top_n} parameter (genes with lowest p-values).
#'   Default: NULL (auto-extract from lm_results).
#'
#' @param condition_col \code{character}. Column name in colData(se) specifying
#'   group assignments (default: "sample_type").
#'
#' @param top_n \code{numeric}. Number of top transcripts to display for each
#'   condition (default: 3).
#'
#' @param output_file \code{character} or \code{NULL}. Path to save the plot as
#'   a PNG file. If NULL, saves to a temporary location (default: NULL).
#'
#' @param metric \code{character}. Method for ranking transcripts within genes.
#'   One of "median", "mean", "variance", or "iqr" (default: "median").
#'
#' @param verbose \code{logical}. If \code{TRUE}, print diagnostic messages
#'   during plotting (default: FALSE).
#'
#' @return A file path (character) to the saved plot PNG file, invisibly.
#'
#' @details
#' This wrapper extracts the following from \code{analysis}:
#' \describe{
#'   \item{SummarizedExperiment}{From \code{analysis@se} containing transcript counts}
#'   \item{LM results}{From \code{analysis@lm_results$lm_interaction} for gene selection}
#' }
#'
#' If no gene is specified, the function automatically selects the top gene from
#' the LM results (lowest p-value). This simplifies visualization of genes with
#' significant q×condition interaction effects.
#'
#' @examples
#' \dontrun{
#'   # Plot top transcripts for the most significant gene
#'   plot_file <- plot_top_transcripts_s4(
#'     analysis,
#'     top_n = 5
#'   )
#'
#'   # Plot specific gene
#'   plot_file <- plot_top_transcripts_s4(
#'     analysis,
#'     gene = "ENSG00000198888",
#'     top_n = 4,
#'     metric = "mean"
#'   )
#' }
#'
#' @seealso \code{\link{plot_top_transcripts}} for the base function,
#' \code{\link{TSENATAnalysis}} for object structure
#'
#' @export
#' @importFrom methods is
plot_top_transcripts_s4 <- function(
    analysis,
    gene = NULL,
    condition_col = NULL,
    top_n = 3,
    output_file = NULL,
    metric = c("median", "mean", "variance", "iqr"),
    verbose = FALSE) {

  # =========================================================================
  # INPUT VALIDATION
  # =========================================================================
  if (!is(analysis, "TSENATAnalysis")) {
    stop("analysis must be a TSENATAnalysis object", call. = FALSE)
  }

  se <- analysis@se
  if (!inherits(se, "SummarizedExperiment")) {
    stop("analysis@se must be a SummarizedExperiment object", call. = FALSE)
  }

  # =========================================================================
  # AUTO-DETECT condition_col
  # =========================================================================
  if (is.null(condition_col)) {
    cd_cols <- colnames(colData(se))

    # Try Priority 1: @config$condition_col
    if ("condition_col" %in% names(analysis@config)) {
      candidate <- analysis@config$condition_col
      if (candidate %in% cd_cols) {
        condition_col <- candidate
      }
    }

    # Try Priority 2: @config$sample_type or direct colData
    if (is.null(condition_col) && "sample_type" %in% cd_cols) {
      condition_col <- "sample_type"
    }

    # Try Priority 3: @config$condition
    if (is.null(condition_col) && "condition" %in% cd_cols) {
      condition_col <- "condition"
    }

    # Fallback: use first column
    if (is.null(condition_col) && length(cd_cols) > 0) {
      condition_col <- cd_cols[1]
    }

    if (is.null(condition_col)) {
      stop("[plot_top_transcripts_s4] Cannot auto-detect condition_col. ",
           "Provide explicitly.", call. = FALSE)
    }

    if (verbose) {
      cat("[plot_top_transcripts_s4] Auto-detected condition_col =", condition_col, "\n")
    }
  }

  # =========================================================================
  # EXTRACT LM RESULTS (for gene ranking if not specified)
  # =========================================================================
  lm_results_df <- NULL
  if (is.null(gene) && !is.null(analysis@lm_results)) {
    if ("lm_interaction" %in% names(analysis@lm_results)) {
      # Extract results data.frame from list structure
      if (is.data.frame(analysis@lm_results$lm_interaction)) {
        lm_results_df <- analysis@lm_results$lm_interaction
      } else if (is.list(analysis@lm_results$lm_interaction) &&
                 "results" %in% names(analysis@lm_results$lm_interaction)) {
        lm_results_df <- analysis@lm_results$lm_interaction$results
      }

      # Auto-select top gene from LM results
      if (!is.null(lm_results_df) && nrow(lm_results_df) > 0) {
        # Find p-value column (handle various naming conventions)
        p_col <- if ("p_interaction" %in% colnames(lm_results_df)) {
          "p_interaction"
        } else if ("pvalue" %in% colnames(lm_results_df)) {
          "pvalue"
        } else if ("p.value" %in% colnames(lm_results_df)) {
          "p.value"
        } else if ("padj" %in% colnames(lm_results_df)) {
          "padj"
        } else if ("p_value" %in% colnames(lm_results_df)) {
          "p_value"
        } else {
          NA
        }

        # Find gene column (handle various naming conventions)
        gene_col <- if ("gene" %in% colnames(lm_results_df)) {
          "gene"
        } else if ("gene_name" %in% colnames(lm_results_df)) {
          "gene_name"
        } else if ("gene_id" %in% colnames(lm_results_df)) {
          "gene_id"
        } else {
          colnames(lm_results_df)[1]  # Fallback to first column
        }

        if (!is.na(p_col) && p_col %in% colnames(lm_results_df)) {
          # Get top genes (sorted by p-value, select top_n)
          top_indices <- order(lm_results_df[[p_col]])[1:min(top_n, nrow(lm_results_df))]
          gene <- as.character(lm_results_df[top_indices, gene_col])

          if (verbose) {
            cat("[plot_top_transcripts_s4] Auto-selected top", length(gene), "genes from LM results\n")
            cat("  ", paste(gene, collapse = ", "), "\n")
          }
        }
      }
    }
  }

  if (is.null(gene)) {
    stop("[plot_top_transcripts_s4] No gene specified and cannot auto-detect from LM results. ",
         "Provide gene explicitly.", call. = FALSE)
  }

  # =========================================================================
  # CALL BASE FUNCTION
  # =========================================================================
  if (verbose) {
    if (length(gene) > 1) {
      cat("[plot_top_transcripts_s4] Calling plot_top_transcripts() for genes:", 
          paste(gene, collapse = ", "), "\n")
    } else {
      cat("[plot_top_transcripts_s4] Calling plot_top_transcripts() for gene:", gene, "\n")
    }
  }

  plot_file <- tryCatch({
    plot_top_transcripts(
      se = se,
      gene = gene,
      condition_col = condition_col,
      res = lm_results_df,
      top_n = top_n,
      output_file = output_file,
      metric = metric[1]  # Use first metric if multiple provided
    )
  }, error = function(e) {
    stop("[plot_top_transcripts_s4] Error in plotting:\n",
         conditionMessage(e), call. = FALSE)
  })

  if (verbose) {
    cat("[plot_top_transcripts_s4] Plot saved to:", plot_file, "\n")
  }

  # Return file path invisibly
  invisible(plot_file)
}
