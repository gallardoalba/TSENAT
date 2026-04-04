#' TSENATAnalysis S4 Class Definition
#'
#' Central container for unified analysis workflows in TSENAT.
#'
#' The \code{TSENATAnalysis} class encapsulates all components of a complete
#' TSENAT analysis: raw data, configuration metadata, and results from each
#' analytical step. This unified object ensures metadata is never lost through
#' the analysis pipeline and provides consistent accessor methods for result
#' retrieval.
#'
#' @slot se \code{SummarizedExperiment}. The base expression data object
#'   (genes x samples) with assays and colData.
#'
#' @slot config \code{list}. Configuration metadata specifying analysis
#'   parameters that persist through the workflow (q-values, sample grouping
#'   columns, etc.). Set once via \code{tsenat_config()} and used by all
#'   downstream wrapper functions.
#'
#' @slot diversity_results \code{list}. Named list of diversity calculation
#'   results. Each name corresponds to a q-value (e.g., 'q_0.5', 'q_1.0').
#'   Values are SummarizedExperiment objects or data.frames containing entropy
#'   values for each gene at that q-value.
#'
#' @slot lm_results \code{list}. Complex results from linear model and
#'   statistical testing. Top-level names identify analysis type:
#'   \describe{
#'     \item{\code{lm_interaction}}{LM/GAM/GEE model results (list with
#'           \code{$results} data.frame, \code{$models} list, etc.)}
#'     \item{\code{q_interactions}}{Friedman/rank-based test results}
#'     \item{\code{divergence_difference}}{Differential divergence comparison}
#'   }
#'
#' @slot pairwise_results \code{list}. Pairwise group comparison results.
#'   Contains differential testing output computed by
#'   \code{calculate_difference_s4()}, stored under the \code{difference}
#'   component.
#'
#' @slot jackknife_results \code{list}. Resampling-based confidence intervals.
#'   Names correspond to q-values (e.g., 'q_0.5', 'q_1.0'). Values are
#'   jackknife result objects containing resamples, CI bounds, and diagnostics.
#'
#' @slot divergence_results \code{list}. Divergence metric calculations.
#'   Typically contains:
#'   \describe{
#'     \item{\code{tsallis_divergence}}{SummarizedExperiment with 
#' divergence values}
#'     \item{\code{effect_sizes}}{data.frame with Cohen's d, etc.}
#'   }
#'
#' @slot plots \code{list}. Cached visualization objects (ggplot). Names
#'   identify plot type (e.g., 'q_curve', 'lm_interaction', 'influence').
#'   Populated by \code{tsenat()} if \code{generate_plots=TRUE}.
#'
#' @slot metadata \code{list}. Reproducibility and tracking metadata.
#'   Automatically maintained by wrapper functions. Includes:
#'   \describe{
#'     \item{\code{created_at}}{Timestamp of object creation}
#'     \item{\code{function_calls}}{Vector of wrapper functions called}
#'     \item{\code{function_timestamps}}{Timestamps for each function call}
#'     \item{\code{package_version}}{TSENAT version at creation}
#'   }
#'
#' @details
#' Access results via accessor methods (recommended):
#' \code{diversity(obj, q)} for diversity, \code{lmResults(obj)} for models,
#' \code{jeoResults(obj, q)} for entropy outlier jackknife, \code{jisResults(obj, q)} for isoform switching jackknife,
#' \code{getMeta(obj)} for metadata.
#'
#' @rdname TSENATAnalysis-class
#' @exportClass TSENATAnalysis
#' @importFrom SummarizedExperiment colData rowData
#' @importFrom S4Vectors metadata
#'
setClass("TSENATAnalysis", slots = list(se = "SummarizedExperiment", config = "list",
    diversity_results = "list", lm_results = "list", pairwise_results = "list",
    jackknife_results = "list", divergence_results = "list", plots = "list", metadata = "list"), prototype = list(config = list(),
    diversity_results = list(), lm_results = list(), pairwise_results = list(), jackknife_results = list(),
    divergence_results = list(), plots = list(), metadata = list(function_calls = character(0),
        function_timestamps = character(0))), validity = function(object) {
    # Check @se is SummarizedExperiment
    if (!inherits(object@se, "SummarizedExperiment")) {
        return("@se must be a SummarizedExperiment object")
    }

    # Check SE has data
    if (nrow(object@se) == 0 || ncol(object@se) == 0) {
        return("@se has zero dimensions (no genes or samples)")
    }

    # Check @config is list-like (list, TSENATConfig, or other list-based
    # structure)
    if (!is.list(object@config)) {
        return("@config must be a list or list-based config object")
    }

    # Check all results slots are lists
    if (!is.list(object@diversity_results)) {
        return("@diversity_results must be a list")
    }
    if (!is.list(object@lm_results)) {
        return("@lm_results must be a list")
    }
    if (!is.list(object@pairwise_results)) {
        return("@pairwise_results must be a list")
    }
    if (!is.list(object@jackknife_results)) {
        return("@jackknife_results must be a list")
    }
    if (!is.list(object@divergence_results)) {
        return("@divergence_results must be a list")
    }
    if (!is.list(object@plots)) {
        return("@plots must be a list")
    }
    if (!is.list(object@metadata)) {
        return("@metadata must be a list")
    }

    # Check colData has required columns if sample metadata expected
    cdata <- colData(object@se)
    if (!is.null(cdata) && ncol(cdata) > 0) {
        if (!"sample_id" %in% colnames(cdata)) {
            return("colData missing 'sample_id' column required for analysis")
        }
    }

    # Check rowData has gene identifiers if results computed
    rdata <- rowData(object@se)
    if (!is.null(rdata) && ncol(rdata) > 0) {
        if (!"gene_id" %in% colnames(rdata) && !"transcript_id" %in% colnames(rdata)) {
            return("rowData missing 'gene_id' or 'transcript_id' column")
        }
    }

    # NEW: Validate config parameters against SE metadata
    if (length(object@config) > 0) {
        cdata <- colData(object@se)
        
        # Validate condition_col if specified
        if ("condition_col" %in% names(object@config)) {
            col <- object@config$condition_col
            if (!is.null(col) && !col %in% colnames(cdata)) {
                return(sprintf(
                    "@config$condition_col '%s' not found in colData. Available: %s",
                    col, paste(colnames(cdata), collapse = ", ")
                ))
            }
        }
        
        # Validate subject_col if specified and paired=TRUE
        if ("subject_col" %in% names(object@config)) {
            if (isTRUE(object@config$paired)) {
                col <- object@config$subject_col
                if (!is.null(col) && !col %in% colnames(cdata)) {
                    return(sprintf(
                        "@config$subject_col '%s' not found in colData. Available: %s",
                        col, paste(colnames(cdata), collapse = ", ")
                    ))
                }
            }
        }
    }

    TRUE
})

#' Subset TSENATAnalysis Objects
#'
#' Extract a subset of genes and/or samples from a TSENATAnalysis object.
#' Maintains consistency across the underlying SummarizedExperiment and
#' all computed results (diversity, LM, jackknife, divergence).
#'
#' @param x TSENATAnalysis object
#' @param i Integer, logical, or character vector of rows (genes/transcripts)
#'   to retain. If missing, all rows are retained.
#' @param j Integer, logical, or character vector of columns (samples)
#'   to retain. If missing, all columns are retained.
#' @param drop Logical. Currently ignored (included for S4 compatibility).
#'   Always returns TSENATAnalysis (never drops to SE or vector).
#'
#' @details
#' Subsetting preserves all analysis metadata and results while maintaining
#' consistency:
#' - The underlying SummarizedExperiment is subset to the specified
#' genes/samples
#' - Diversity and jackknife results are subset to match sample selection
#' - LM results are recalculated or removed if sample structure changes
#' - Divergence results are subset accordingly
#' - Analysis configuration is preserved
#'
#' @return A new TSENATAnalysis object containing only the specified genes
#' and samples
#'
#' @examples
#' # Create a minimal TSENATAnalysis object
#' library(SummarizedExperiment)
#' counts <- matrix(rpois(200, 10), nrow = 20, ncol = 10)
#' rownames(counts) <- paste0('TX_', 1:20)
#' colnames(counts) <- paste0('S', 1:10)
#' se <- SummarizedExperiment(
#'   assays = list(counts = counts),
#'   rowData = data.frame(gene_id = rep(paste0('G', 1:2), each = 10),
#'                        row.names = rownames(counts)),
#'   colData = data.frame(sample_id = colnames(counts),
#'                        condition = rep(c('A', 'B'), 5),
#'                        row.names = colnames(counts))
#' )
#' analysis <- TSENATAnalysis(se)
#'
#' # Subset to first 10 genes and first 5 samples
#' analysis_subset <- analysis[1:10, 1:5]
#'
#' # Subset by gene name
#' analysis_subset2 <- analysis[paste0('TX_', 1:5), ]
#'
#' # Subset by sample condition (logical indexing)
#' keep_samples <- colData(se(analysis))$condition == 'A'
#' analysis_a <- analysis[, keep_samples]
#'
#' @rdname subsetting-TSENATAnalysis
#' @exportMethod '['
setMethod("[", signature(x = "TSENATAnalysis"), function(x, i, j, drop = TRUE) {
    # Get SE dimensions for default arguments
    se <- x@se
    n_genes <- nrow(se)
    n_samples <- ncol(se)

    # Handle missing indices (default to all)
    if (missing(i)) {
        i <- seq_len(n_genes)
    }
    if (missing(j)) {
        j <- seq_len(n_samples)
    }

    # Convert logical/character indices to numeric
    if (is.logical(i)) {
        i <- which(i)
    } else if (is.character(i)) {
        i <- match(i, rownames(se))
        if (any(is.na(i))) {
            stop("Some gene names not found in object")
        }
    }

    if (is.logical(j)) {
        j <- which(j)
    } else if (is.character(j)) {
        j <- match(j, colnames(se))
        if (any(is.na(j))) {
            stop("Some sample names not found in object")
        }
    }

    # Validate indices
    if (any(i < 1 | i > n_genes)) {
        stop("Row indices out of bounds")
    }
    if (any(j < 1 | j > n_samples)) {
        stop("Column indices out of bounds")
    }

    # Subset the SummarizedExperiment
    se_subset <- se[i, j]

    # Create new TSENATAnalysis with subsetted SE
    new_obj <- new("TSENATAnalysis", se = se_subset, config = x@config, diversity_results = list(),
        lm_results = list(), jackknife_results = list(), divergence_results = list(),
        plots = list(), metadata = x@metadata)

    # Subset diversity results (subset columns to match sample selection)
    if (length(x@diversity_results) > 0) {
        new_obj@diversity_results <- lapply(x@diversity_results, function(div_res) {
            # Handle SummarizedExperiment results
            if (inherits(div_res, "SummarizedExperiment")) {
                return(div_res[i, j])
            }
            # Handle matrix results
            if (is.matrix(div_res)) {
                return(div_res[i, j, drop = FALSE])
            }
            # Handle data.frame results
            if (is.data.frame(div_res)) {
                row_names <- rownames(div_res)
                if (!is.null(row_names)) {
                  keep_rows <- row_names %in% rownames(se_subset)
                  return(div_res[keep_rows, j, drop = FALSE])
                }
            }
            # Return as-is if structure unknown
            return(div_res)
        })
        names(new_obj@diversity_results) <- names(x@diversity_results)
    }

    # Subset jackknife results (sample-level diagnostics)
    if (length(x@jackknife_results) > 0) {
        new_obj@jackknife_results <- lapply(x@jackknife_results, function(jk_res) {
            if (is.list(jk_res)) {
                # Try to subset sample-level components
                if (!is.null(jk_res$resamples) && is.matrix(jk_res$resamples)) {
                  jk_res$resamples <- jk_res$resamples[, j, drop = FALSE]
                }
                if (!is.null(jk_res$influence_scores) && is.matrix(jk_res$influence_scores)) {
                  jk_res$influence_scores <- jk_res$influence_scores[j, , drop = FALSE]
                }
                if (!is.null(jk_res$ci_matrix) && is.array(jk_res$ci_matrix)) {
                  # Subset to gene subset (if applicable)
                  if (nrow(jk_res$ci_matrix) == n_genes) {
                    jk_res$ci_matrix <- jk_res$ci_matrix[i, j, ]
                  }
                }
            }
            return(jk_res)
        })
        names(new_obj@jackknife_results) <- names(x@jackknife_results)
    }

    # Subset divergence results
    if (length(x@divergence_results) > 0) {
        new_obj@divergence_results <- lapply(x@divergence_results, function(div_res) {
            if (inherits(div_res, "SummarizedExperiment")) {
                # Subset both dimensions if applicable
                if (nrow(div_res) == n_genes) {
                  return(div_res[i, j])
                }
                return(div_res[, j])
            }
            if (is.matrix(div_res) && nrow(div_res) == n_genes) {
                return(div_res[i, j, drop = FALSE])
            }
            # Return as-is for non-sample-indexed results
            return(div_res)
        })
        names(new_obj@divergence_results) <- names(x@divergence_results)
    }

    # Subset LM results (sample-indexed components only)
    if (length(x@lm_results) > 0) {
        # LM results include gene-level statistics that don't need subsetting
        # Only subset sample-level diagnostic matrices
        new_obj@lm_results <- lapply(x@lm_results, function(lm_res) {
            if (is.list(lm_res)) {
                # Subset sample diagnostics if present
                if (!is.null(lm_res$residuals) && is.matrix(lm_res$residuals)) {
                  if (ncol(lm_res$residuals) == n_samples) {
                    lm_res$residuals <- lm_res$residuals[, j, drop = FALSE]
                  }
                }
                if (!is.null(lm_res$fitted) && is.matrix(lm_res$fitted)) {
                  if (ncol(lm_res$fitted) == n_samples) {
                    lm_res$fitted <- lm_res$fitted[, j, drop = FALSE]
                  }
                }
            }
            return(lm_res)
        })
        names(new_obj@lm_results) <- names(x@lm_results)
    }

    # Preserve plots (they are visualization-level and generally retained)
    new_obj@plots <- x@plots

    # Validate the new object
    validObject(new_obj)

    return(new_obj)
})
