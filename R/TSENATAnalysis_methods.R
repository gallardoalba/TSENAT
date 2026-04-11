setMethod("getSE", "TSENATAnalysis", function(object) {
    object@se
})

setMethod("getMeta", "TSENATAnalysis", function(object, key = NULL) {
    # Returns ONLY essential metadata (timestamps, version, workflow type)
    # Large result tables, sample stats, function logs are accessed via results(), etc.
    
    # Filter metadata to essential fields only
    essential_meta <- list()
    
    if ("created_at" %in% names(object@metadata)) {
        essential_meta$created_at <- object@metadata$created_at
    }
    if ("ended_at" %in% names(object@metadata)) {
        essential_meta$ended_at <- object@metadata$ended_at
    }
    if ("package_version" %in% names(object@metadata)) {
        essential_meta$package_version <- object@metadata$package_version
    }
    if ("tsenat_version" %in% names(object@metadata)) {
        essential_meta$tsenat_version <- object@metadata$tsenat_version
    }
    if ("workflow_type" %in% names(object@metadata)) {
        essential_meta$workflow_type <- object@metadata$workflow_type
    }
    if ("workflow" %in% names(object@metadata)) {
        # Extract only essential workflow info
        workflow_info <- object@metadata$workflow
        if (is.list(workflow_info)) {
            essential_meta$workflow <- list(
                type = workflow_info$workflow_type,
                completion_time = workflow_info$completion_time
            )
        }
    }
    
    # If key specified, return specific field
    if (!is.null(key)) {
        return(essential_meta[[key]])
    }
    
    # Return filtered essential metadata only
    essential_meta
})

setMethod("getConfig", "TSENATAnalysis", function(object, key = NULL) {
    if (is.null(key)) {
        return(object@config)
    }
    object@config[[key]]
})

setMethod("getPlot", "TSENATAnalysis", function(object, type = NULL) {
    if (is.null(type)) {
        return(object@plots)
    }
    object@plots[[type]]
})

setMethod("addPlot", "TSENATAnalysis", function(object, type, plot, replace = FALSE) {
    if (!replace && type %in% names(object@plots)) {
        warning("Plot type '", type, "' already exists. Set replace=TRUE to overwrite.",
            call. = FALSE)
        return(object)
    }
    object@plots[[type]] <- plot
    object
})

setMethod("show", "TSENATAnalysis", function(object) {
    n_genes <- nrow(object@se)
    n_samples <- ncol(object@se)
    has_diversity <- length(object@diversity_results) > 0
    has_lm <- length(object@lm_results) > 0
    has_jackknife <- length(object@jackknife_results) > 0
    has_divergence <- length(object@divergence_results) > 0
    status_parts <- c()
    if (has_diversity) status_parts <- c(status_parts, "DIVERSITY")
    if (has_lm) status_parts <- c(status_parts, "LM")
    if (has_jackknife) status_parts <- c(status_parts, "JACKKNIFE")
    if (has_divergence) status_parts <- c(status_parts, "DIVERGENCE")
    status_str <- if (length(status_parts) > 0) {
        paste(status_parts, collapse = " | ")
    } else {
        "EMPTY"
    }
    message("TSENATAnalysis Object")
    message("====================")
    message(sprintf("Genes:       %d", n_genes))
    message(sprintf("Samples:     %d", n_samples))
    message(sprintf("Configuration: %d parameters", length(object@config)))
    message(sprintf("Analysis status: %s", status_str))
    if (length(object@metadata) > 0 && !is.null(object@metadata$created_at)) {
        message(sprintf("Created: %s", object@metadata$created_at))
    }
    message("")
})

setMethod("summary", "TSENATAnalysis", function(object) {
    message("TSENAT Analysis Summary")
    message("=======================\n")
    message("DATA STRUCTURE:")
    message(sprintf("  Genes:        %d", nrow(object@se)))
    message(sprintf("  Samples:      %d", ncol(object@se)))
    assay_names <- paste(SummarizedExperiment::assayNames(object@se), collapse = ", ")
    message(sprintf("  Assays:       %s", assay_names))
    message("\nCONFIGURATION:")
    if (length(object@config) > 0) {
        message(sprintf("  Parameters: %d", length(object@config)))
    }
    message("\nANALYSIS RESULTS:")
    if (length(object@diversity_results) > 0) {
        message("  [OK] Diversity: computed")
    } else {
        message("  [--] Diversity: not computed")
    }
    if (length(object@lm_results) > 0) {
        message("  [OK] LM Interaction: computed")
    } else {
        message("  [--] LM Interaction: not computed")
    }
    if (length(object@jackknife_results) > 0) {
        message("  [OK] Jackknife: computed")
    } else {
        message("  [--] Jackknife: not computed")
    }
    if (length(object@divergence_results) > 0) {
        message("  [OK] Divergence: computed")
    } else {
        message("  [--] Divergence: not computed")
    }
    invisible(object)
})

setMethod("setConfig", "TSENATAnalysis", function(object, value) {
    if (inherits(value, "TSENATConfig")) {
        value <- unclass(value)
    }
    if (!is.list(value)) {
        stop("Configuration must be a list or TSENATConfig object", call. = FALSE)
    }
    object@config <- value
    validObject(object)
    object
})

setMethod("setConfigValue", "TSENATAnalysis", function(object, key, value) {
    config <- getConfig(object)
    if (is.null(config)) {
        config <- list()
    }
    config[[key]] <- value
    setConfig(object, config)
})

#' Extract SummarizedExperiment from TSENATAnalysis
#'
#' @title Extract SummarizedExperiment Slot
#'
#' @param object \code{TSENATAnalysis} object.
#'
#' @return The \code{SummarizedExperiment} object stored in \code{@se} slot,
#'   containing raw count data and sample metadata.
#'
#' @details
#' Provides read-only access to the underlying SummarizedExperiment. To modify
#' the SummarizedExperiment, use the configuration methods or directly access
#' \code{object@se}.
#'
#' @aliases se,TSENATAnalysis-method
#' @rdname TSENATAnalysis-se
setMethod("se", "TSENATAnalysis", function(object) {
    object@se
})

#' Get or Set Metadata
#'
#' @title Metadata Accessor Methods
#'
#' @description Access or set the metadata list stored in a TSENATAnalysis object.
#' The getter function retrieves all metadata or a specific key-value.
#' The setter function replaces the entire metadata list.
#'
#' @param x A \code{\linkS4class{TSENATAnalysis}} object
#' @param key Optional character string specifying a metadata key to retrieve
#' @param value A list of metadata to assign
#'
#' @return
#' \describe{
#'   \item{\code{metadata}}{Returns the full metadata list, or a single value if \code{key} is specified}
#'   \item{\code{metadata<-}}{Returns the modified \code{TSENATAnalysis} object}
#' }
#'
#' @examples
#' # Create a TSENATAnalysis object
#' library(SummarizedExperiment)
#' se <- SummarizedExperiment(assays = list(counts = matrix(1:100, nrow = 10)))
#' analysis <- new('TSENATAnalysis', se = se, config = list())
#'
#' # Get metadata (empty by default)
#' metadata(analysis)
#'
#' \donttest{
#' # Set metadata
#' metadata(analysis) <- list(processing_date = Sys.Date(), method = "test")
#'
#' # Retrieve all metadata
#' metadata(analysis)
#'
#' # Retrieve specific metadata key
#' metadata(analysis, key = "method")
#' }
#'
#' @rdname TSENATAnalysis-metadata
#' @exportMethod metadata
setMethod("metadata", "TSENATAnalysis", function(x, key = NULL) {
    if (is.null(key)) {
        return(x@metadata)
    }
    x@metadata[[key]]
})

#' @aliases metadata<-,TSENATAnalysis-method
#' @rdname TSENATAnalysis-metadata
#' @exportMethod 'metadata<-'
setReplaceMethod("metadata", "TSENATAnalysis", function(x, value) {
    x@metadata <- value
    x
})

#' @param i Gene indices (numeric, logical, or character vector). Defaults to all genes.
#' @param j Sample indices (numeric, logical, or character vector). Defaults to all samples.
#' @param drop Ignored for TSENATAnalysis objects; included for S4 method signature compatibility.
#' @param x A \code{\linkS4class{TSENATAnalysis}} object to subset.
#'
#' @rdname TSENATAnalysis-class
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
