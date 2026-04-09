#' TSENATAnalysis S4 Class
#'
#' Central container for unified analysis workflows in TSENAT.
#'
#' The \code{TSENATAnalysis} class encapsulates all components of a complete
#' TSENAT analysis: raw data, configuration metadata, and results from each
#' analytical step. This unified object ensures metadata is never lost through
#' the analysis pipeline and provides consistent accessor methods for result
#' retrieval.
#'
#' @return An S4 object of class \code{TSENATAnalysis} containing:
#'   \itemize{
#'     \item Raw data (SummarizedExperiment)
#'     \item Analysis configuration
#'     \item Results from diversity,  linear model,  jackknife,  and 
#' divergence analyses
#'     \item Cached visualization objects
#'     \item Reproducibility metadata and function history
#'   }
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
#'     \item{\code{rank_test}}{Friedman/rank-based test results}
#'     \item{\code{divergence_difference}}{Differential divergence comparison}
#'   }
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
#' @section Accessor Methods:
#'   \describe{
#'     \item{\code{results(object, type, ...)}}{Unified interface to extract all analysis results (diversity, divergence, lm, jackknife, rank_test, effect_sizes_divergence, switching_tables)}
#'     \item{\code{getPlot(object, type=NULL)}}{Retrieve cached plot}
#'     \item{\code{addPlot(object, type, plot)}}{Add/cache a new plot}
#'     \item{\code{show(object)}}{Display object summary}
#'     \item{\code{summary(object)}}{Get detailed analysis summary}
#'   }
#'
#' @section Validation:
#'   Validity is checked at object construction. Ensures @se is a
#'   SummarizedExperiment and all slots are correct types.
#'
#' @examples
#' # Load real TSENAT data
#' data(readcounts)
#' metadata_df <- read.table(system.file('extdata', 'metadata.tsv', package
#' = 'TSENAT'),
#'   header = TRUE, sep = '\t')
#' gff3_file <- system.file('extdata', 'annotation.gff3.gz', package = 'TSENAT')
#' config <- tsenat_config(sample_col = 'sample', condition_col = 'condition')
#' analysis <- build_analysis_s4(readcounts = readcounts, tx2gene =
#' gff3_file, metadata = metadata_df, config = config,
#'   tpm = tpm, effective_length = effective_length)
#' analysis <- filter_analysis_s4(analysis, min_samples = 1, subset_n_genes
#' = 200)
#' @name TSENATAnalysis-class
#' @rdname TSENATAnalysis-class
#' @exportClass TSENATAnalysis
#' 
NULL

#' Constructor for TSENATAnalysis objects
#'
#' Creates a new TSENATAnalysis object with a SummarizedExperiment base
#' and optional initial configuration.
#'
#' @param se \code{SummarizedExperiment}. The base expression data object.
#' @param config \code{list}. Optional initial configuration (usually set
#'   via \code{tsenat_config()} instead).
#'
#' @return A new \code{TSENATAnalysis} object.
#'
#' @details
#' The constructor initializes all slots with empty lists except @se,
#' which must be provided. The @metadata slot automatically records:
#' - creation timestamp
#' - TSENAT package version
#' - initial function call
#'
#' @importFrom SummarizedExperiment SummarizedExperiment colData rowData
#' @importFrom S4Vectors metadata
#'
#' @examples
#' # Load real TSENAT data
#' data(readcounts)
#' metadata_df <- read.table(system.file('extdata', 'metadata.tsv', package
#' = 'TSENAT'),
#'   header = TRUE, sep = '\t')
#' gff3_file <- system.file('extdata', 'annotation.gff3.gz', package = 'TSENAT')
#' config <- tsenat_config(sample_col = 'sample', condition_col = 'condition')
#' analysis <- build_analysis_s4(readcounts = readcounts, tx2gene =
#' gff3_file, metadata = metadata_df, config = config,
#'   tpm = tpm, effective_length = effective_length)
#' analysis <- filter_analysis_s4(analysis, min_samples = 1, subset_n_genes
#' = 200)
#' @export
TSENATAnalysis <- function(se, config = list()) {
    # Validate input
    if (!inherits(se, "SummarizedExperiment")) {
        stop("se must be a SummarizedExperiment object", call. = FALSE)
    }

    # Convert TSENATConfig S4 object to list if needed
    if (inherits(config, "TSENATConfig")) {
        config <- unclass(config)
    }

    # NEW: Validate config parameters early (before object creation)
    if (length(config) > 0 && !is.null(config)) {
        cdata <- SummarizedExperiment::colData(se)

        # Validate condition_col
        if ("condition_col" %in% names(config)) {
            col <- config$condition_col
            if (!is.null(col) && !col %in% colnames(cdata)) {
                stop(sprintf("Invalid condition_col '%s': column not found in colData.\\nAvailable columns: %s",
                  col, paste(colnames(cdata), collapse = ", ")), call. = FALSE)
            }
        }

        # Validate subject_col if paired
        if ("paired" %in% names(config) && isTRUE(config$paired)) {
            if ("subject_col" %in% names(config)) {
                col <- config$subject_col
                if (!is.null(col) && !col %in% colnames(cdata)) {
                  stop(sprintf("Invalid subject_col '%s': column not found in colData.\\nAvailable columns: %s",
                    col, paste(colnames(cdata), collapse = ", ")), call. = FALSE)
                }
            }
        }
    }

    # Ensure sample_id column exists in colData (required by validator)
    if (!"sample_id" %in% colnames(SummarizedExperiment::colData(se))) {
        SummarizedExperiment::colData(se)$sample_id <- colnames(se)
    }

    # Create new object with all slots initialized
    new("TSENATAnalysis", se = se, config = if (length(config) > 0)
        config else list(), diversity_results = list(), lm_results = list(), jackknife_results = list(),
        divergence_results = list(), plots = list(), metadata = list(created_at = Sys.time(),
            package_version = as.character(utils::packageVersion("TSENAT")), function_calls = character()))
}

# ============================================================================
# PLOT ACCESSORS
# ============================================================================

#' Get cached plot
#'
#' @param object \code{TSENATAnalysis} object.
#' @param type \code{character}. Plot type: 'q_curve', 'lm_interaction',
#'   'divergence', 'influence', 'volcano', etc.
#'   If NULL, returns all cached plots.
#'
#' @return ggplot object or list of plots.
#'
#' @examples
#' # Load real TSENAT data
#' data(readcounts)
#' metadata_df <- read.table(system.file('extdata', 'metadata.tsv', package
#' = 'TSENAT'),
#'   header = TRUE, sep = '\t')
#' gff3_file <- system.file('extdata', 'annotation.gff3.gz', package = 'TSENAT')
#' config <- tsenat_config(sample_col = 'sample', condition_col = 'condition')
#' analysis <- build_analysis_s4(readcounts = readcounts, tx2gene =
#' gff3_file, metadata = metadata_df, config = config,
#'   tpm = tpm, effective_length = effective_length)
#' analysis <- filter_analysis_s4(analysis, min_samples = 1, subset_n_genes
#' = 200)
#' all_plots <- getPlot(analysis)
#'
#' @rdname TSENATAnalysis-methods
#' @export
setMethod("getPlot", "TSENATAnalysis", function(object, type = NULL) {
    if (is.null(type)) {
        # Return all plots (empty list if none exist)
        return(object@plots)
    }

    if (length(object@plots) == 0) {
        warning("No plots found. Run tsenat() with generate_plots=TRUE.")
        return(NULL)
    }

    if (!(type %in% names(object@plots))) {
        warning("Plot type '", type, "' not found.\n", "Available: ", paste(names(object@plots),
            collapse = ", "))
        return(NULL)
    }

    object@plots[[type]]
})

#' Add plot to cache
#'
#' @param object \code{TSENATAnalysis} object.
#' @param type \code{character}. Name/type for this plot.
#' @param plot Object to cache (ggplot, etc.).
#' @param replace \code{logical}. If TRUE, replace existing plot of same type.
#'
#' @return Modified TSENATAnalysis object with plot cached.
#'
#' @examples
#' # Demonstrates adding a plot to analysis object (requires ggplot2)
#' library(SummarizedExperiment)
#' se <- SummarizedExperiment(
#'   assays = list(counts = matrix(rpois(100, 10), nrow = 10, ncol = 10)),
#'   colData = data.frame(sample_id = paste0('S', 1:10))
#' )
#' analysis <- TSENATAnalysis(se)
#' # analysis <- addPlot(analysis, type = 'example', plot = NULL)
#'
#' @rdname TSENATAnalysis-methods
#' @export
setMethod("addPlot", "TSENATAnalysis", function(object, type, plot, replace = FALSE) {
    if (!replace && type %in% names(object@plots)) {
        warning("Plot type '", type, "' already exists. Set replace=TRUE to overwrite.",
            call. = FALSE)
        return(object)
    }

    object@plots[[type]] <- plot
    object
})

# ============================================================================
# SHOW METHOD
# ============================================================================

#' Display TSENATAnalysis object
#'
#' @param object \code{TSENATAnalysis} object to display.
#'
#' @return Invisibly returns the \code{TSENATAnalysis} object (called for its
#'   side effect of printing a formatted summary to the console).
#'
#' @details
#' Provides a concise summary of: data dimensions, completed analyses,
#' number of results, and metadata tracking.
#'
#' @examples
#' # Load real TSENAT data
#' data(readcounts)
#' metadata_df <- read.table(system.file('extdata', 'metadata.tsv', package
#' = 'TSENAT'),
#'   header = TRUE, sep = '\t')
#' gff3_file <- system.file('extdata', 'annotation.gff3.gz', package = 'TSENAT')
#' config <- tsenat_config(sample_col = 'sample', condition_col = 'condition')
#' analysis <- build_analysis_s4(readcounts = readcounts, tx2gene =
#' gff3_file, metadata = metadata_df, config = config,
#'   tpm = tpm, effective_length = effective_length)
#' analysis <- filter_analysis_s4(analysis, min_samples = 1, subset_n_genes
#' = 200)
#' show(analysis)
#'
#' @rdname TSENATAnalysis-methods
#' @export
setMethod("show", "TSENATAnalysis", function(object) {
    n_genes <- nrow(object@se)
    n_samples <- ncol(object@se)
    
    # Quick status check - all calculation functions populate their respective slots
    has_diversity <- length(object@diversity_results) > 0
    has_lm <- length(object@lm_results) > 0
    has_jackknife <- length(object@jackknife_results) > 0
    has_divergence <- length(object@divergence_results) > 0
    
    # Build status string
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

# ============================================================================
# SUMMARY METHOD
# ============================================================================

#' Detailed summary of TSENATAnalysis
#'
#' @param object \code{TSENATAnalysis} object.
#'
#' @return Prints detailed summary (invisibly returns object).
#'
#' @details
#' Provides comprehensive analysis summary including dimensions, results
#' counts, validation status, and metadata tracking.
#'
#' @examples
#' # Load real TSENAT data
#' data(readcounts)
#' metadata_df <- read.table(system.file('extdata', 'metadata.tsv', package
#' = 'TSENAT'),
#'   header = TRUE, sep = '\t')
#' gff3_file <- system.file('extdata', 'annotation.gff3.gz', package = 'TSENAT')
#' config <- tsenat_config(sample_col = 'sample', condition_col = 'condition')
#' analysis <- build_analysis_s4(readcounts = readcounts, tx2gene =
#' gff3_file, metadata = metadata_df, config = config,
#'   tpm = tpm, effective_length = effective_length)
#' analysis <- filter_analysis_s4(analysis, min_samples = 1, subset_n_genes
#' = 200)
#' summary(analysis)
#'
#' @rdname TSENATAnalysis-methods
#' @export
setMethod("summary", "TSENATAnalysis", function(object) {
    message("TSENAT Analysis Summary")
    message("=======================\n")

    # === DATA SECTION ===
    message("DATA STRUCTURE:")
    message(sprintf("  Genes:        %d", nrow(object@se)))
    message(sprintf("  Samples:      %d", ncol(object@se)))
    assay_names <- paste(SummarizedExperiment::assayNames(object@se), collapse = ", ")
    message(sprintf("  Assays:       %s", assay_names))

    # === CONFIGURATION SECTION (Grouped) ===
    message("\nCONFIGURATION:")
    if (length(object@config) > 0) {
        cfg <- object@config
        
        # Design and filtering
        message("  Design & Filtering:")
        if ("paired" %in% names(cfg)) {
            message(sprintf("    - Design: %s", if (isTRUE(cfg$paired)) "paired" else "unpaired"))
        }
        if ("stringency" %in% names(cfg)) {
            message(sprintf("    - Stringency: %s", cfg$stringency %||% "default"))
        }
        if ("norm" %in% names(cfg)) {
            message(sprintf("    - Normalization: %s", if (isTRUE(cfg$norm)) "enabled" else "disabled"))
        }
        
        # Diversity computation
        message("  Diversity Metrics:")
        if ("q_values" %in% names(cfg)) {
            q_vals <- cfg$q_values
            message(sprintf("    - Q-spectrum: %.2f to %.2f (%d values)", min(q_vals), max(q_vals), length(q_vals)))
        }
        if ("pseudocount" %in% names(cfg)) {
            pc_val <- if (cfg$pseudocount == 0) "disabled" else as.character(cfg$pseudocount)
            message(sprintf("    - Pseudocount: %s", pc_val))
        }
        
        # Statistical methods
        message("  Statistical Methods:")
        if ("lm_method" %in% names(cfg)) {
            message(sprintf("    - LM fitting: %s", toupper(cfg$lm_method %||% "GAM")))
        }
        if ("lm_pcorr" %in% names(cfg)) {
            message(sprintf("    - P-value correction: %s", toupper(cfg$lm_pcorr %||% "BH")))
        }
        if ("jis_use_lm_fdr" %in% names(cfg)) {
            message(sprintf("    - Jackknife filtering: %s", if (isTRUE(cfg$jis_use_lm_fdr)) "LM-based" else "all genes"))
        }
        
        # Bootstrap configuration
        if (isTRUE(cfg$bootstrap)) {
            message("  Bootstrap & Confidence Intervals:")
            message(sprintf("    - Method: %s", toupper(cfg$bootstrap_method %||% "PERCENTILE")))
            message(sprintf("    - Replicates: %d", cfg$nboot %||% 1000))
            message(sprintf("    - CI level: %.2f", cfg$bootstrap_ci %||% 0.95))
        }
    }

    # === ANALYSIS RESULTS SECTION ===
    message("\nANALYSIS RESULTS:")
    
    # Check diversity results
    if (length(object@diversity_results) > 0) {
        if (is.matrix(object@diversity_results)) {
            n_genes <- ncol(object@diversity_results)
            n_q <- nrow(object@diversity_results)
            message(sprintf("  ✓ Diversity: %d genes × %d q-values", n_genes, n_q))
        } else if (is.list(object@diversity_results) && length(object@diversity_results) > 0) {
            message(sprintf("  ✓ Diversity: %d q-value(s)", length(object@diversity_results)))
        }
    } else {
        message("  ✗ Diversity: not computed")
    }
    
    # Check LM results
    if (length(object@lm_results) > 0) {
        n_genes <- 0
        for (component in names(object@lm_results)) {
            if (is.list(object@lm_results[[component]]) && 
                !is.null(object@lm_results[[component]]$pvalue_results)) {
                pvals <- object@lm_results[[component]]$pvalue_results
                if (is.data.frame(pvals) && nrow(pvals) > 0) {
                    n_genes <- nrow(pvals)
                    break
                }
            }
        }
        if (n_genes > 0) {
            message(sprintf("  ✓ LM Interaction: %d genes", n_genes))
        } else {
            message("  ✓ LM Interaction: results stored")
        }
    } else {
        message("  ✗ LM Interaction: not computed")
    }
    
    # Check jackknife results
    if (length(object@jackknife_results) > 0) {
        n_jackknife <- 0
        if (is.list(object@jackknife_results)) {
            if (!is.null(object@jackknife_results$switching_summary) &&
                is.data.frame(object@jackknife_results$switching_summary)) {
                n_jackknife <- nrow(object@jackknife_results$switching_summary)
            }
        }
        if (n_jackknife > 0) {
            message(sprintf("  ✓ Jackknife Switching: %d genes", n_jackknife))
        } else {
            message(sprintf("  ✓ Jackknife Switching: %d q-value(s)", length(object@jackknife_results)))
        }
    } else {
        message("  ✗ Jackknife Switching: not computed")
    }
    
    # Check divergence results
    if (length(object@divergence_results) > 0) {
        if (is.data.frame(object@divergence_results)) {
            n_div <- nrow(object@divergence_results)
            message(sprintf("  ✓ Divergence Metrics: %d genes", n_div))
        } else if (is.list(object@divergence_results)) {
            message(sprintf("  ✓ Divergence Metrics: %s", paste(names(object@divergence_results), collapse = ", ")))
        } else {
            message("  ✓ Divergence Metrics: computed")
        }
    } else {
        message("  ✗ Divergence Metrics: not computed")
    }
    
    # Check visualizations
    if (length(object@plots) > 0) {
        message(sprintf("  ✓ Visualizations: %d plot(s)", length(object@plots)))
    } else {
        message("  ✗ Visualizations: not generated")
    }

    # === METADATA & PROCESSING SECTION ===
    message("\nPROCESSING & METADATA:")
    if (!is.null(object@metadata$created_at)) {
        message(sprintf("  Created: %s", format(object@metadata$created_at, "%Y-%m-%d %H:%M:%S")))
    }
    if (!is.null(object@metadata$ended_at)) {
        message(sprintf("  Completed: %s", format(object@metadata$ended_at, "%Y-%m-%d %H:%M:%S")))
    }
    if (!is.null(object@metadata$package_version)) {
        message(sprintf("  Package version: %s", object@metadata$package_version))
    }
    if (length(object@metadata$function_calls) > 0) {
        # Compress workflow: group repeated function calls by function name
        calls <- object@metadata$function_calls
        
        # Extract function names and build compact summary
        compressed_calls <- character()
        i <- 1
        while (i <= length(calls)) {
            # Extract function name (before '[')
            current_call <- calls[i]
            func_name <- sub("\\[.*", "", current_call)
            
            # Count consecutive calls of same function
            j <- i
            while (j <= length(calls) && sub("\\[.*", "", calls[j]) == func_name) {
                j <- j + 1
            }
            
            # Group summary: if >3 consecutive calls, show compacted form
            n_consecutive <- j - i
            if (n_consecutive > 3 && grepl("=", current_call)) {
                # Extract first parameter from first call
                first_param <- sub(".*\\[([^\\]]+).*", "\\1", calls[i])
                last_call <- calls[j - 1]
                last_param <- sub(".*\\[([^\\]]+).*", "\\1", last_call)
                
                # Check if parameters are q-values
                if (grepl("q=", first_param) && grepl("q=", last_param)) {
                    first_q <- as.numeric(sub(".*q=([0-9.]+).*", "\\1", first_param))
                    last_q <- as.numeric(sub(".*q=([0-9.]+).*", "\\1", last_param))
                    compressed_calls <- c(compressed_calls, 
                        sprintf("%s (%d q-values: %.1f-%.1f)", func_name, n_consecutive, first_q, last_q))
                } else {
                    compressed_calls <- c(compressed_calls, 
                        sprintf("%s (×%d)", func_name, n_consecutive))
                }
            } else {
                # Keep individual calls if 3 or fewer
                for (k in seq(i, j - 1)) {
                    compressed_calls <- c(compressed_calls, calls[k])
                }
            }
            
            i <- j
        }
        
        # Format workflow vertically for better readability
        message("  Workflow:")
        for (step in compressed_calls) {
            message(sprintf("    → %s", step))
        }
    }

    message("")

    invisible(object)
})

# ============================================================================
# CONFIGURATION ACCESSORS (GAP 3 FIX)
# ============================================================================

#' Get configuration from TSENATAnalysis
#'
#' Retrieve the configuration parameters stored in a TSENATAnalysis object.
#' These parameters control analysis behavior including q-values, normalization,
#' and output settings.
#'
#' @param object \code{TSENATAnalysis} object.
#'
#' @return \code{list} containing configuration parameters (q_values, method, etc.)
#'
#' @details
#' Configuration is stored in the @config slot and controls how downstream
#' analyses are performed. Use \code{\link{setConfig}} to replace the entire
#' configuration or \code{\link{setConfigValue}} for targeted updates.
#'
#' @seealso
#' \code{\link{setConfig}} for replacing configuration,
#' \code{\link{setConfigValue}} for single value updates,
#' \code{\link{tsenat_config}} for creating configuration objects
#'
#' @examples
#' data(readcounts)
#' metadata_df <- read.table(system.file('extdata', 'metadata.tsv', 
#'   package = 'TSENAT'), header = TRUE, sep = '\t')
#' gff3_file <- system.file('extdata', 'annotation.gff3.gz', package = 'TSENAT')
#' config <- tsenat_config(sample_col = 'sample', condition_col = 'condition')
#' analysis <- build_analysis_s4(readcounts = readcounts, tx2gene = gff3_file,
#'   metadata = metadata_df, config = config, tpm = tpm, effective_length = effective_length)
#' config <- getConfig(analysis)
#' print(config$q_values)
#'
#' @rdname TSENATAnalysis-methods
#' @export
setGeneric("getConfig", function(object) {
    standardGeneric("getConfig")
})

#' @rdname TSENATAnalysis-methods
#' @export
setMethod("getConfig", "TSENATAnalysis", function(object) {
    object@config
})

#' Replace analysis configuration
#'
#' @param object \code{TSENATAnalysis} object.
#' @param value \code{list}. New configuration. Should contain q_values, 
#'   column name specifications, and/or other parameters.
#'
#' @return Modified \code{TSENATAnalysis} object with updated @config.
#'
#' @details
#' Provides type-safe replacement of @config slot. Typically called once
#' at the start of an analysis via \code{tsenat_config()} rather than directly.
#'
#' Replace entire configuration in TSENATAnalysis
#'
#' @param object \code{TSENATAnalysis} object.
#' @param value List or 
#' \code{TSENATConfig} object containing configuration settings.
#'
#' @return Updated \code{TSENATAnalysis} object with replaced configuration.
#'
#' @details
#' Replaces the entire configuration of a TSENATAnalysis object. This method
#' is useful when you need to apply a new set of configuration parameters to
#' an existing analysis object. All previous configuration values are replaced
#' with those in the new value object.
#'
#' @examples
#' # Load real TSENAT data
#' data(readcounts)
#' metadata_df <- read.table(system.file('extdata', 'metadata.tsv', package
#' = 'TSENAT'),
#'   header = TRUE, sep = '\t')
#' gff3_file <- system.file('extdata', 'annotation.gff3.gz', package = 'TSENAT')
#' 
#' # Build and subset analysis
#' config <- tsenat_config(sample_col = 'sample', condition_col = 'condition')
#' analysis <- build_analysis_s4(readcounts = readcounts, tx2gene =
#' gff3_file, metadata = metadata_df, config = config,
#'   tpm = tpm, effective_length = effective_length)
#' analysis <- filter_analysis_s4(analysis, min_samples = 1, subset_n_genes
#' = 200)
#' 
#' # Replace configuration with new settings
#' new_config <- tsenat_config(q_values = c(0.5, 1.0, 1.5), seed = 42)
#' analysis <- setConfig(analysis, new_config)
#' 
#' # Verify the new configuration was applied
#' current_config <- getConfig(analysis)
#' print(current_config$q_values)  # Shows c(0.5, 1.0, 1.5)
#'
#' @rdname TSENATAnalysis-methods
#' @export
setGeneric("setConfig", function(object, value) {
    standardGeneric("setConfig")
})

#' @rdname TSENATAnalysis-methods
#' @aliases setConfig,TSENATAnalysis-method
#' @keywords internal
#'
#' @examples
#' # Load real TSENAT data
#' data(readcounts)
#' metadata_df <- read.table(system.file('extdata', 'metadata.tsv', package
#' = 'TSENAT'),
#'   header = TRUE, sep = '\t')
#' gff3_file <- system.file('extdata', 'annotation.gff3.gz', package = 'TSENAT')
#' 
#' # Build and subset analysis
#' config <- tsenat_config(sample_col = 'sample', condition_col = 'condition')
#' analysis <- build_analysis_s4(readcounts = readcounts, tx2gene =
#' gff3_file, metadata = metadata_df, config = config,
#'   tpm = tpm, effective_length = effective_length)
#' analysis <- filter_analysis_s4(analysis, min_samples = 1, subset_n_genes
#' = 200)
#' 
#' # Use setConfig via the method (called by setConfig generic)
#' new_config <- tsenat_config(q_values = c(0.5, 1.0, 1.5, 2.0), seed = 123)
#' analysis <- setConfig(analysis, new_config)
#' 
#' # Verify and display
#' summary(analysis)
setMethod("setConfig", "TSENATAnalysis", function(object, value) {
    # Convert TSENATConfig S4 object to list if needed
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

#' Set a single configuration value in TSENATAnalysis
#'
#' @param object \code{TSENATAnalysis} object.
#' @param key Character. Name of the config item to set.
#' @param value Value to set for the config item.
#'
#' @return Updated \code{TSENATAnalysis} object with the new config value.
#'
#' @details
#' Convenience method for setting a single configuration value without
#' needing to retrieve, merge, and set the entire config list. This preserves
#' all other configuration values while updating only the specified key.
#'
#' Unlike \code{\link{setConfig}}, which replaces the entire configuration,
#' \code{setConfigValue} performs a targeted update. It retrieves the current
#' config, updates one key-value pair, and stores the modified config back.
#'
#' @seealso
#' \code{\link{setConfig}} for replacing entire configuration,
#' \code{\link{getConfig}} for retrieving configuration,
#' \code{\link{tsenat_config}} for creating configuration objects
#'
#' @examples
#' # Load real TSENAT data
#' data(readcounts)
#' metadata_df <- read.table(system.file('extdata', 'metadata.tsv', package
#' = 'TSENAT'),
#'   header = TRUE, sep = '\t')
#' gff3_file <- system.file('extdata', 'annotation.gff3.gz', package = 'TSENAT')
#' 
#' # Build and subset analysis
#' config <- tsenat_config(sample_col = 'sample', condition_col = 'condition')
#' analysis <- build_analysis_s4(readcounts = readcounts, tx2gene =
#' gff3_file, metadata = metadata_df, config = config,
#'   tpm = tpm, effective_length = effective_length)
#' analysis <- filter_analysis_s4(analysis, min_samples = 1, subset_n_genes
#' = 200)
#' 
#' # Update a single configuration value while preserving others
#' analysis <- setConfigValue(analysis, 'q_values', c(0.5, 1.0, 1.5))
#' 
#' # Verify the update
#' config <- getConfig(analysis)
#' print(config$q_values)  # Shows c(0.5, 1.0, 1.5)
#'
#' @rdname TSENATAnalysis-methods
#' @export
setGeneric("setConfigValue", function(object, key, value) {
    standardGeneric("setConfigValue")
})

#' @rdname TSENATAnalysis-methods
#' @aliases setConfigValue,TSENATAnalysis-method
#' @keywords internal
#'
#' @examples
#' # Load real TSENAT data
#' data(readcounts)
#' metadata_df <- read.table(system.file('extdata', 'metadata.tsv', package
#' = 'TSENAT'),
#'   header = TRUE, sep = '\t')
#' gff3_file <- system.file('extdata', 'annotation.gff3.gz', package = 'TSENAT')
#' 
#' # Build and subset analysis with initial configuration
#' config <- tsenat_config(sample_col = 'sample', condition_col = 'condition')
#' analysis <- build_analysis_s4(readcounts = readcounts, tx2gene =
#' gff3_file, metadata = metadata_df, config = config,
#'   tpm = tpm, effective_length = effective_length)
#' analysis <- filter_analysis_s4(analysis, min_samples = 1, subset_n_genes
#' = 200)
#' 
#' # Use setConfigValue to update single configuration values
#' # This preserves all other config values
#' analysis <- setConfigValue(analysis, 'q_values', c(0.5, 1.5, 2.0))
#' analysis <- setConfigValue(analysis, 'seed', 456)
#' 
#' # Verify the updates
#' summary(analysis)
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
#' @param object \code{TSENATAnalysis} object.
#'
#' @return The \code{SummarizedExperiment} containing transcript/gene counts.
#'
#' @details
#' Provides type-safe accessor for the embedded \code{SummarizedExperiment}.
#'
#' @seealso
#' Other TSENATAnalysis accessors:  \code{\link{diversity}},
#'  \code{\link{divergence}},
#' \code{\link{lmResults}},
#'  \code{\link[S4Vectors]{metadata}}
#'
#' @examples
#' # Load real TSENAT data
#' data(readcounts)
#' metadata_df <- read.table(system.file('extdata', 'metadata.tsv', package
#' = 'TSENAT'),
#'   header = TRUE, sep = '\t')
#' gff3_file <- system.file('extdata', 'annotation.gff3.gz', package = 'TSENAT')
#' config <- tsenat_config(sample_col = 'sample', condition_col = 'condition')
#' analysis <- build_analysis_s4(readcounts = readcounts, tx2gene =
#' gff3_file, metadata = metadata_df, config = config,
#'   tpm = tpm, effective_length = effective_length)
#' analysis <- filter_analysis_s4(analysis, min_samples = 1, subset_n_genes
#' = 200)
#' @rdname TSENATAnalysis-methods
#' @export
setGeneric("se", function(object) {
    standardGeneric("se")
})

#' @noRd
if (!isGeneric("metadata")) {
    setGeneric("metadata", function(x, key = NULL) {
        standardGeneric("metadata")
    })
}

#' @rdname TSENATAnalysis-methods
#' @export
setMethod("se", "TSENATAnalysis", function(object) {
    object@se
})

#' Extract metadata from TSENATAnalysis
#'
#' @param x \code{TSENATAnalysis} object.
#' @param key \code{character} (optional). Specific metadata key to extract.
#'   If NULL, returns entire metadata list.
#'
#' @return The metadata list, or a specific metadata element if key is provided.
#'
#' @details
#' Provides type-safe accessor for analysis metadata (timestamps, function
#' calls,
#' intermediate results, etc.).
#'
#' @seealso
#' Other TSENATAnalysis accessors:  \code{\link{diversity}},
#'  \code{\link{divergence}},
#' \code{\link{lmResults}}, \code{\link{se}}
#'
#' @examples
#' # Load real TSENAT data
#' data(readcounts)
#' metadata_df <- read.table(system.file('extdata', 'metadata.tsv', package
#' = 'TSENAT'),
#'   header = TRUE, sep = '\t')
#' gff3_file <- system.file('extdata', 'annotation.gff3.gz', package = 'TSENAT')
#' analysis <- build_analysis_s4(readcounts = readcounts, tx2gene =
#' gff3_file, metadata = metadata_df,
#'   tpm = tpm, effective_length = effective_length)
#' analysis <- filter_analysis_s4(analysis, min_samples = 1, subset_n_genes
#' = 200)
#' @noRd
#' @rdname TSENATAnalysis-methods

setMethod("metadata", "TSENATAnalysis", function(x, key = NULL) {
    if (is.null(key)) {
        return(x@metadata)
    }

    if (key %in% names(x@metadata)) {
        return(x@metadata[[key]])
    }

    NULL
})

#' Test rank-based method assumptions
#'
#' @param analysis \code{TSENATAnalysis} object with diversity results.
#' @param q \code{numeric}. Q-value for diversity data to test.
#' @param checks \code{character}. Which assumptions to check.
#' @param alpha \code{numeric}. Significance level (default: 0.05).
#' @param ... Additional arguments.
#'
#' @return \code{TSENATAnalysis} object with stored assumption test results.
#'
#' @details
#' Tests whether rank-based analysis methods are appropriate for the data
#' by evaluating exchangeability, monotonicity, and consistency assumptions.
#'
#' @export
setGeneric("test_rankbased_assumptions_s4", function(analysis, q = NULL, checks = c("exchangeability",
    "monotonicity", "consistency"), alpha = 0.05, ...) {
    standardGeneric("test_rankbased_assumptions_s4")
})
