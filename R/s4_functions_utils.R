#' Resolve Parameter Value with Priority: User > Config > Default
#'
#' @param user_value User-provided parameter value (checked first)
#' @param config_list Configuration list (typically analysis@config)
#' @param config_key Key name in config_list to check
#' @param default_value Default value if not found elsewhere
#' @param description Parameter name for error messages (optional)
#' @param allow_null Logical. If FALSE, stop if result is NULL. Default: TRUE.
#'
#' @return Resolved parameter value following priority: user > config > default
#'
#' @details
#' Priority resolution order:
#' 1. If \code{user_value} is not NULL, return it
#' 2. If \code{config_key} exists in \code{config_list}, return config value
#' 3. Return \code{default_value}
#'
#' This function centralizes the common pattern of parameter resolution
#' used throughout S4 wrapper functions.
#'
#' @examples
#' config <- list(nthreads = 4, verbose = TRUE, q_values = seq(0.01, 2, by =
#' 0.05))
#' 
#' # User-provided value takes priority
#' resolve_slot_param(user_value = 8, config_list = config, 
#'                    config_key = 'nthreads', default_value = 1)
#' # Returns: 8
#' 
#' # Falls back to config if user_value is NULL
#' resolve_slot_param(user_value = NULL, config_list = config, 
#'                    config_key = 'nthreads', default_value = 1)
#' # Returns: 4
#' 
#' # Uses default if neither user nor config
#' resolve_slot_param(user_value = NULL, config_list = config, 
#'                    config_key = 'nonexistent', default_value = 1)
#' # Returns: 1
#'

#' @noRd
resolve_slot_param <- function(user_value, config_list, config_key, default_value,
    description = NULL, allow_null = TRUE) {

    # Priority 1: User-provided value
    if (!is.null(user_value)) {
        return(user_value)
    }

    # Priority 2: Config value
    if (!is.null(config_list) && is.list(config_list) && config_key %in% names(config_list)) {
        config_val <- config_list[[config_key]]
        if (!is.null(config_val)) {
            return(config_val)
        }
    }

    # Priority 3: Default value
    if (!is.null(default_value)) {
        return(default_value)
    }

    # No value found
    if (!allow_null) {
        param_desc <- if (is.null(description))
            config_key else description
        stop("Could not resolve parameter '", param_desc, "'. ", "Please provide either: user argument, config['",
            config_key, "'], or update allow_null=TRUE", call. = FALSE)
    }

    return(NULL)
}


#' Auto-Detect Column from Candidate Names
#'
#' Search for a column by trying multiple candidate names with fallback options.
#'
#' @param available_cols Character vector of available column names
#' @param config_list Configuration list with optional priority column name
#' @param config_key Key name in config_list to check first
#' @param priority_candidates Character vector of candidate column names to
#' try in order
#' @param default_fallback Character string. Fallback if no matches found.
#' Default: NULL (no fallback).
#' @param verbose Logical. Print auto-detection messages. Default: FALSE.
#' @param param_name Parameter name for messages (e.g., 'condition_col').
#' Required for verbose.
#'
#' @return Character string with detected column name, or NULL if not found
#' and no fallback
#'
#' @details
#' Priority resolution order:
#' 1.  If \code{config_key} exists in \code{config_list},  check if 
#' it's in available columns
#' 2. Check each \code{priority_candidates} in order against available columns
#' 3. Return \code{default_fallback} if provided
#' 4. Return NULL otherwise
#'
#' This function centralizes the complex auto-detection logic used by multiple
#' S4 wrapper functions for finding condition_col, gene_col, isoform_col, etc.
#'
#' @examples
#' available <- c('sample_type', 'group', 'condition', 'other')
#' 
#' # Config takes priority
#' col <- auto_detect_column(available, list(col = 'sample_type'),
#'                           config_key = 'col', 
#'                           priority_candidates = c('condition', 'group'))
#' # Returns: 'sample_type'
#' 
#' # Falls back to priority candidates
#' col <- auto_detect_column(available, list(other = 'value'),
#'                           config_key = 'col', 
#'                           priority_candidates = c('condition', 'group'))
#' # Returns: 'group' (first match)
#'

#' @noRd
auto_detect_column <- function(available_cols, config_list = NULL, config_key = NULL,
    priority_candidates = character(0), default_fallback = NULL, verbose = FALSE,
    param_name = NULL) {

    # Priority 1: Check config
    if (!is.null(config_list) && !is.null(config_key) && config_key %in% names(config_list)) {
        cand <- config_list[[config_key]]
        if (!is.null(cand) && cand %in% available_cols) {
            return(cand)
        }
    }

    # Priority 2-N: Check priority candidates in order
    if (length(priority_candidates) > 0) {
        idx <- match(priority_candidates, available_cols)
        # Find first non-NA match
        first_match_idx <- which.min(is.na(idx))
        if (!is.na(first_match_idx) && first_match_idx <= length(idx) && !is.na(idx[first_match_idx])) {
            detected <- available_cols[idx[first_match_idx]]

            if (verbose && !is.null(param_name)) {
                message(sprintf("[auto_detect_column] Auto-detected %s = %s", param_name,
                  detected))
            }

            return(detected)
        }
    }

    # Fallback: return default_fallback
    if (!is.null(default_fallback)) {
        if (verbose && !is.null(param_name)) {
            message(sprintf("[auto_detect_column] Using fallback for %s = %s", param_name,
                default_fallback))
        }
        return(default_fallback)
    }

    return(NULL)
}


#' Centralized Output File Saving for Analysis Results
#'
#' Save analysis results to file with automatic format detection and handling.
#'
#' @param data Object to save (data.frame, matrix, ggplot,
#' SummarizedExperiment, or TSENATAnalysis)
#' @param output_file File path. Format auto-detected from extension (.rds,
#' .tsv, .csv, .txt, .pdf, .png, .jpg)
#' @param object Optional parent object (e.g., TSENATAnalysis) for
#' context-aware saving. Default: NULL.
#' @param verbose Logical. Print saving status. Default: FALSE.
#' @param create_dir Logical. Create directory if it doesn't exist. Default:
#' TRUE.
#' @param func_name Function name for error messages (e.g.,
#' 'calculate_diversity_s4'). Default: 'wrapper_function'
#' @param width Numeric or NULL. Plot width in inches (for ggplot/PDF/PNG
#' output). Default: NULL (use ggplot defaults).
#' @param height Numeric or NULL. Plot height in inches (for ggplot/PDF/PNG
#' output). Default: NULL (use ggplot defaults).
#'
#' @return Invisibly returns TRUE if successful, FALSE if error suppressed
#'
#' @details
#' **Supported formats:**
#' \itemize{
#'   \item \code{.rds}: R serialized object (for S4 objects, TSENATAnalysis)
#'   \item \code{. tsv,  . csv,  . txt}:  Tables (for  data. frame,  matrix,
#'  extracted tables)
#'   \item \code{.pdf}: PDF (for ggplot objects)
#'   \item \code{.png, .jpg}: Raster image (for ggplot objects)
#' }
#'
#' **Auto-detection logic:**
#' \itemize{
#'   \item ggplot object → saves as .pdf or .png/.jpg if extension matches
#'   \item data.frame → saves as .tsv/.csv/.txt
#'   \item SummarizedExperiment/matrix → converts to data.frame, saves as table
#'   \item TSENATAnalysis → saves as .rds (S4 object preservation)
#'   \item Other → attempts default save based on extension
#' }
#'
#' Errors in saving only warn (not stop) to prevent entire analysis failure
#' if output
#' file cannot be written.
#'
#' @examples
#' # Save a data frame
#' df <- data.frame(gene = c('GENE1', 'GENE2'), value = c(1.5, 2.3))
#' save_analysis_output(df, 'results/output.tsv')
#' 
#' # Save a ggplot object
#' p <- ggplot(data.frame(x = 1:10, y = rnorm(10)), aes(x, y)) + geom_point()
#' save_analysis_output(p, 'results/plot.pdf')
#'

#' @noRd
save_analysis_output <- function(data, output_file, object = NULL, verbose = FALSE,
    create_dir = TRUE, func_name = "wrapper_function", width = NULL, height = NULL) {

    if (is.null(output_file)) {
        return(invisible(FALSE))
    }

    # Create directory if requested
    if (create_dir) {
        output_dir <- dirname(output_file)
        if (output_dir != "." && !dir.exists(output_dir)) {
            tryCatch({
                dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
            }, error = function(e) {
                warning("[", func_name, "] Could not create output directory: ",
                  output_dir, " (", conditionMessage(e), ")", call. = FALSE)
            })
        }
    }

    # Determine format from extension
    ext <- tolower(sub("^.*\\.", ".", output_file))

    # Save based on format
    tryCatch({
        if (ext %in% c(".tsv", ".csv", ".txt")) {
            # Table format
            if (is.data.frame(data)) {
                write_data <- data
            } else if (inherits(data, "SummarizedExperiment")) {
                write_data <- as.data.frame(SummarizedExperiment::assay(data))
            } else if (is.matrix(data)) {
                write_data <- as.data.frame(data)
            } else {
                write_data <- as.data.frame(data)
            }

            # Auto-detect separator
            sep <- if (ext == ".tsv")
                "\t" else ifelse(ext == ".csv", ",", "\t")
            write.table(write_data, file = output_file, sep = sep, quote = FALSE,
                row.names = TRUE, col.names = NA)

            if (verbose) {
                message("[", func_name, "] Saved table to: ", output_file)
            }

        } else if (ext %in% c(".pdf", ".png", ".jpg", ".jpeg")) {
            # Plot format (ggplot)
            if (!requireNamespace("ggplot2", quietly = TRUE)) {
                warning("[", func_name, "] ggplot2 required for plot output", call. = FALSE)
                return(invisible(FALSE))
            }

            if (!inherits(data, "ggplot")) {
                warning("[", func_name, "] Cannot save non-ggplot object as ", ext,
                  ". Skipping plot output.", call. = FALSE)
                return(invisible(FALSE))
            }

            # Use provided dimensions or defaults
            plot_width <- if (is.null(width))
                8 else width
            plot_height <- if (is.null(height))
                6 else height

            if (ext == ".pdf") {
                grDevices::pdf(output_file, width = plot_width, height = plot_height)
                print(data)
                grDevices::dev.off()
            } else {
                grDevices::png(output_file, width = plot_width, height = plot_height,
                  units = "in", res = 300)
                print(data)
                grDevices::dev.off()
            }

            if (verbose) {
                message("[", func_name, "] Saved plot to: ", output_file)
            }

        } else if (ext == ".rds") {
            # R serialized object
            saveRDS(data, file = output_file)

            if (verbose) {
                message("[", func_name, "] Saved RDS object to: ", output_file)
            }

        } else {
            # Unknown format
            warning("[", func_name, "] Unknown output format: ", ext, ". Trying default RDS save.",
                call. = FALSE)
            saveRDS(data, file = output_file)
        }

        return(invisible(TRUE))

    }, error = function(e) {
        warning("[", func_name, "] Could not save to ", output_file, ": ", conditionMessage(e),
            call. = FALSE)
        return(invisible(FALSE))
    })
}


#' Extract Multi-Q Tabular Results with Q-Value Column
#'
#' Helper to combine results from multiple q-values into a single data.frame
#' with q_value column.
#'
#' @param result Result object (may be list with multi-q structure or single
#' result)
#' @param is_multiq Logical. Whether result has multi-q structure. Default:
#' auto-detect.
#' @param extract_fn Function to extract table from each result element.
#'   Signature:  \code{function(result_element,  q_key)}.  Default:
#'  extracts 'summary_table'.
#' @param q_value_col Name of q-value column to add. Default: 'q_value'
#'
#' @return data.frame with combined results and q_value column
#'
#' @examples
#' # Create sample multi-q result structure
#' q_result <- list(
#'   q_0.5 = list(summary_table = data.frame(
#'     gene_id = c('GENE1', 'GENE2'),
#'     entropy = c(1.2, 1.5),
#'     psi_mean = c(0.3, 0.7)
#'   )),
#'   q_1.0 = list(summary_table = data.frame(
#'     gene_id = c('GENE1', 'GENE2'),
#'     entropy = c(1.1, 1.4),
#'     psi_mean = c(0.32, 0.68)
#'   ))
#' )
#' 
#' # Extract and combine results across q-values
#' combined_results <- extract_multiq_table(q_result)
#' head(combined_results)
#' 
#' @export
extract_multiq_table <- function(result, is_multiq = NULL, extract_fn = NULL, q_value_col = "q_value") {

    # Auto-detect multi-q if not specified
    if (is.null(is_multiq)) {
        is_multiq <- inherits(result, "tsenat_isoform_switching_multiq") || (is.list(result) &&
            any(grepl("^q_", names(result))))
    }

    # Default extraction function
    if (is.null(extract_fn)) {
        extract_fn <- function(result_element, q_key) {
            if (!is.null(result_element$summary_table)) {
                return(result_element$summary_table)
            }
            return(NULL)
        }
    }

    if (!is_multiq) {
        # Single q-value: return as-is or wrap
        table <- extract_fn(result, "q_1_00")
        if (!is.null(table) && is.data.frame(table)) {
            return(table)
        }
        return(NULL)
    }

    # Multi-q case: combine results
    combined <- do.call(rbind, lapply(names(result), function(q_key) {
        q_result <- result[[q_key]]
        table <- extract_fn(q_result, q_key)

        if (!is.null(table) && is.data.frame(table)) {
            table[[q_value_col]] <- sub("^q_", "", q_key)
            return(table)
        }
        return(NULL)
    }))

    if (!is.null(combined)) {
        rownames(combined) <- NULL  # Reset row names
    }

    return(combined)
}
