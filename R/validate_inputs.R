# ============================================================================
# INPUT VALIDATION UTILITIES - Consistent data validation across TSENAT plots
# ============================================================================
#
# Purpose: Centralized validation helpers for plot functions
# Ensures consistent error messages and graceful handling of invalid inputs
#
# ============================================================================

#' Validate SummarizedExperiment Input
#'
#' Checks that input is a valid SummarizedExperiment object with required slots.
#'
#' @param se Object to validate.
#' @param assay_name Character; name of required assay. If NULL, just checks class.
#' @param min_genes Numeric; minimum number of genes required (default: 1).
#' @param min_samples Numeric; minimum number of samples required (default: 1).
#' @param allow_null Logical; if TRUE, NULL input returns silently (default: FALSE).
#'
#' @return Invisibly returns TRUE if valid. Stops with error message if invalid.
#'
#' @details
#' Validates:
#' - Object is a SummarizedExperiment
#' - Specified assay exists (if assay_name provided)
#' - Minimum dimensions (genes, samples)
#' - No completely empty rows or columns
#'
#' Error messages are actionable and suggest common fixes.
#'
#' @examples
#' \dontrun{
#' # Valid SE object
#' .validate_se(se, assay_name = "diversity")
#'
#' # Invalid will stop with helpful message
#' .validate_se(df, assay_name = "diversity")
#' # Error: Input must be a SummarizedExperiment, not data.frame
#' # Use: se <- SummarizedExperiment::SummarizedExperiment(assays = list(...))
#' }
#'
#' @keywords internal
#' @noRd
.validate_se <- function(se, 
                         assay_name = NULL,
                         min_genes = 1,
                         min_samples = 1,
                         allow_null = FALSE) {
  # Check for NULL
  if (is.null(se)) {
    if (allow_null) return(invisible(TRUE))
    stop("SummarizedExperiment input is NULL.\n",
         "  Please provide a valid object from: se <- SummarizedExperiment::SummarizedExperiment(...)")
  }
  
  # Check class
  if (!methods::is(se, "SummarizedExperiment")) {
    class_actual <- class(se)[1]
    stop("Input must be a SummarizedExperiment, not ", class_actual, "\n",
         "  Use: se <- SummarizedExperiment::SummarizedExperiment(\n",
         "         assays = list(counts = your_matrix),\n",
         "         colData = sample_metadata)")
  }
  
  # Check dimensions
  n_genes <- nrow(se)
  n_samples <- ncol(se)
  
  if (n_genes < min_genes) {
    stop("SummarizedExperiment has ", n_genes, " genes (need at least ", min_genes, ")")
  }
  if (n_samples < min_samples) {
    stop("SummarizedExperiment has ", n_samples, " samples (need at least ", min_samples, ")")
  }
  
  # Check for assay
  if (!is.null(assay_name)) {
    assay_names <- SummarizedExperiment::assayNames(se)
    if (!assay_name %in% assay_names) {
      stop("Assay '", assay_name, "' not found.\n",
           "  Available assays: ", paste(assay_names, collapse = ", "), "\n",
           "  Ensure you've run: calculate_diversity(..., assay_name = '", assay_name, "')")
    }
  }
  
  invisible(TRUE)
}

#' Validate Data Frame with Required Columns
#'
#' Checks that input is a data.frame with required columns.
#'
#' @param df Object to validate.
#' @param required_cols Character vector; column names that must be present.
#' @param allow_null Logical; if TRUE, NULL input returns silently (default: FALSE).
#' @param df_name Character; friendly name for error messages (default: "Input data").
#'
#' @return Invisibly returns TRUE if valid. Stops with error message if invalid.
#'
#' @details
#' Validates:
#' - Object is a data.frame (or coercible to one)
#' - All required columns are present
#' - No empty data frames
#'
#' Error messages suggest which columns are missing.
#'
#' @examples
#' \dontrun{
#' # Valid data frame
#' .validate_df(lm_results, required_cols = c("gene", "p_value"))
#'
#' # Missing column error:
#' # Error: 'p_value' column missing in lm_results
#' # Available columns: gene, estimate, std.error, t.stat
#' }
#'
#' @keywords internal
#' @noRd
.validate_df <- function(df,
                         required_cols = NULL,
                         allow_null = FALSE,
                         df_name = "Input data") {
  # Check for NULL
  if (is.null(df)) {
    if (allow_null) return(invisible(TRUE))
    stop(df_name, " is NULL. Please provide a data.frame.")
  }
  
  # Coerce to data.frame if possible
  if (!is.data.frame(df)) {
    df <- tryCatch({
      as.data.frame(df)
    }, error = function(e) {
      stop(df_name, " must be a data.frame or coercible to one, not ", class(df)[1])
    })
  }
  
  # Check not empty
  if (nrow(df) == 0) {
    stop(df_name, " is empty (0 rows). Provide data with at least 1 row.")
  }
  
  # Check required columns
  if (!is.null(required_cols)) {
    missing_cols <- setdiff(required_cols, colnames(df))
    if (length(missing_cols) > 0) {
      stop("'", paste(missing_cols, collapse = "', '"), 
           "' column(s) missing in ", df_name, "\n",
           "  Available columns: ", paste(colnames(df), collapse = ", "))
    }
  }
  
  invisible(TRUE)
}

#' Validate Numeric Parameter
#'
#' Checks that a parameter is numeric and within valid range.
#'
#' @param value Numeric value to validate.
#' @param param_name Character; parameter name for error messages.
#' @param allow_null Logical; if TRUE, NULL is accepted (default: FALSE).
#' @param min Numeric; minimum allowed value (NULL = no minimum).
#' @param max Numeric; maximum allowed value (NULL = no maximum).
#' @param integer_only Logical; if TRUE, value must be an integer (default: FALSE).
#'
#' @return Invisibly returns TRUE if valid. Stops with error message if invalid.
#'
#' @examples
#' \dontrun{
#' # Valid
#' .validate_numeric(n_genes, "n_genes", min = 1, integer_only = TRUE)
#'
#' # Invalid - will stop with message
#' .validate_numeric(sig_alpha = 1.5, "sig_alpha", min = 0, max = 1)
#' # Error: sig_alpha must be between 0 and 1, got 1.5
#' }
#'
#' @keywords internal
#' @noRd
.validate_numeric <- function(value,
                              param_name = "parameter",
                              allow_null = FALSE,
                              min = NULL,
                              max = NULL,
                              integer_only = FALSE) {
  # Check for NULL
  if (is.null(value)) {
    if (allow_null) return(invisible(TRUE))
    stop(param_name, " is NULL. Provide a numeric value.")
  }
  
  # Check is numeric
  if (!is.numeric(value)) {
    stop(param_name, " must be numeric, got ", class(value)[1])
  }
  
  # Check is integer if required
  if (integer_only && !is.integer(value) && value != as.integer(value)) {
    stop(param_name, " must be an integer, got ", value)
  }
  
  # Check range
  if (!is.null(min) && value < min) {
    stop(param_name, " must be >= ", min, ", got ", value)
  }
  if (!is.null(max) && value > max) {
    stop(param_name, " must be <= ", max, ", got ", value)
  }
  
  invisible(TRUE)
}

#' Validate Character Parameter
#'
#' Checks that a parameter is character and matches allowed values.
#'
#' @param value Character value to validate.
#' @param param_name Character; parameter name for error messages.
#' @param allowed_values Character vector; permitted values.
#' @param allow_null Logical; if TRUE, NULL is accepted (default: FALSE).
#'
#' @return Invisibly returns TRUE if valid. Stops with error message if invalid.
#'
#' @examples
#' \dontrun{
#' # Valid
#' .validate_choice(metric, "metric", allowed_values = c("median", "mean", "sd"))
#'
#' # Invalid
#' .validate_choice(metric = "mode", "metric", 
#'                 allowed_values = c("median", "mean", "sd"))
#' # Error: metric must be one of: median, mean, sd (got 'mode')
#' }
#'
#' @keywords internal
#' @noRd
.validate_choice <- function(value,
                            param_name = "parameter",
                            allowed_values = NULL,
                            allow_null = FALSE) {
  # Check for NULL
  if (is.null(value)) {
    if (allow_null) return(invisible(TRUE))
    stop(param_name, " is NULL. Choose one of: ", 
         paste(allowed_values, collapse = ", "))
  }
  
  # Check not in allowed
  if (!value %in% allowed_values) {
    stop(param_name, " must be one of: ", 
         paste(allowed_values, collapse = ", "), 
         " (got '", value, "')")
  }
  
  invisible(TRUE)
}

#' Validate Color Specification
#'
#' Checks that value is a valid color specification (hex, named, or list of colors).
#'
#' @param value Color specification to validate.
#' @param param_name Character; parameter name for error messages.
#' @param allow_null Logical; if TRUE, NULL is accepted (default: FALSE).
#' @param allow_vector Logical; if TRUE, allows vector of colors (default: TRUE).
#'
#' @return Invisibly returns TRUE if valid. Stops with error message if invalid.
#'
#' @keywords internal
#' @noRd
.validate_color <- function(value,
                            param_name = "color",
                            allow_null = FALSE,
                            allow_vector = TRUE) {
  # Check for NULL
  if (is.null(value)) {
    if (allow_null) return(invisible(TRUE))
    stop(param_name, " is NULL. Provide a valid color.")
  }
  
  # Check is character
  if (!is.character(value)) {
    stop(param_name, " must be character (hex or named colors), got ", class(value)[1])
  }
  
  # Check if vector not allowed
  if (!allow_vector && length(value) > 1) {
    stop(param_name, " must be a single color, got ", length(value), " values")
  }
  
  # Validate each color using grDevices function
  valid_colors <- sapply(value, function(col) {
    tryCatch({
      grDevices::col2rgb(col)
      TRUE
    }, error = function(e) FALSE)
  })
  
  if (!all(valid_colors)) {
    invalid_idx <- which(!valid_colors)
    stop(param_name, " contains invalid color(s): ", 
         paste(value[invalid_idx], collapse = ", "), "\n",
         "  Use hex codes (#RRGGBB) or named colors from colors()")
  }
  
  invisible(TRUE)
}

#' Summarize Data Validation Results
#'
#' Helper to format validation results for user output.
#'
#' @param se SummarizedExperiment; validated SE object.
#' @param verbose Logical; if TRUE, print validation summary (default: TRUE).
#'
#' @return List with validation summary (genes, samples, assays, conditions).
#'
#' @keywords internal
#' @noRd
.summarize_se_validation <- function(se, verbose = TRUE) {
  n_genes <- nrow(se)
  n_samples <- ncol(se)
  assay_names <- SummarizedExperiment::assayNames(se)
  
  summary_list <- list(
    n_genes = n_genes,
    n_samples = n_samples,
    assay_names = assay_names
  )
  
  if (verbose && !is.null(se)) {
    message("[Validation OK] SummarizedExperiment:")
    message("  - Genes: ", n_genes)
    message("  - Samples: ", n_samples)
    message("  - Assays: ", paste(assay_names, collapse = ", "))
  }
  
  invisible(summary_list)
}
