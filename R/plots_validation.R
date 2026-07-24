# ============================================================================
# VALIDATION, FORMATTING & GENE FILTERING
# Extracted from plots_helpers.R — July 2026 refactoring (I11)
# ============================================================================


# ============================================================================
# DATA VALIDATION & QUALITY CHECKS
# ============================================================================

#' Validate Diversity SummarizedExperiment
#'
#' Checks that SE has required structure for diversity visualization.
#'
#' @param se A \code{SummarizedExperiment} object.
#' @param check_metadata Logical: also validate metadata? (default: TRUE)
#'
#' @return Logical TRUE if valid, else error with message.
#'

#' @noRd

.validate_diversity_se <- function(se, check_metadata = TRUE) {

    if (!inherits(se, "SummarizedExperiment")) {
        stop("se must be a SummarizedExperiment object", call. = FALSE)
    }

    if (nrow(se) == 0) {
        stop("SummarizedExperiment has no rows (samples)", call. = FALSE)
    }

    if (ncol(se) == 0) {
        stop("SummarizedExperiment has no columns (genes)", call. = FALSE)
    }

    # Check for diversity assay
    assay_names <- SummarizedExperiment::assayNames(se)
    if (!("diversity" %in% assay_names)) {
        stop("Required 'diversity' assay not found. ", "Available: ", paste(assay_names,
            collapse = ", "), call. = FALSE)
    }

    # Check for valid data
    div_mat <- SummarizedExperiment::assay(se, "diversity")
    if (all(is.na(div_mat))) {
        stop("All diversity values are NA", call. = FALSE)
    }

    if (check_metadata) {
        # Check for at least one q-value
        meta <- S4Vectors::metadata(se)
        if (!("q" %in% names(meta)) || length(meta$q) == 0) {
            warning("q-values not found in SE metadata", call. = FALSE)
        }
    }

    return(TRUE)
}

#' Validate Results Data Frame for Gene Selection
#'
#' Checks that results DataFrame has required columns.
#'
#' @param results Data frame (SAIT results, effect sizes, etc.).
#' @param require_pvalue Logical: check for p-value column? (default: TRUE)
#'
#' @return Logical TRUE if valid, else error.
#'

#' @noRd

.validate_results_df <- function(results, require_pvalue = TRUE) {

    if (!is.data.frame(results)) {
        stop("results must be a data frame", call. = FALSE)
    }

    if (nrow(results) == 0) {
        stop("results data frame is empty", call. = FALSE)
    }

    # Check for gene column
    gene_cols <- c("gene_id", "gene", "gene_name")
    has_gene <- any(gene_cols %in% colnames(results))
    if (!has_gene) {
        stop("No gene identifier column found. ", "Expected one of: ", paste(gene_cols,
            collapse = ", "), call. = FALSE)
    }

    # Check for p-value column
    if (require_pvalue) {
        p_cols <- c("adj_p_interaction", "p_interaction", "padj", "pvalue")
        has_pval <- any(p_cols %in% colnames(results))
        if (!has_pval) {
            stop("No p-value column found. ", "Expected one of: ", paste(p_cols,
                collapse = ", "), call. = FALSE)
        }
    }

    return(TRUE)
}

# ============================================================================
# FORMATTING & UTILITY FUNCTIONS
# ============================================================================

#' Format P-Value for Display
#'
#' Converts p-value to formatted string (scientific or threshold).
#'
#' @param pval Numeric p-value.
#' @param threshold Numeric: cutoff for '< threshold' format (default: 0.001).
#' @param digits Integer: decimal places for scientific notation (default: 2).
#'
#' @return Character string formatted p-value.
#'

#' @noRd

.format_pvalue <- function(pval, threshold = 0.001, digits = 2) {

    if (is.na(pval)) {
        return("NA")
    }

    if (pval < threshold) {
        return(paste0("< ", threshold))
    }

    return(format(pval, scientific = TRUE, digits = digits))
}

#' Format Q-Value Label
#'
#' Converts numeric q-value to display label (e.g., 'q = 1.0').
#'
#' @param q_val Numeric q-value.
#' @param prefix Character: prefix for label (default: 'q').
#'
#' @return Character string label.
#'

#' @noRd

.format_q_label <- function(q_val, prefix = "q") {
    if (is.na(q_val)) {
        return("NA")
    }
    return(sprintf("%s = %.2f", prefix, as.numeric(q_val)))
}

#' Format Label for Display
#'
#' @noRd

