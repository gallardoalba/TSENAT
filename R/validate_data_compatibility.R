#' Validate Gene Name Alignment Between Objects
#'
#' Checks that gene identifiers are consistent between two data objects
#' (e.g., SummarizedExperiment and sait_results data.frame).
#'
#' @param se \code{SummarizedExperiment} or  \code{NULL}.  Source object with 
#' rownames as gene IDs.
#' @param results \code{data. frame} or  \code{NULL}.  Results object with 
#' rownames as gene IDs.
#' @param se_name \code{character}.  Display name for  first object (e. g. ,
#'  'SummarizedExperiment').
#' @param results_name \code{character}.  Display name for  second object (e.
#' g. ,  'SAIT results').
#' @param verbose \code{logical}.  If TRUE,  print alignment report (default:
#'  TRUE).
#'
#' @return List with validation results:
#' \itemize{
#'   \item \code{is_aligned}: logical, TRUE if genes match
#'   \item \code{n_se_genes}: number of genes in SE
#'   \item \code{n_results_genes}: number of genes in results
#'   \item \code{n_matched}: number of genes in both
#'   \item \code{n_se_only}: genes only in SE
#'   \item \code{n_results_only}: genes only in results
#'   \item \code{mismatch_details}: character vector with specifics
#' }
#'
#' @details
#' Performs strict matching: rownames(se) vs rownames(results).
#' Flags mismatches and provides suggestion for resolution.
#'
#' @examples
#' library(SummarizedExperiment)
#' tx_counts <- matrix(sample(10:100, 400, replace = TRUE), nrow = 40, ncol
#' = 10,
#'   dimnames = list(paste0('TX', 1:40), paste0('Sample', 1:10)))
#' se <- SummarizedExperiment(assays = list(counts = tx_counts))
#' S4Vectors::metadata(se)$tx2gene <- data.frame(
#'   Transcript = paste0('TX', 1:40), Gen = rep(paste0('GENE', 1:10), each = 4))
#' sait_results <- data.frame(gene = paste0('GENE', 1:10), p_value = runif(10))
#' validation <- .validate_gene_names(
#'   se, sait_results, se_name = 'Input SE', results_name = 'LM Results'
#' )
#'
#' @noRd

.validate_gene_names <- function(se = NULL, results = NULL, se_name = "Object1",
    results_name = "Object2", verbose = TRUE) {

    # Extract gene identifiers
    se_genes <- if (!is.null(se) && methods::is(se, "SummarizedExperiment")) {
        rownames(se)
    } else {
        NULL
    }

    # Extract results genes - try column first, then rownames
    results_genes <- NULL
    if (!is.null(results) && is.data.frame(results)) {
        # Try to get genes from column first
        if ("gene" %in% colnames(results)) {
            results_genes <- results$gene
        } else if ("gene_id" %in% colnames(results)) {
            results_genes <- results$gene_id
        } else if ("gene_name" %in% colnames(results)) {
            results_genes <- results$gene_name
        } else {
            # Fallback to rownames
            results_genes <- rownames(results)
        }
    } else if (!is.null(results) && is.matrix(results)) {
        results_genes <- rownames(results)
    }

    # Basic validation
    if (is.null(se_genes) || length(se_genes) == 0) {
        return(list(is_aligned = FALSE, error = paste0(se_name, " has no rownames")))
    }

    if (is.null(results_genes) || length(results_genes) == 0) {
        return(list(is_aligned = FALSE, error = paste0(results_name, " has no rownames")))
    }

    # Compare gene sets
    n_se <- length(se_genes)
    n_results <- length(results_genes)
    matched <- intersect(se_genes, results_genes)
    n_matched <- length(matched)
    se_only <- setdiff(se_genes, results_genes)
    results_only <- setdiff(results_genes, se_genes)

    is_aligned <- identical(sort(se_genes), sort(results_genes))

    # Build report
    mismatch_details <- c()

    if (!is_aligned) {
        mismatch_details <- c(paste0("Gene name mismatch between ", se_name, " and ",
            results_name), paste0("  ", se_name, ": ", n_se, " genes"), paste0("  ",
            results_name, ": ", n_results, " genes"), paste0("  Matched: ", n_matched,
            " genes"), paste0("  Only in ", se_name, ": ", length(se_only), " genes"),
            paste0("  Only in ", results_name, ": ", length(results_only), " genes"))

        if (length(se_only) > 0 && length(se_only) <= 5) {
            mismatch_details <- c(mismatch_details, paste0("    Examples: ", paste(head(se_only,
                3), collapse = ", ")))
        }

        if (length(results_only) > 0 && length(results_only) <= 5) {
            mismatch_details <- c(mismatch_details, paste0("    Examples: ", paste(head(results_only,
                3), collapse = ", ")))
        }
    }

    if (verbose) {
        if (is_aligned) {
            message("[validate_gene_names] Gene names aligned perfectly")
            message("  ", n_se, " genes in both ", se_name, " and ", results_name)
        } else {
            message("[validate_gene_names] Gene name mismatch!")
            for (line in mismatch_details) message(line)
            message("\nSOLUTION: Ensure both objects were derived from the same analysis.")
        }
    }

    return(list(is_aligned = is_aligned, n_se_genes = n_se, n_results_genes = n_results,
        n_matched = n_matched, n_se_only = length(se_only), n_results_only = length(results_only),
        mismatch_details = paste(mismatch_details, collapse = "\n")))
}


#' Validate Data Dimension Compatibility
#'
#' Checks that assay dimensions and colData align across objects.
#'
#' @param se \code{SummarizedExperiment}. Object to validate.
#' @param expected_n_samples \code{integer} or  \code{NULL}.
#'  Expected number of samples.
#'   If NULL, no check performed (default: NULL).
#' @param expected_assays \code{character} vector or  \code{NULL}.
#'  Expected assay names
#'   (e.g., c('diversity', 'entropy')). If NULL, no check performed.
#' @param verbose \code{logical}.  If TRUE,  print validation report (default:
#'  TRUE).
#'
#' @return List with validation results:
#' \itemize{
#'   \item \code{is_valid}: logical, TRUE if all checks pass
#'   \item \code{n_samples}: actual number of samples
#'   \item \code{n_genes}: actual number of genes
#'   \item \code{assays_present}: character vector of assay names
#'   \item \code{issues}: character vector of problems found
#' }
#'
#' @noRd

.validate_se_dimensions <- function(se, expected_n_samples = NULL, expected_assays = NULL,
    verbose = TRUE) {

    if (!methods::is(se, "SummarizedExperiment")) {
        stop("'se' must be a SummarizedExperiment object", call. = FALSE)
    }

    issues <- c()
    n_genes <- nrow(se)
    n_samples <- ncol(se)
    assays_present <- names(SummarizedExperiment::assays(se))

    # Check sample count
    if (!is.null(expected_n_samples)) {
        if (n_samples != expected_n_samples) {
            issues <- c(issues, paste0("Sample count mismatch: expected ", expected_n_samples,
                ", got ", n_samples))
        }
    }

    # Check assays
    if (!is.null(expected_assays)) {
        missing_assays <- setdiff(expected_assays, assays_present)
        if (length(missing_assays) > 0) {
            issues <- c(issues, paste0("Missing assays: ", paste(missing_assays,
                collapse = ", "), ". Available: ", paste(assays_present, collapse = ", ")))
        }
    }

    # Check colData consistency
    cd <- SummarizedExperiment::colData(se)
    if (nrow(cd) != n_samples) {
        issues <- c(issues, paste0("colData row count (", nrow(cd), ") != ncol(se) (",
            n_samples, ")"))
    }

    # Check rowData consistency
    rd <- SummarizedExperiment::rowData(se)
    if (!is.null(rd) && nrow(rd) > 0 && nrow(rd) != n_genes) {
        issues <- c(issues, paste0("rowData row count (", nrow(rd), ") != nrow(se) (",
            n_genes, ")"))
    }

    is_valid <- length(issues) == 0

    if (verbose) {
        if (is_valid) {
            message("[validate_se_dimensions] All dimension checks passed")
            message("  Samples: ", n_samples, ", Genes: ", n_genes)
            if (length(assays_present) > 0) {
                message("  Assays: ", paste(assays_present, collapse = ", "))
            }
        } else {
            message("[validate_se_dimensions] Dimension mismatches found!")
            for (issue in issues) message("  ", issue)
        }
    }

    return(list(is_valid = is_valid, n_samples = n_samples, n_genes = n_genes, assays_present = assays_present,
        issues = issues))
}


#' Validate LM Results Data Structure
#'
#' Checks that SAIT results object has expected structure and content.
#'
#' @param sait_results \code{data.frame} or \code{list}. SAIT results to validate.
#' @param expected_genes \code{character} vector or  \code{NULL}.
#'  Gene names expected
#'   in results (default: NULL, no check).
#' @param required_columns \code{character} vector.  Column names that 
#' must be present
#' (default: c('gene', 'p_interaction', 'effect_size')). For flexibility
#' with different
#'   LM implementations, only ONE of the p-value columns needs to be present:
#'   'p_value', 'p_interaction', 'p_raw', or 'adj_p_interaction'.
#' @param verbose \code{logical}.  If TRUE,  print validation report (default:
#'  TRUE).
#'
#' @return List with validation results:
#' \itemize{
#'   \item \code{is_valid}: logical, TRUE if structure is valid
#'   \item \code{n_results}: number of result rows
#'   \item \code{columns_present}: character vector of column names
#'   \item \code{issues}: character vector of problems found
#' }
#'
#' @noRd

.validate_sait_results <- function(sait_results, expected_genes = NULL, required_columns = c("gene"),
    verbose = TRUE) {

    issues <- c()

    # Check if it's a data.frame
    if (!is.data.frame(sait_results)) {
        if (is.list(sait_results)) {
            # Could be list with $results and $model_data
            if ("results" %in% names(sait_results)) {
                sait_results <- sait_results$results
            } else {
                issues <- c(issues, "SAIT results must be data.frame or list with 'results' element")
                return(list(is_valid = FALSE, n_results = 0, columns_present = c(),
                  issues = issues))
            }
        } else {
            issues <- c(issues, "SAIT results must be data.frame or list")
            return(list(is_valid = FALSE, n_results = 0, columns_present = c(), issues = issues))
        }
    }

    n_results <- nrow(sait_results)
    cols_present <- colnames(sait_results)

    # Check basic required columns (only 'gene' is truly required, but
    # flexible)
    for (col in required_columns) {
        if (!(col %in% cols_present) && col != "p_value") {
            # 'gene' column not strictly required if genes in rownames
            if (col == "gene" && !is.null(rownames(sait_results))) {
                # OK - genes are in rownames
            } else if (col == "gene") {
                issues <- c(issues, paste0("Missing gene identifier. Need 'gene' column or genes in rownames. Available: ",
                  paste(cols_present, collapse = ", ")))
            }
        }
    }

    # Check for at least one p-value column (flexible acceptance)
    pval_columns <- c("p_value", "p_interaction", "p_raw", "adj_p_interaction", "padj",
        "p")
    has_pval <- any(pval_columns %in% cols_present)

    if (!has_pval) {
        issues <- c(issues, paste0("Missing p-value column. Need at least one of: ",
            paste(pval_columns, collapse = ", "), ". Available: ", paste(cols_present,
                collapse = ", ")))
    }

    # Check for effect size or coefficient column (if available)
    effect_cols <- c("coef", "effect_size", "log2FoldChange", "logFC", "beta")
    has_effect <- any(effect_cols %in% cols_present)
    # Note: effect size/coef is nice to have, but not required for validation

    # Check gene names if expected
    if (!is.null(expected_genes)) {
        # Try to get gene names from column or rownames
        result_genes <- NULL
        if ("gene" %in% cols_present || "gene_id" %in% cols_present || "gene_name" %in%
            cols_present) {
            gene_col <- intersect(c("gene", "gene_id", "gene_name"), cols_present)[1]
            result_genes <- sait_results[[gene_col]]
        } else if (!is.null(rownames(sait_results))) {
            result_genes <- rownames(sait_results)
        }

        if (!is.null(result_genes)) {
            unmatched <- setdiff(expected_genes, result_genes)
            if (length(unmatched) > 0) {
                issues <- c(issues, paste0("Gene mismatch: ", length(unmatched),
                  "/", length(expected_genes), " genes not in results"))
            }
        }
    }

    # Check for p-value validity (only if column exists)
    pval_col <- intersect(pval_columns, cols_present)[1]
    if (!is.na(pval_col) && pval_col %in% cols_present) {
        n_na <- sum(is.na(sait_results[[pval_col]]))
        # Note: having NAs is OK for partial results, just informational
        if (n_na > 0 && verbose) {
            # Don't add to issues - just note it
        }

        invalid_pvals <- sum(sait_results[[pval_col]] < 0 | sait_results[[pval_col]] >
            1, na.rm = TRUE)
        if (invalid_pvals > 0) {
            issues <- c(issues, paste0("Invalid p-values in '", pval_col, "' (",
                invalid_pvals, " outside [0,1])"))
        }
    }

    is_valid <- length(issues) == 0

    if (verbose) {
        if (is_valid) {
            message("[validate_sait_results] SAIT results structure valid")
            message("  Results: ", n_results, " genes")
            message("  Columns: ", paste(cols_present, collapse = ", "))
            if (has_pval)
                message("  P-value column: ", pval_col)
            if (has_effect)
                message("  Effect size: available")
        } else {
            message("[validate_sait_results] SAIT results structure issues!")
            for (issue in issues) message("  ", issue)
        }
    }

    return(list(is_valid = is_valid, n_results = n_results, columns_present = cols_present,
        issues = issues))
}


#' Check All Data Compatibility for Plotting
#'
#' Comprehensive validation for plotting functions.
#' Runs all alignment checks between SE, SAIT results, and metadata.
#'
#' @param se \code{SummarizedExperiment}. Input data object.
#' @param sait_results \code{data.frame} or \code{NULL}. SAIT interaction results.
#' @param stop_on_error \code{logical}.  If TRUE,
#'  stop execution on validation failure.
#'   If FALSE, return issues and continue (default: TRUE).
#' @param verbose \code{logical}. Print validation report (default: TRUE).
#'
#' @return Invisible list with all validation results. Returns NULL if all
#' checks pass.
#'   If validation fails and stop_on_error=TRUE, throws error before returning.
#'
#' @details
#' Runs in sequence:
#' 1. SE dimension check
#' 2. SAIT results structure check
#' 3. Gene name alignment between SE and SAIT results
#'
#' If any check fails, provides consolidated error message.
#'
#' @noRd

.validate_plot_data <- function(se, sait_results = NULL, stop_on_error = TRUE, verbose = TRUE) {

    all_issues <- list()

    # 1. Check SE dimensions
    se_check <- .validate_se_dimensions(se, verbose = verbose)
    if (!se_check$is_valid) {
        all_issues$se_dimensions <- se_check$issues
    }

    # 2. Check SAIT results if provided
    if (!is.null(sait_results)) {
        sait_check <- .validate_sait_results(sait_results, verbose = verbose)
        if (!sait_check$is_valid) {
            all_issues$sait_results <- sait_check$issues
        }

        # 3. Check gene alignment
        if (nrow(se) > 0 && nrow(sait_results) > 0) {
            gene_check <- .validate_gene_names(se = se, results = sait_results, se_name = "SummarizedExperiment",
                results_name = "LM Results", verbose = verbose)
            if (!gene_check$is_aligned) {
                all_issues$gene_alignment <- gene_check$mismatch_details
            }
        }
    }

    # Report consolidated results
    if (length(all_issues) > 0) {
        error_msg <- "Data validation failed:\n\n"
        for (check_name in names(all_issues)) {
            error_msg <- paste0(error_msg, "  [", check_name, "]\n")
            for (issue in all_issues[[check_name]]) {
                error_msg <- paste0(error_msg, "    ", issue, "\n")
            }
        }

        if (stop_on_error) {
            stop(error_msg, call. = FALSE)
        } else {
            warning(error_msg, call. = FALSE)
        }
    }

    if (verbose && length(all_issues) == 0) {
        message("[validate_plot_data] All data compatibility checks passed!")
    }

    return(invisible(all_issues))
}
