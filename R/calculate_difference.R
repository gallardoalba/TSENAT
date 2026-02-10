#' Calculate splicing diversity changes between two conditions.
#' @param x A \code{SummarizedExperiment} with splicing diversity values for
#' each gene in each sample or a \code{data.frame} with gene names in the  first
#' column and splicing diversity values for each sample in additional  columns.
#' @param samples A vector of length one, specifying the column name of the
#' \code{colData} annotation column from the \code{SummarizedExperiment}
#' object, that should be used as the category column or a character vector
#' with an equal length to the number of columns in the input dataset,
#' specifying the category of each sample in the case of a \code{data.frame}
#' input.
#' @param control Name of the control sample category, defined in the
#' \code{samples} vector, e.g. \code{control = 'Normal'} or \code{control =
#' 'WT'}.
#' @param method Method to use for calculating the average splicing diversity
#' value in a condition. Can be \code{'mean'} or \code{'median'}.
#' @param test Method to use for p-value calculation: use \code{'wilcoxon'} for
#' Wilcoxon rank sum test or \code{'shuffle'} for a label shuffling test.
#' @param randomizations Number of random shuffles, used for the label shuffling
#' test (default = 100).
#' @param pcorr P-value correction method applied to the Wilcoxon rank sum test
#' or label shuffling test results, as defined in the \code{p.adjust}  function.
#' @param assayno An integer value. In case of multiple assays in a
#' \code{SummarizedExperiment} input, the argument specifies the assay number
#' to use for difference calculations.
#' @param verbose If \code{TRUE}, the function will print additional diagnostic
#' messages.
#' @param pseudocount Numeric scalar. Passed to \code{calculate_fc} and used to
#' add a small value to non-positive group summaries before computing
#' differences and log2 fold-changes. Default \code{1e-6}. Rows excluded for
#' low sample counts remain \code{NA}.
#' @param paired Logical; if `TRUE`, run paired versions of tests when
#'   supported (default: `FALSE`).
#' @param exact Logical; passed to the Wilcoxon test to request exact p-values
#'   when supported (default: `FALSE`).
#' @return A \code{data.frame} with the mean or median values of splicing
#' diversity across sample categories and all samples, log2(fold change) of  the
#' two different conditions, raw and corrected p-values.
#' @import methods
#' @importFrom SummarizedExperiment SummarizedExperiment assays assay colData
#' @export
#' @details The function calculates diversity changes between two sample
#' conditions. It uses the output of the diversity calculation function, which
#' is a \code{SummarizedExperiment} object of splicing diversity values.
#' Additionally, it can use a \code{data.frame} as input, where the first column
#' contains gene names, and all additional columns contain splicing diversity
#' values for each sample. A vector of sample conditions also serves as input,
#' used for aggregating the samples by condition.   It calculates the mean or
#' median of the splicing diversity data per sample  condition, the difference
#' of these values and the log2 fold change of the two  conditions. Furthermore,
#' the user can select a statistical method to  calculate the significance of
#' the changes. The p-values and adjusted p-values  are calculated using a
#' Wilcoxon sum rank test or label shuffling test.   The function will exclude
#' genes of low sample size from the significance  calculation, depending on
#' which statistical test is applied.
#' @examples
#' x <- data.frame(Genes = letters[seq_len(10)], matrix(runif(80), ncol = 8))
#' samples <- c(rep("Healthy", 4), rep("Pathogenic", 4))
#' calculate_difference(x, samples,
#'     control = "Healthy", method = "mean", test =
#'         "wilcoxon"
#' )
calculate_difference <- function(x, samples = NULL, control, method = "mean", test = "wilcoxon",
  randomizations = 100, pcorr = "BH", assayno = 1, verbose = TRUE, paired = FALSE,
  exact = FALSE, pseudocount = 0) {
    # internal small helpers (kept here to avoid adding new files)
    .tsenat_prepare_df <- function(x, samples, assayno) {
        if (inherits(x, "RangedSummarizedExperiment") || inherits(x, "SummarizedExperiment")) {
            # allow samples to be NULL (use default 'sample_type' col)
            if (is.null(samples)) {
                if ("sample_type" %in% colnames(SummarizedExperiment::colData(x))) {
                    samples_col <- "sample_type"
                } else {
                    stop("When providing a SummarizedExperiment, supply 'samples' as a colData column name or call map_metadata() to populate 'sample_type'",
                        call. = FALSE)
                }
            } else {
                if (length(samples) != 1) {
                    stop("'samples' must be a single colData column.", call. = FALSE)
                }
                samples_col <- samples
            }
            samples_vec <- SummarizedExperiment::colData(x)[[samples_col]]
            if (!is.numeric(assayno) || length(SummarizedExperiment::assays(x)) <
                assayno) {
                stop("Invalid 'assayno'.", call. = FALSE)
            }
            df <- as.data.frame(SummarizedExperiment::assays(x)[[assayno]])
            genes <- rownames(df)
            df <- cbind(genes = genes, df)
            list(df = df, samples = samples_vec)
        } else {
            df <- as.data.frame(x)
            list(df = df, samples = samples)
        }
    }

    .tsenat_sample_matrix <- function(dfr) {
        as.matrix(dfr[, -c(1, ncol(dfr) - 1, ncol(dfr)), drop = FALSE])
    }

    # Validate input container Reject matrices explicitly (tests expect this
    # error for matrix input)
    if (is.matrix(x)) {
        stop("Input type unsupported; see ?calculate_difference.", call. = FALSE)
    }
    if (!(is.data.frame(x) || inherits(x, "RangedSummarizedExperiment") || inherits(x,
        "SummarizedExperiment"))) {
        stop("Input data type not supported; see ?calculate_difference.", call. = FALSE)
    }

    # prepare data.frame and sample vector (handles SummarizedExperiment)
    pd <- .tsenat_prepare_df(x, samples, assayno)
    df <- pd$df
    samples <- pd$samples

    # Partition and validate inputs
    part <- .tsenat_calculate_difference_partition(df = df, samples = samples, control = control, method = method, test = test, pcorr = pcorr, randomizations = randomizations, verbose = verbose)
    df <- part$df
    samples <- part$samples
    groups <- part$groups
    idx1 <- part$idx1
    idx2 <- part$idx2
    df_keep <- part$df_keep
    df_small <- part$df_small

    result_list <- list()

    # Helper to extract the numeric matrix of sample columns (keeps original
    # order)
    sample_matrix <- .tsenat_sample_matrix

    if (nrow(df_keep) > 0) {
        if (nrow(df_small) > 0 && verbose) {
            message(sprintf("Note: %d genes excluded due to low sample counts.",
                nrow(df_small)))
        }
        ymat <- sample_matrix(df_keep)
        # p-value calculation
        if (test == "wilcoxon") {
            ptab <- wilcoxon(ymat, samples, pcorr = pcorr, paired = paired, exact = exact)
            # wilcoxon should return a data.frame of p-values named
            # appropriately
        } else {
            ptab <- label_shuffling(ymat, samples, control, method, randomizations = randomizations,
                pcorr = pcorr, paired = paired)
        }
        result_list$tested <- data.frame(genes = df_keep[, 1], calculate_fc(ymat,
            samples, control, method, pseudocount = pseudocount), ptab, stringsAsFactors = FALSE)
    }

    if (nrow(df_small) > 0) {
        small_mat <- sample_matrix(df_small)
        result_list$small <- data.frame(genes = df_small[, 1], calculate_fc(small_mat,
            samples, control, method, pseudocount = pseudocount), raw_p_values = NA,
        adjusted_p_values = NA, stringsAsFactors = FALSE)
    }

    # Combine results preserving tested rows first
    if (length(result_list) == 0) {
        return(data.frame())
    }
    res <- do.call(rbind, result_list)
    rownames(res) <- NULL
    if (!("log2_fold_change" %in% colnames(res)) && ("log2FC" %in% colnames(res))) {
        res$log2_fold_change <- res$log2FC
    }
    return(res)
}


#' Linear-model interaction test for Tsallis entropy   For each gene, fit a
#' linear model of the form `entropy ~ q * group` and  extract the p-value for
#' the interaction term (whether the effect of `q`  differs between groups). The
#' function expects a `SummarizedExperiment`  produced by
#' `calculate_diversity()` when multiple `q` values have been  computed (column
#' names contain `_q=`).
#' @param se A `SummarizedExperiment` containing a `diversity` assay produced
#' by `calculate_diversity(..., q = <vector>)`.
#' @param sample_type_col Optional column name in `colData(se)` that contains
#' a grouping factor for samples (character). If `NULL`, the function will
#' attempt to infer group from column names (suffix `_N` interpreted as
#' 'Normal').
#' @param min_obs Minimum number of non-NA observations required to fit a
#' model for a gene (default: 10).
#' @param method Modeling method to use for interaction testing: one of
#' \code{c('linear', 'lmm', 'gam', 'fpca')} (default: 'linear').
#' @param pvalue Type of p-value to compute for linear mixed models: one of
#' \code{c('satterthwaite', 'lrt', 'both')} (default: 'satterthwaite').
#' @param subject_col Optional column name in `colData(se)` that contains
#' subject/individual identifiers for paired or repeated-measures designs
#' (character). If provided with `method = 'lmm'`, used as random effect.
#' @param paired Logical; whether samples are paired (default: FALSE).
#' @param nthreads Number of threads (mc.cores) to use for parallel processing
#' (default: 1).
#' @param assay_name Name of the assay in the SummarizedExperiment to use
#' (default: 'diversity').
#' @param verbose Logical; whether to print progress messages during execution
#' (default: FALSE).
#' @return A data.frame with columns `gene`, `p_interaction`, and
#' `adj_p_interaction`, ordered by ascending `p_interaction`.
#' @export
#' @examples
#' data("tcga_brca_luma_dataset", package = "TSENAT")
#' rc <- as.matrix(tcga_brca_luma_dataset[1:20, -1, drop = FALSE])
#' gs <- tcga_brca_luma_dataset$genes[1:20]
#' se <- calculate_diversity(rc, gs, q = c(0.1, 1), norm = TRUE)
#' # Provide a minimal sample-type mapping so the example runs during checks
#' SummarizedExperiment::colData(se) <- S4Vectors::DataFrame(
#'     sample_type = rep(c("Normal", "Tumor"), length.out = ncol(se)),
#'     row.names = colnames(se)
#' )
#' calculate_lm_interaction(se, sample_type_col = "sample_type")
calculate_lm_interaction <- function(se, sample_type_col = NULL, min_obs = 10, method = c("linear",
      "lmm", "gam", "fpca"), pvalue = c("satterthwaite", "lrt", "both"), subject_col = NULL,
  paired = FALSE, nthreads = 1, assay_name = "diversity", verbose = FALSE) {
    method <- match.arg(method)
    pvalue <- match.arg(pvalue)
    if (verbose) {
        message("[calculate_lm_interaction] method=", method)
    }
    # internal flags: keep these internal to avoid documenting them in Rd
    suppress_lme4_warnings <- TRUE
    progress <- FALSE
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
        stop("SummarizedExperiment required")
    }

    mat <- SummarizedExperiment::assay(se, assay_name)
    if (is.null(mat)) {
        stop(sprintf("Assay '%s' not found in SummarizedExperiment", assay_name))
    }

    sample_q <- colnames(mat)
    if (is.null(sample_q) || length(sample_q) == 0) {
        stop("No column names found on diversity assay")
    }

    # parse sample names and q values from column names like 'Sample_q=0.01'
    sample_names <- sub("_q=.*", "", sample_q)
    has_q <- grepl("_q=", sample_q)
    if (!any(has_q)) {
        stop("Could not parse q values; expected '_q=' in column names", call. = FALSE)
    }
    if (!all(has_q)) {
        stop("Some column names are missing '_q='; ensure all diversity columns include a q value",
            call. = FALSE)
    }
    q_vals <- as.numeric(sub(".*_q=", "", sample_q))

    # determine group for each sample
    sample_type_in_coldata <- !is.null(sample_type_col) && sample_type_col %in% colnames(SummarizedExperiment::colData(se))
    if (sample_type_in_coldata) {
        st <- as.character(SummarizedExperiment::colData(se)[, sample_type_col])
        names(st) <- SummarizedExperiment::colData(se)$samples %||% colnames(mat)
        # when user provides a sample_type_col, index by the sample names
        # (strip q suffix)
        group_vec <- unname(st[sample_names])
    } else {
        stop("No sample grouping found: please supply `sample_type_col` or map sample",
            "types into `colData(se)` before calling calculate_lm_interaction().",
            call. = FALSE)
    }

    if (verbose && progress) {
        message("[calculate_lm_interaction] parsed samples and groups")
    }
    all_results <- list()
    fit_one <- function(g) {
        .tsenat_fit_one_interaction(g = g, se = se, mat = mat, q_vals = q_vals, sample_names = sample_names, group_vec = group_vec, method = method, pvalue = pvalue, subject_col = subject_col, paired = paired, min_obs = min_obs, verbose = verbose, suppress_lme4_warnings = suppress_lme4_warnings, progress = progress)
    }

    if (nthreads > 1 && .Platform$OS.type == "unix") {
        if (!requireNamespace("parallel", quietly = TRUE)) {
            stop("parallel package required for multi-threading")
        }
        res_list <- parallel::mclapply(rownames(mat), fit_one, mc.cores = nthreads)
    } else {
        res_list <- lapply(rownames(mat), fit_one)
    }
    all_results <- Filter(Negate(is.null), res_list)

    if (length(all_results) == 0) {
        return(data.frame())
    }
    res <- do.call(rbind, all_results)
    res$adj_p_interaction <- stats::p.adjust(res$p_interaction, method = "BH")
    # Sort first by adjusted p-values, then by raw p-values for stable ordering
    res <- res[order(res$adj_p_interaction, res$p_interaction), , drop = FALSE]
    rownames(res) <- NULL

    .tsenat_report_fit_summary(res, verbose = verbose)

    # Return the result data.frame (do not attach to or return a SummarizedExperiment)
    return(res)
}

# small helper (replacement for `%||%`) to provide default when NULL
`%||%` <- function(a, b) if (is.null(a)) b else a
