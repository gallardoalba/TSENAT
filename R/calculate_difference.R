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
calculate_difference <- function(x, samples, control, method = "mean",
                                    test = "wilcoxon", randomizations = 100,
                                    pcorr = "BH", assayno = 1, verbose = FALSE,
                                    paired = FALSE, exact = FALSE, pseudocount = 0) {
    # internal small helpers (kept here to avoid adding new files)
    .tsenat_prepare_df <- function(x, samples, assayno) {
        if (inherits(
            x,
            "RangedSummarizedExperiment"
        ) || inherits(
            x,
            "SummarizedExperiment"
        )) {
            if (length(samples) != 1) {
                stop(
                    "'samples' must be a single colData column.",
                    call. = FALSE
                )
            }
            samples <- SummarizedExperiment::colData(x)[[samples]]
            if (!is.numeric(assayno) ||
                length(SummarizedExperiment::assays(x)) < assayno) {
                stop("Invalid 'assayno'.", call. = FALSE)
            }
            df <- as.data.frame(SummarizedExperiment::assays(x)[[assayno]])
            genes <- rownames(df)
            df <- cbind(genes = genes, df)
        } else {
            df <- as.data.frame(x)
        }
        list(df = df, samples = samples)
    }

    .tsenat_sample_matrix <- function(dfr) {
        as.matrix(dfr[,
            -c(
                1,
                ncol(dfr) - 1,
                ncol(dfr)
            ),
            drop = FALSE
        ])
    }

    # Validate input container
    # Reject matrices explicitly (tests expect this error for matrix input)
    if (is.matrix(x)) {
        stop(
            "Input type unsupported; see ?calculate_difference.",
            call. = FALSE
        )
    }
    if (!(is.data.frame(x) || inherits(
        x,
        "RangedSummarizedExperiment"
    ) || inherits(
        x,
        "SummarizedExperiment"
    ))) {
        stop(
            "Input data type not supported; see ?calculate_difference.",
            call. = FALSE
        )
    }

    # prepare data.frame and sample vector (handles SummarizedExperiment)
    pd <- .tsenat_prepare_df(x, samples, assayno)
    df <- pd$df
    samples <- pd$samples

    # Basic consistency checks
    if (ncol(df) - 1 != length(samples)) {
        stop("Column count doesn't match length(samples).", call. = FALSE)
    }
    groups <- levels(as.factor(samples))
    # capture dots early and initialize 'case' to avoid unbound variable errors
    # Note: any extra arguments are forwarded via `...` where used below.
    if (length(groups) > 2) {
        stop("More than two conditions; provide exactly two.", call. = FALSE)
    }
    if (length(groups) < 2) {
        stop("Fewer than two conditions; provide exactly two.", call. = FALSE)
    }
    if (!(control %in% samples)) {
        stop("Control sample type not found in samples.", call. = FALSE)
    }
    if (!(method %in% c(
        "mean",
        "median"
    ))) {
        stop("Invalid method; see ?calculate_difference.", call. = FALSE)
    }
    if (!(test %in% c(
        "wilcoxon",
        "shuffle"
    ))) {
        stop("Invalid test method; see ?calculate_difference.", call. = FALSE)
    }
    valid_pcorr <- c(
        "holm",
        "hochberg",
        "hommel",
        "bonferroni",
        "BH",
        "BY",
        "fdr",
        "none"
    )
    if (!(pcorr %in% valid_pcorr)) {
        stop(
            "Invalid p-value correction; see ?calculate_difference.",
            call. = FALSE
        )
    }

    # Informational warnings about sample size
    tab <- table(samples)
    if (test == "wilcoxon") {
        if (randomizations != 100 && verbose) {
            message("'randomizations' ignored for wilcoxon.")
        }
        if (any(tab < 3) || sum(tab) < 8) {
            warning("Low sample size for wilcoxon.", call. = FALSE)
        }
    }
    if (test == "shuffle") {
        if (sum(tab) <= 5) {
            warning("Low sample size for label shuffling.", call. = FALSE)
        }
        if (sum(tab) > 5 && sum(tab) < 10) {
            warning("Label shuffling may be unreliable.", call. = FALSE)
        }
    }

    # Count non-NA observations per gene per group
    idx1 <- which(samples == groups[1])
    idx2 <- which(samples == groups[2])
    df$cond_1 <- rowSums(!is.na(df[, idx1 + 1, drop = FALSE]))
    df$cond_2 <- rowSums(!is.na(df[, idx2 + 1, drop = FALSE]))

    # Partition genes by sufficient observations for the chosen test
    if (test == "wilcoxon") {
        keep_mask <- (
            df$cond_1 >= 3 &
                df$cond_2 >= 3 &
                (df$cond_1 + df$cond_2) >= 8
        )
    } else {
        keep_mask <- (df$cond_1 + df$cond_2) >= 5
    }

    df_keep <- df[keep_mask, , drop = FALSE]
    df_small <- df[!keep_mask, , drop = FALSE]

    result_list <- list()

    # Helper to extract the numeric matrix of sample columns (keeps original
    # order)
    sample_matrix <- .tsenat_sample_matrix

    if (nrow(df_keep) > 0) {
        if (nrow(df_small) > 0 && verbose) {
            message(sprintf(
                "Note: %d genes excluded due to low sample counts.",
                nrow(df_small)
            ))
        }
        ymat <- sample_matrix(df_keep)
        # p-value calculation
        if (test == "wilcoxon") {
            ptab <- wilcoxon(ymat, samples, pcorr = pcorr, paired = paired, exact = exact)
            # wilcoxon should return a data.frame of p-values named appropriately
        } else {
            ptab <- label_shuffling(
                ymat,
                samples,
                control,
                method,
                randomizations = randomizations,
                pcorr = pcorr,
                paired = paired
            )
        }
        result_list$tested <- data.frame(
            genes = df_keep[
                ,
                1
            ],
            calculate_fc(
                ymat,
                samples,
                control,
                method,
                pseudocount = pseudocount
            ),
            ptab,
            stringsAsFactors = FALSE
        )
    }

    if (nrow(df_small) > 0) {
        small_mat <- sample_matrix(df_small)
        result_list$small <- data.frame(
            genes = df_small[
                ,
                1
            ],
            calculate_fc(
                small_mat,
                samples,
                control,
                method,
                pseudocount = pseudocount
            ),
            raw_p_values = NA,
            adjusted_p_values = NA,
            stringsAsFactors = FALSE
        )
    }

    # Combine results preserving tested rows first
    if (length(result_list) == 0) {
        return(data.frame())
    }
    res <- do.call(rbind, result_list)
    rownames(res) <- NULL
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
#' "Normal").
#' @param min_obs Minimum number of non-NA observations required to fit a
#' model for a gene (default: 10).
#' @param method Modeling method to use for interaction testing: one of
#' \code{c("linear", "gam", "fpca")}.
#' @param nthreads Number of threads (mc.cores) to use for parallel processing
#' (default: 1).
#' @param assay_name Name of the assay in the SummarizedExperiment to use
#' (default: "diversity").
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
#'   sample_type = rep(c("Normal", "Tumor"), length.out = ncol(se)),
#'   row.names = colnames(se)
#' )
#' calculate_lm_interaction(se, sample_type_col = "sample_type")
calculate_lm_interaction <- function(se, sample_type_col = NULL, min_obs = 10,
                                        method = c(
                                            "linear",
                                            "lmm",
                                            "gam",
                                            "fpca"
                                        ),
                                        nthreads = 1,
                                        assay_name = "diversity") {
    method <- match.arg(method)
    message("[calculate_lm_interaction] method=", method)
    if (!requireNamespace("SummarizedExperiment",
        quietly = TRUE
    )) {
        stop("SummarizedExperiment required")
    }

    mat <- SummarizedExperiment::assay(se, assay_name)
    if (is.null(mat)) {
        stop(sprintf(
            "Assay '%s' not found in SummarizedExperiment",
            assay_name
        ))
    }

    sample_q <- colnames(mat)
    if (is.null(sample_q) || length(sample_q) == 0) {
        stop("No column names found on diversity assay")
    }

    # parse sample names and q values from column names like 'Sample_q=0.01'
    sample_names <- sub("_q=.*", "", sample_q)
    q_vals <- as.numeric(sub(".*_q=", "", sample_q))
    if (all(is.na(q_vals))) {
        stop("Could not parse q values; expected '_q=' in column names")
    }

    # determine group for each sample
    sample_type_in_coldata <- !is.null(sample_type_col) &&
        sample_type_col %in% colnames(SummarizedExperiment::colData(se))
    if (sample_type_in_coldata) {
        st <- as.character(SummarizedExperiment::colData(se)[, sample_type_col])
        names(st) <- SummarizedExperiment::colData(se)$samples %||% colnames(mat)
        # when user provides a sample_type_col, index by the sample names (strip q
        # suffix)
        group_vec <- unname(st[sample_names])
    } else {
        stop(
            "No sample grouping found: please supply `sample_type_col` or map sample", 
            "types into `colData(se)` before calling calculate_lm_interaction().",
            call. = FALSE
        )
    }

    message("[calculate_lm_interaction] parsed samples and groups")
    all_results <- list()
    fit_one <- function(g) {
        vals <- as.numeric(mat[g, ])
        df <- data.frame(entropy = vals, q = q_vals, group = factor(group_vec))
        if (sum(!is.na(df$entropy)) < min_obs) {
            return(NULL)
        }
        if (length(unique(na.omit(df$group))) < 2) {
            return(NULL)
        }

        if (method == "linear") {
            fit <- try(stats::lm(entropy ~ q * group, data = df), silent = TRUE)
            if (inherits(fit, "try-error")) {
                return(NULL)
            }
            coefs <- summary(fit)$coefficients
            ia_idx <- grep("^q:group", rownames(coefs))
            if (length(ia_idx) == 0) {
                return(NULL)
            }
            p_interaction <- coefs[ia_idx[1], "Pr(>|t|)"]
            return(data.frame(
                gene = g,
                p_interaction = p_interaction,
                stringsAsFactors = FALSE
            ))
        }

        if (method == "lmm") {
            if (!requireNamespace("lme4", quietly = TRUE)) {
                stop("Package 'lme4' is required for method = 'lmm'")
            }
            # subject identifier (base sample names without _q=...)
            subject <- sample_names
            df$subject <- factor(subject)
            # require at least two subjects and at least two groups represented
            if (length(unique(na.omit(df$subject))) < 2) return(NULL)
            if (length(unique(na.omit(df$group))) < 2) return(NULL)
            # fit null (no interaction) and alternative (with q:group interaction)
            fit0 <- try(lme4::lmer(entropy ~ q + group + (1 | subject), data = df, REML = FALSE), silent = TRUE)
            fit1 <- try(lme4::lmer(entropy ~ q * group + (1 | subject), data = df, REML = FALSE), silent = TRUE)
            if (inherits(fit0, "try-error") || inherits(fit1, "try-error")) return(NULL)
            an <- try(stats::anova(fit0, fit1), silent = TRUE)
            if (inherits(an, "try-error")) return(NULL)
            # p-value typically in column 'Pr(>Chisq)' in the second row
            p_interaction <- NA_real_
            if (nrow(an) >= 2) {
                pcol <- grep("Pr\\(>Chisq\\)|Pr\\(>F\\)|Pr\\(>Chi\\)", colnames(an), value = TRUE)
                if (length(pcol) == 0) {
                    # fallback to second column if available
                    p_interaction <- as.numeric(an[2, ncol(an)])
                } else {
                    p_interaction <- as.numeric(an[2, pcol[1]])
                }
            }
            return(data.frame(
                gene = g,
                p_interaction = p_interaction,
                stringsAsFactors = FALSE
            ))
        }

        if (method == "gam") {
            if (!requireNamespace("mgcv",
                quietly = TRUE
            )) {
                stop("Package 'mgcv' is required for method = 'gam'")
            }
            # choose smoothing basis dimension k based on number of unique q values
            uq_len <- length(unique(na.omit(q_vals)))
            k_q <- max(2, min(10, uq_len - 1))
            fit_null <- try(
                mgcv::gam(
                    entropy ~ group + s(q,
                        k = k_q
                    ),
                    data = df
                ),
                silent = TRUE
            )
            fit_alt <- try(
                mgcv::gam(
                    entropy ~ group + s(q,
                        by = group,
                        k = k_q
                    ),
                    data = df
                ),
                silent = TRUE
            )
            if (inherits(fit_null, "try-error") || inherits(fit_alt, "try-error")) {
                return(NULL)
            }
            an <- try(mgcv::anova.gam(fit_null, fit_alt, test = "F"), silent = TRUE)
            if (inherits(an, "try-error")) {
                return(NULL)
            }
            # anova.gam returns a table; p-value typically in second row, column
            # 'Pr(F)'
            p_interaction <- NA_real_
            if (nrow(an) >= 2) {
                if ("Pr(F)" %in% colnames(an)) {
                    p_interaction <- an[2, "Pr(F)"]
                } else if ("Pr(>F)" %in% colnames(an)) {
                    p_interaction <- an[2, "Pr(>F)"]
                } else if ("p-value" %in% colnames(an)) {
                    p_interaction <- an[
                        2,
                        "p-value"
                    ]
                }
            }
            return(data.frame(
                gene = g,
                p_interaction = p_interaction,
                stringsAsFactors = FALSE
            ))
        }

        if (method == "fpca") {
            message("[fpca] processing gene: ", g)
            # Build per-sample curves across q: rows = samples, cols = unique q values
            uq <- sort(unique(q_vals))
            samples_u <- unique(sample_names)
            message(sprintf(
                "[fpca] uq=%s samples_u=%s",
                paste(uq,
                    collapse = ","
                ),
                paste(samples_u,
                    collapse = ","
                )
            ))
            curve_mat <- matrix(NA_real_, nrow = length(samples_u), ncol = length(uq))
            # assign rownames defensively; avoid assigning colnames to prevent dimname
            # mismatch
            if (length(samples_u) > 0) rownames(curve_mat) <- samples_u
            message(sprintf(
                "[fpca] created curve_mat with dims: %s",
                paste(dim(curve_mat),
                    collapse = ","
                )
            ))
            for (i in seq_along(sample_names)) {
                s <- sample_names[i]
                qv <- q_vals[i]
                qi <- match(qv, uq)
                message("[fpca] i=", i, " s=", s, " qv=", qv, " qi=", qi)
                if (is.na(qi)) next
                # protect in case sample name not present in rownames
                if (s %in% rownames(curve_mat)) {
                    curve_mat[
                        s,
                        qi
                    ] <- as.numeric(mat[
                        g,
                        i
                    ])
                }
            }
            # keep samples with at least half of q points present
            good_rows <- which(rowSums(!is.na(curve_mat)) >= max(
                2,
                ceiling(ncol(curve_mat) / 2)
            ))
            if (length(good_rows) < min_obs) {
                return(NULL)
            }
            mat_sub <- curve_mat[good_rows, , drop = FALSE]
            # simple imputation for remaining NAs using column means
            col_means <- apply(mat_sub, 2, function(col) mean(col, na.rm = TRUE))
            for (r in seq_len(nrow(mat_sub))) {
                mat_sub[
                    r,
                    is.na(mat_sub[r, ])
                ] <- col_means[is.na(mat_sub[r, ])]
            }
            # perform PCA across q (observations = samples)
            pca <- try(
                stats::prcomp(mat_sub,
                    center = TRUE,
                    scale. = FALSE
                ),
                silent = TRUE
            )
            if (inherits(pca, "try-error")) {
                return(NULL)
            }
            if (ncol(pca$x) < 1) {
                return(NULL)
            }
            pc1 <- pca$x[, 1]
            # map groups to the rows used
            used_samples <- rownames(mat_sub)
            grp_vals <- group_vec[match(used_samples, sample_names)]
            # require two groups
            if (length(unique(na.omit(grp_vals))) < 2) {
                return(NULL)
            }
            # simple t-test between two groups
            g1 <- unique(na.omit(grp_vals))[1]
            g2 <- unique(na.omit(grp_vals))[2]
            x1 <- pc1[grp_vals == g1]
            x2 <- pc1[grp_vals == g2]
            if (length(x1) < 2 || length(x2) < 2) {
                return(NULL)
            }
            t_res <- try(stats::t.test(x1, x2), silent = TRUE)
            if (inherits(t_res, "try-error")) {
                return(NULL)
            }
            pval <- as.numeric(t_res$p.value)
            return(data.frame(
                gene = g,
                p_interaction = pval,
                stringsAsFactors = FALSE
            ))
        }
        return(NULL)
    }

    if (nthreads > 1 && .Platform$OS.type == "unix") {
        if (!requireNamespace("parallel",
            quietly = TRUE
        )) {
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
    return(res)
}

# small helper (replacement for `%||%`) to provide default when NULL
`%||%` <- function(a, b) if (is.null(a)) b else a
