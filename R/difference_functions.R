#' Calculate splicing diversity changes between two conditions.
#'
#' @param x A \code{matrix} with the splicing diversity values.
#' @param samples Character vector with an equal length to the number of columns
#' in the input dataset, specifying the category of each sample.
#' @param control Name of the control sample category, defined in the
#' \code{samples} vector, e.g. \code{control = 'Normal'} or \code{control =
#' 'WT'}.
#' @param method Method to use for calculating the average splicing diversity
#' value in a condition. Can be \code{'mean'} or \code{'median'}.
#' @param pseudocount Numeric scalar. Small value added to non-positive
#' observed group summaries to avoid zeros when computing differences and
#' log2 fold-changes. If \code{pseudocount <= 0} the function will automatically
#' choose a scale-aware value equal to half the smallest positive observed
#' group summary (i.e. half the smallest observed mean/median across groups);
#' if no positive values are present the fallback is \code{1e-6}. Rows with
#' insufficient observations remain \code{NA} and are not imputed.
#' @return A \code{data.frame} with mean or median value of splicing diversity
#' across sample categories, the difference between these values and the log2
#' fold change values.
#' @details The function uses a matrix of splicing diversity values in order to
#' calculate mean or median differences and log2 fold changes between two
#' conditions.
#' @import stats
calculate_fc <- function(x, samples, control, method = "mean", pseudocount = 0) {
    # validate control and samples inputs
    if (is.null(control) || !nzchar(control)) {
        stop("`control` must be provided to calculate_fc", call. = FALSE)
    }
    if (length(samples) != ncol(x)) {
        stop("Length of 'samples' must equal number of columns in 'x'", call. = FALSE)
    }
    if (!(control %in% samples)) {
        stop("Control sample type not found in samples.", call. = FALSE)
    }
    agg <- .tsenat_aggregate_fc_values(x = x, samples = samples, method = method,
        control = control)
    value <- agg$value
    sorted <- agg$sorted

    # Defensive numeric coercion: ensure group means are numeric and mark any
    # non-finite or non-positive values as NA. This prevents Inf/NaN when
    # computing log2 fold changes downstream.
    value <- matrix(as.numeric(value), nrow = nrow(value), ncol = ncol(value), dimnames = dimnames(value))
    value[!is.finite(value)] <- NA

    # compute and apply pseudocount based on observed group summaries
    value <- .tsenat_apply_pseudocount(value, pseudocount)

    # compute difference and log2 fold-change with NA-safe handling
    diff_vec <- value[, 1] - value[, 2]
    na_mask <- is.na(value[, 1]) | is.na(value[, 2])
    diff_vec[na_mask] <- NA

    log2fc_vec <- log2(value[, 1]/value[, 2])
    log2fc_vec[na_mask] <- NA

    result <- data.frame(value, difference = diff_vec, log2_fold_change = log2fc_vec,
        check.names = FALSE, stringsAsFactors = FALSE)
    colnames(result) <- c(paste(sorted[1, 1], "_", method, sep = ""), paste(sorted[2,
        1], "_", method, sep = ""), paste(method, "_difference", sep = ""), "log2_fold_change")
    return(result)
}

#' Calculate p-values using Wilcoxon rank sum test.
#'
#' @param x A \code{matrix} with the splicing diversity values.
#' @param samples Character vector with an equal length to the number of columns
#' in the input dataset, specifying the category of each sample.
#' @param pcorr P-value correction method applied to the results, as defined in
#' the \code{p.adjust} function.
#' @param paired If \code{TRUE}, the Wilcox-test will be paired, and therefore
#' it will be a signed rank test instead of the rank sum test.
#' @param exact If \code{TRUE}, an exact p-value will be computed.
#' @return Raw and corrected p-values in a matrix.
#' @import stats
wilcoxon <- function(x, samples, pcorr = "BH", paired = FALSE, exact = FALSE) {
    # Determine group indices (two groups expected)
    groups <- unique(sort(samples))
    if (length(groups) != 2) {
        stop("`samples` must contain exactly two groups for Wilcoxon tests.")
    }
    g1_idx <- as.numeric(which(samples %in% groups[1]))
    g2_idx <- as.numeric(which(samples %in% groups[2]))

    if (isTRUE(paired)) {
        if (length(g1_idx) != length(g2_idx)) {
            stop("Paired Wilcoxon requires equal numbers of samples in each ", "group.")
        }
        # Paired tests assume columns are already ordered/aligned by the caller
        # (e.g., via `map_metadata()`); do not attempt to infer pairing here.
    }

    p_values <- vector("list", nrow(x))
    for (i in seq_len(nrow(x))) {
        p_values[i] <- tryCatch({
            wilcox.test(x[i, g1_idx], x[i, g2_idx], paired = paired, exact = exact)$p.value
        }, error = function(e) {
            NA_real_
        }, warning = function(w) {
            # swallow specific warnings but return NA on unusual states
            NA_real_
        })
    }

    raw_p_values <- ifelse(is.na(vapply(p_values, c, numeric(1))), 1, vapply(p_values,
        c, numeric(1)))
    adjusted_p_values <- p.adjust(raw_p_values, method = pcorr)
    out <- cbind(raw_p_values, adjusted_p_values)
    colnames(out) <- c("raw_p_values", "adjusted_p_values")
    return(out)
}

#' Calculate p-values using label shuffling.
#'
#' @param x A \code{matrix} with the splicing diversity values.
#' @param samples Character vector with an equal length to the number of columns
#' in the input dataset, specifying the category of each sample.
#' @param control Name of the control sample category, defined in the
#' \code{samples} vector, e.g. \code{control = 'Normal'} or \code{control =
#' 'WT'}.
#' @param method Method to use for calculating the average splicing diversity
#' value in a condition. Can be \code{'mean'} or \code{'median'}.
#' @param randomizations The number of random shuffles.
#' @param pcorr P-value correction method applied to the results, as defined in
#' the \code{p.adjust()} function.
#' @param paired Logical; if \code{TRUE} perform a paired permutation scheme
#'   (default: \code{FALSE}). When paired is \code{TRUE}, permutations
#'   should preserve pairing between samples; the function currently permutes
#'   sample labels and therefore paired analyses are only meaningful when the
#'   caller has arranged \code{samples} accordingly.
#' @param paired_method Character; method for paired permutations. One of
#'   \code{'swap'} (randomly swap labels within pairs) or \code{'signflip'}
#'   (perform sign-flip permutations; can enumerate all 2^n_pairs combinations
#'   for an exact test when \code{randomizations = 0} or \code{randomizations >= 2^n_pairs}).
#' @return Raw and corrected p-values.
#' @details The permutation p-values are computed two-sided as the proportion
#' of permuted log2 fold-changes at least as extreme as the observed value,
#' with a pseudocount added: (count + 1) / (n_perm + 1).
#' @import stats
#' @note The permutation test returns two-sided empirical p-values using a
#' pseudocount to avoid zero p-values for small numbers of permutations. See
#' the function documentation for details.
label_shuffling <- function(x, samples, control, method, randomizations = 100, pcorr = "BH",
    paired = FALSE, paired_method = c("swap", "signflip")) {
    paired_method <- match.arg(paired_method)
    # observed log2 fold changes
    log2_fc <- calculate_fc(x, samples, control, method)[, 4]

    # build permutation/null distribution of log2 fold changes
    if (isTRUE(paired)) {
        perm_mat <- .tsenat_permute_paired(x = x, samples = samples, control = control,
            method = method, randomizations = randomizations, paired_method = paired_method)
    } else {
        permuted <- replicate(randomizations, calculate_fc(x, sample(samples), control,
            method), simplify = FALSE)
        # each element is a data.frame/matrix; extract log2_fold_change column
        # (4th column)
        perm_mat <- vapply(permuted, function(z) as.numeric(z[, 4]), numeric(nrow(x)))
    }

    # compute two-sided permutation p-value with pseudocount: (count >= |obs| +
    # 1) / (n_perm + 1)
    raw_p_values <- vapply(seq_len(nrow(perm_mat)), function(i) {
        obs <- log2_fc[i]
        nulls <- perm_mat[i, ]
        if (is.na(obs) || all(is.na(nulls))) {
            return(1)
        }
        nulls_non_na <- nulls[!is.na(nulls)]
        n_non_na <- length(nulls_non_na)
        if (n_non_na == 0) {
            return(1)
        }
        cnt <- sum(abs(nulls_non_na) >= abs(obs))
        pval <- (cnt + 1)/(n_non_na + 1)
        return(pval)
    }, numeric(1))

    adjusted_p_values <- p.adjust(raw_p_values, method = pcorr)
    out <- cbind(raw_p_values, adjusted_p_values)
    colnames(out) <- c("raw_p_values", "adjusted_p_values")
    return(out)
}


#' Run a differential test by name: Wilcoxon or label-shuffle
#'
#' Thin wrapper that selects between `wilcoxon()` and `label_shuffling()`.
#'
#' @param x A matrix with splicing diversity values (rows = features).
#' @param samples Character vector of sample group labels (length = ncol(x)).
#' @param control Name of the control group (required for label-shuffle).
#' @param method Character; one of `'wilcoxon'` or `'shuffle'`.
#' @param fc_method Character; aggregation method used by the permutation
#'   test when `method = 'shuffle'` ('mean' or 'median').
#' @param paired Logical passed to `wilcoxon()` when using the Wilcoxon test.
#' @param exact Logical passed to `wilcoxon()` to request exact p-values.
#' @param randomizations Integer number of permutations for `label_shuffling()`.
#' @param pcorr P-value adjustment method (passed to `p.adjust`).
#' @param seed Integer seed used to make permutations reproducible (default
#'   123). The function sets a temporary RNG seed via `withr::local_seed(seed)`
#'   before running `label_shuffling()` when `method = 'shuffle'`.
#' @return A two-column matrix with raw and adjusted p-values (as returned by
#'   the underlying functions).
#' @export
#' @examples
#' mat <- matrix(rnorm(20), nrow = 5)
#' samples <- rep(c('A', 'B'), length.out = ncol(mat))
#' test_differential(mat, samples, control = 'A', method = 'wilcoxon')
#'
#' @param paired_method Character; forwarded to `label_shuffling()` when
#'   `method = 'shuffle'`. See `label_shuffling()` for details.
test_differential <- function(x, samples, control = NULL, method = c("wilcoxon",
    "shuffle"), fc_method = "mean", paired = FALSE, exact = FALSE, randomizations = 100,
    pcorr = "BH", seed = 123L, paired_method = c("swap", "signflip")) {
    paired_method <- match.arg(paired_method)
    method <- match.arg(method)
    if (method == "wilcoxon") {
        return(wilcoxon(x, samples, pcorr = pcorr, paired = paired, exact = exact))
    }

    # shuffle / permutation-based test
    if (is.null(control) || !nzchar(control)) {
        stop("`control` must be provided when method = 'shuffle'", call. = FALSE)
    }
    if (!is.null(seed)) {
        if (!is.numeric(seed) || length(seed) != 1) {
            stop("`seed` must be a single numeric value", call. = FALSE)
        }
        # Use withr::local_seed to temporarily set RNG state for reproducible
        # permutation tests without permanently modifying the global RNG.
        withr::local_seed(as.integer(seed))
    }
    return(label_shuffling(x, samples, control, fc_method, randomizations = randomizations,
        pcorr = pcorr, paired = paired, paired_method = paired_method))
}
