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
#'   should preserve pairing between samples.
#' @param paired_method Character; method for paired permutations. One of
#'   \code{'swap'} (randomly swap labels within pairs) or \code{'signflip'}
#'   (perform sign-flip permutations; can enumerate all 2^n_pairs combinations
#'   for  an exact test when  \code{randomizations = 0} or 
#' \code{randomizations >= 2^n_pairs}).
#' @param pairs Optional character vector with an equal length to the number of 
#' columns in the input dataset, specifying the pairing identifier for each
#' sample.
#' When provided with \code{paired = TRUE}, samples are matched based on this 
#' pairing information.
#' @param nthreads Number of threads for parallel processing (default: 1).
#'   Set to > 1 to parallelize per-feature p-value computation.
#' @param robust_loss_type Character; loss function for M-estimation (Tukey,
#' Huber, or other).
#'   Used when non-parametric tests switch to robust parametric alternatives.
#'   Default: 'huber'.
#' @param robust_scale_method Character; scale selection method for M-estimation
#'   (e.g., 'mad' for median absolute deviation). Default: 'mad'.
#' @return Raw and corrected p-values.
#' @details
#' \strong{S019 Implementation: Phipson & Smyth (2010) Bias Correction}
#'
#' This function implements the critical p-value correction from Phipson &
#' Smyth (2010):
#' \deqn{p = \frac{b + 1}{m + 1}}{p = (b + 1) / (m + 1)}
#'
#' Instead of the traditional formula p = b/m,
#'  where \code{b} is the count of permutations
#' with  |test_statistic| >= |observed_statistic| and 
#' \code{m} is the total number of
#' permutations.
#'
#' \strong{Why This Correction Matters:}
#' - \strong{Prevents p = 0:} Traditional formula produces p = 0 when observed
#' statistic is more extreme than all m permutations. This is statistically
#' incorrect.
#' - \strong{Proper Calibration:} The pseudocount ensures valid Type I error
#' control
#' and proper coverage properties, especially important with small
#' permutation counts.
#' - \strong{Minimum P-Value:} With m permutations, p_min = 1/(m+1), not 0.
#'   Example: With m = 1000, p_min ~= 0.000999 (not 0).
#' - \strong{Standard Practice:} This correction is now implemented in
#' limma, edgeR,
#'   DESeq2, and other standard bioinformatics packages.
#'
#' The permutation p-values are computed two-sided as the proportion
#' of permuted log2 fold-changes at least as extreme as the observed value,
#' with the pseudocount applied: (count + 1) / (n_perm + 1).
#' 
#' For paired designs, the function supports two permutation schemes:
#' - \code{'swap'}: Randomly swaps sample labels within pairs
#' - \code{'signflip'}: Performs sign-flip permutations (Pesarin & Salmaso 2010)
#' @note The permutation test returns two-sided empirical p-values using the
#' Phipson & Smyth (2010) pseudocount correction to avoid zero p-values.
#' This ensures proper statistical calibration regardless of the number of
#' permutations.
#' @references
#' Phipson, B., and Smyth, G. K. (2010). Permutation p-values should never
#' be zero:
#' calculating exact p-values when permutations are randomly drawn.
#' Statistical Applications in Genetics and Molecular Biology, 9(1), 39.
#' DOI: 10.2202/1544-6115.1585
#'
#' Pesarin, F., and Salmaso, L. (2010). Permutation Tests for Complex Data:
#' Theory, Applications and Software. John Wiley & Sons.
#'
#' Good, P. I. (2005). Permutation, Parametric and Bootstrap Tests of Hypotheses
#' (3rd ed.). Springer Series in Statistics.
#' @noRd
#' @examples
#' set.seed(123)
#' # Create a matrix of splicing diversity values (2 genes x 4 samples)
#' mat <- matrix(rnorm(8), nrow = 2)
#' samples <- c('Normal', 'Normal', 'Tumor', 'Tumor')
#' 
#' # Run label shuffling test with S019 correction (100 permutations)
#' # P-values will follow (b+1)/(m+1) formula with m=100
#' result <- .label_shuffling(mat, samples, control = 'Normal', 
#'                           method = 'mean',  randomizations = 100,
#'  pcorr = 'BH')
#' head(result)
.label_shuffling <- function(x, samples, control, method, randomizations = 100, pcorr = "BH",
    paired = FALSE, paired_method = c("swap", "signflip"), nthreads = 1, pairs = NULL,
    robust_loss_type = "huber", robust_scale_method = "mad") {
    paired_method <- match.arg(paired_method)

    # Validate inputs and get observed statistics
    .validate_label_shuffling_inputs(samples, control, pairs, ncol(x), paired)
    fc_result <- .calculate_fc(x, samples, control, method)
    log2_fc <- fc_result[, 4]
    group_means <- fc_result[, seq_len(2)]

    # Pre-compute pseudocount and group identifiers (computed once, reused)
    pseudocount_val <- .compute_pseudocount(fc_result)
    groups <- sort(unique(samples))  # OPTIMIZATION: Compute once, use for both perm and effect sizes
    case_group <- setdiff(groups, control)

    # Generate permutation null distribution
    perm_mat <- .generate_permutation_matrix(x, samples, control, method, 
        randomizations, paired, paired_method, pairs, pseudocount_val, 
        robust_loss_type, robust_scale_method)

    # Compute p-values from permutation distribution (with S019 correction)
    raw_p_values <- .compute_pvalues_from_permutations(log2_fc, perm_mat, nthreads)
    adjusted_p_values <- p.adjust(raw_p_values, method = pcorr)

    # Compute effect sizes (independent of permutation distribution)
    effect_stats <- .compute_all_effect_sizes(x, samples, pairs, groups, nthreads)

    # Format and return results
    .format_pvalue_output(raw_p_values, adjusted_p_values, log2_fc, 
        effect_stats, group_means)
}

#' Validate inputs for .label_shuffling()
#'
#' @param samples Character vector of sample group labels
#' @param control Control group name
#' @param pairs Optional pairing vector
#' @param n_samples Number of samples (ncol of data matrix)
#' @param paired Logical; whether paired design
#'
#' @noRd
.validate_label_shuffling_inputs <- function(samples, control, pairs, n_samples, paired) {
    if (!(control %in% samples)) {
        stop("Control group '", control, "' not found in unique sample types: ",
            paste(unique(samples), collapse = ", "), call. = FALSE)
    }

    unique_groups <- unique(samples)
    if (length(unique_groups) != 2) {
        stop(".label_shuffling() requires exactly 2 sample groups; found ",
            length(unique_groups), ": ", paste(unique_groups, collapse = ", "), 
            call. = FALSE)
    }

    # Check for even number of samples when paired design
    if (isTRUE(paired) && (n_samples %% 2 != 0)) {
        stop("Paired permutation requires an even number of samples",
            call. = FALSE)
    }

    # Check pairs parameter is provided when paired design
    if (isTRUE(paired) && is.null(pairs)) {
        stop("paired=TRUE requires `pairs` parameter to be provided.", call. = FALSE)
    }

    if (isTRUE(paired) && !is.null(pairs) && length(pairs) != n_samples) {
        stop("`pairs` must have length equal to n_samples.", call. = FALSE)
    }
}

#' Compute pseudocount from fold-change results
#'
#' @param fc_result Data frame from .calculate_fc()
#'
#' @noRd
.compute_pseudocount <- function(fc_result) {
    pos_vals <- as.matrix(fc_result[, seq_len(2)])
    pos_vals <- pos_vals[!is.na(pos_vals) & pos_vals > 0]
    if (length(pos_vals) > 0) min(pos_vals, na.rm = TRUE) / 2 else 1e-06
}

#' Generate permutation matrix for label shuffling test
#'
#' @noRd
.generate_permutation_matrix <- function(x, samples, control, method, randomizations, 
    paired, paired_method, pairs, pseudocount_val, robust_loss_type, robust_scale_method) {
    
    perm_mat <- matrix(NA_real_, nrow = nrow(x), ncol = randomizations)
    unique_groups <- unique(samples)
    case_group <- setdiff(unique_groups, control)

    if (isTRUE(paired) && !is.null(pairs)) {
        perm_mat <- .generate_paired_permutations(x, samples, control, case_group, 
            method, randomizations, pairs, paired_method, pseudocount_val, 
            robust_loss_type, robust_scale_method)
    } else {
        perm_mat <- .generate_unpaired_permutations(x, samples, control, case_group, 
            method, randomizations, pseudocount_val, robust_loss_type, robust_scale_method)
    }

    perm_mat
}

#' Generate paired permutations with pair structure preservation
#'
#' @noRd
.generate_paired_permutations <- function(x, samples, control, case_group, method, 
    randomizations, pairs, paired_method, pseudocount_val, robust_loss_type, robust_scale_method) {
    
    perm_mat <- matrix(NA_real_, nrow = nrow(x), ncol = randomizations)
    unique_pairs <- unique(pairs)
    pair_indices <- .prepare_pair_indices(pairs, unique_pairs)
    n_pairs <- length(unique_pairs)
    
    # OPTIMIZATION: Pre-generate all random decisions for all randomizations at once
    # This avoids repeated sample() calls inside the loop
    if (paired_method == "swap") {
        # For swap: generate (randomizations × n_pairs) boolean matrix
        swap_matrix <- matrix(sample(c(TRUE, FALSE), size = randomizations * n_pairs, replace = TRUE), 
                              nrow = randomizations, ncol = n_pairs)
    } else if (paired_method == "signflip") {
        # For signflip: generate (randomizations × n_pairs) boolean matrix
        swap_matrix <- matrix(sample(c(TRUE, FALSE), size = randomizations * n_pairs, replace = TRUE), 
                              nrow = randomizations, ncol = n_pairs)
    }

    for (r in seq_len(randomizations)) {
        perm_samples <- samples
        
        # Apply pre-computed permutations for this randomization
        for (p_idx in seq_len(n_pairs)) {
            pair_idx <- pair_indices[[p_idx]]
            if (length(pair_idx) == 2 && swap_matrix[r, p_idx]) {
                perm_samples[pair_idx] <- perm_samples[rev(pair_idx)]
            }
        }

        perm_case_idx <- which(perm_samples == case_group)
        perm_ctrl_idx <- which(perm_samples == control)
        perm_mat[, r] <- .fast_log2fc_permutation(x, perm_case_idx, perm_ctrl_idx,
            method, pseudocount_val, robust_loss_type, robust_scale_method)
    }

    perm_mat
}

#' Prepare pair indices for vectorized pair operations
#'
#' @noRd
.prepare_pair_indices <- function(pairs, unique_pairs) {
    # Vectorized version using split() to avoid repeated which() calls
    pair_split <- split(seq_along(pairs), pairs)
    pair_split[match(unique_pairs, names(pair_split))]
}

#' Generate unpaired permutations
#'
#' @noRd
.generate_unpaired_permutations <- function(x, samples, control, case_group, method, 
    randomizations, pseudocount_val, robust_loss_type, robust_scale_method) {
    
    perm_mat <- matrix(NA_real_, nrow = nrow(x), ncol = randomizations)
    
    # Note: permutation tests use current RNG state set by user.

    for (r in seq_len(randomizations)) {
        perm_samples <- sample(samples)
        perm_case_idx <- which(perm_samples == case_group)
        perm_ctrl_idx <- which(perm_samples == control)
        perm_mat[, r] <- .fast_log2fc_permutation(x, perm_case_idx, perm_ctrl_idx,
            method, pseudocount_val, robust_loss_type, robust_scale_method)
    }

    perm_mat
}

#' Compute p-values from permutation distribution (S019 Phipson & Smyth correction)
#'
#' @noRd
.compute_pvalues_from_permutations <- function(log2_fc, perm_mat, nthreads) {
    # OPTIMIZATION: Pre-compute row counts of non-NA values to avoid repeated filtering
    # Each element [i] = count of non-NA permutation values for feature i
    n_non_na_perm <- rowSums(!is.na(perm_mat))
    
    .compute_pval <- function(i) {
        obs <- log2_fc[i]
        if (is.na(obs)) return(1)
        
        nulls <- perm_mat[i, ]
        n_non_na <- n_non_na_perm[i]
        if (n_non_na == 0) return(1)
        
        # Avoid creating intermediate nulls_non_na vector for memory efficiency
        cnt <- sum(abs(nulls[!is.na(nulls)]) >= abs(obs))
        (cnt + 1) / (n_non_na + 1)  # S019: Phipson & Smyth (2010)
    }

    unlist(.bplapply(seq_len(nrow(perm_mat)), .compute_pval, nthreads = nthreads))
}


#' Compute effect sizes (U and r) for all features
#'
#' @noRd
.compute_all_effect_sizes <- function(x, samples, pairs, groups, nthreads) {
    effect_sizes <- .bplapply(seq_len(nrow(x)), function(i) {
        .compute_one_effect_size(x, i, samples, pairs, groups)
    }, nthreads = nthreads)

    # Simplified extraction: convert list of named vectors to separate vectors
    list(
        U = vapply(effect_sizes, "[", i = "U", FUN.VALUE = numeric(1)),
        r = vapply(effect_sizes, "[", i = "r", FUN.VALUE = numeric(1))
    )
}

#' Compute Wilcoxon U and effect size r for a single feature
#'
#' @noRd
.compute_one_effect_size <- function(x, feature_idx, samples, pairs, groups) {
    tryCatch({
        if (!is.null(pairs)) {
            .wilcox_effect_sizes_paired(x, feature_idx, pairs, samples, groups)
        } else {
            .wilcox_effect_sizes_unpaired(x, feature_idx, samples, groups)
        }
    }, error = function(e) {
        c(U = NA_real_, r = NA_real_)
    })
}

#' Wilcoxon effect sizes for paired design
#'
#' @noRd
.wilcox_effect_sizes_paired <- function(x, feature_idx, pairs, samples, groups) {
    unique_pairs <- unique(pairs)
    pair_indices <- .prepare_pair_indices(pairs, unique_pairs)
    
    # Pre-allocate difference vector (more efficient than growing with c())
    diffs_list <- vector("list", length(unique_pairs))
    
    for (p_idx in seq_along(unique_pairs)) {
        pair_idx <- pair_indices[[p_idx]]
        # Use pre-computed indices instead of which() for each pair
        g1_mask <- samples[pair_idx] == groups[1]
        g2_mask <- samples[pair_idx] == groups[2]
        
        # Extract indices for this pair's groups
        g1_idx <- pair_idx[g1_mask]
        g2_idx <- pair_idx[g2_mask]
        
        if (length(g1_idx) == 1 && length(g2_idx) == 1) {
            diffs_list[[p_idx]] <- x[feature_idx, g1_idx] - x[feature_idx, g2_idx]
        }
    }
    
    all_diffs <- unlist(diffs_list)

    if (length(all_diffs) < 2) return(c(U = NA_real_, r = NA_real_))

    wt <- wilcox.test(all_diffs, mu = 0, exact = FALSE)
    U <- as.numeric(wt$statistic)
    n <- length(all_diffs)
    
    # r = Z / sqrt(n) for paired test
    expected_U <- n * (n + 1) / 4
    var_U <- (n * (n + 1) * (2 * n + 1)) / 24
    
    # Guard against zero variance (can occur with very small n)
    if (var_U <= 0) return(c(U = U, r = NA_real_))
    
    Z <- (U - expected_U) / sqrt(var_U)
    r <- Z / sqrt(n)
    
    c(U = U, r = pmax(-1, pmin(1, r)))  # Clamp r to [-1, 1]
}

#' Wilcoxon effect sizes for unpaired design
#'
#' @noRd
.wilcox_effect_sizes_unpaired <- function(x, feature_idx, samples, groups) {
    g1_idx <- which(samples == groups[1])
    g2_idx <- which(samples == groups[2])

    if (length(g1_idx) == 0 || length(g2_idx) == 0) {
        return(c(U = NA_real_, r = NA_real_))
    }

    wt <- wilcox.test(x[feature_idx, g1_idx], x[feature_idx, g2_idx], 
        paired = FALSE, exact = FALSE)
    U <- as.numeric(wt$statistic)
    n1 <- length(g1_idx)
    n2 <- length(g2_idx)
    n <- n1 + n2
    
    # r = Z / sqrt(n) for unpaired test
    expected_U <- n1 * n2 / 2
    var_U <- (n1 * n2 * (n1 + n2 + 1)) / 12
    
    # Guard against zero variance
    if (var_U <= 0) return(c(U = U, r = NA_real_))
    
    Z <- (U - expected_U) / sqrt(var_U)
    r <- Z / sqrt(n)
    
    c(U = U, r = pmax(-1, pmin(1, r)))  # Clamp r to [-1, 1]
}

#' Format output as data frame with all statistics
#'
#' @noRd
.format_pvalue_output <- function(raw_p_values, adjusted_p_values, log2_fc, 
    effect_stats, group_means) {
    
    out <- data.frame(
        pvalue = raw_p_values,
        padj = adjusted_p_values,
        log2FC = log2_fc,
        U = effect_stats$U,
        r = effect_stats$r,
        group_means,
        check.names = FALSE,
        stringsAsFactors = FALSE
    )

    # Extract group names by removing method suffix (e.g., "Normal_mean" -> "Normal")
    group_col_names <- colnames(group_means)
    group_names <- gsub("_(mean|median|m_estimate)$", "", group_col_names)
    colnames(out) <- c("pvalue", "padj", "log2FC", "U", "r", group_names)
    
    out
}
