# Helper utilities for calculate_difference

.tsenat_calculate_difference_partition <- function(df, samples, control, method, test, pcorr, randomizations, verbose) {
    if (ncol(df) - 1 != length(samples)) stop("Column count doesn't match length(samples).", call. = FALSE)
    uniq_groups <- unique(as.character(samples))
    if (length(uniq_groups) > 2) stop("More than two conditions; provide exactly two.", call. = FALSE)
    if (length(uniq_groups) < 2) stop("Fewer than two conditions; provide exactly two.", call. = FALSE)
    if (!(control %in% uniq_groups)) stop("Control sample type not found in samples.", call. = FALSE)

    case_label <- setdiff(uniq_groups, control)
    groups <- c(case_label, control)

    if (!(method %in% c("mean", "median"))) stop("Invalid method; see ?calculate_difference.", call. = FALSE)
    if (!(test %in% c("wilcoxon", "shuffle"))) stop("Invalid test method; see ?calculate_difference.", call. = FALSE)
    valid_pcorr <- c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
    if (!(pcorr %in% valid_pcorr)) stop("Invalid p-value correction; see ?calculate_difference.", call. = FALSE)

    tab <- table(samples)
    if (test == "wilcoxon") {
        if (randomizations != 100 && verbose) message("'randomizations' ignored for wilcoxon.")
        if (any(tab < 3) || sum(tab) < 8) warning("Low sample size for wilcoxon.", call. = FALSE)
    }
    if (test == "shuffle") {
        if (sum(tab) <= 5) warning("Low sample size for label shuffling.", call. = FALSE)
        if (sum(tab) > 5 && sum(tab) < 10) warning("Label shuffling may be unreliable.", call. = FALSE)
    }

    idx_case <- which(samples == groups[1])
    idx_control <- which(samples == groups[2])
    idx1 <- idx_case
    idx2 <- idx_control

    df$cond_1 <- rowSums(!is.na(df[, idx1 + 1, drop = FALSE]))
    df$cond_2 <- rowSums(!is.na(df[, idx2 + 1, drop = FALSE]))

    if (test == "wilcoxon") {
        keep_mask <- (df$cond_1 >= 3 & df$cond_2 >= 3 & (df$cond_1 + df$cond_2) >= 8)
    } else {
        keep_mask <- (df$cond_1 + df$cond_2) >= 5
    }

    df_keep <- df[keep_mask, , drop = FALSE]
    df_small <- df[!keep_mask, , drop = FALSE]

    list(df = df, samples = samples, groups = groups, idx1 = idx1, idx2 = idx2, df_keep = df_keep, df_small = df_small)
}

# Helpers for calculate_fc
.tsenat_aggregate_fc_values <- function(x, samples, method, control) {
    if (method == "mean") {
        value <- aggregate(t(x), by = list(samples), mean, na.rm = TRUE)
    } else if (method == "median") {
        value <- aggregate(t(x), by = list(samples), median, na.rm = TRUE)
    } else {
        stop("Invalid method; must be 'mean' or 'median'")
    }

    sorted <- value[value$Group.1 != control, ]
    sorted[2, ] <- value[value$Group.1 == control, ]
    value <- t(sorted[, -1])
    value[is.na(value[, 1]), c(1)] <- NA
    value[is.na(value[, 2]), c(2)] <- NA
    return(list(value = value, sorted = sorted))
}

.tsenat_apply_pseudocount <- function(value, pseudocount) {
    if (!is.numeric(pseudocount) || length(pseudocount) != 1) {
        pseudocount <- 0
    }
    if (pseudocount <= 0) {
        pos_vals <- value[!is.na(value) & value > 0]
        if (length(pos_vals) > 0) {
            pc <- min(pos_vals, na.rm = TRUE) / 2
        } else {
            pc <- 1e-06
        }
    } else {
        pc <- pseudocount
    }
    replace_idx <- !is.na(value) & value <= 0
    value[replace_idx] <- pc
    return(value)
}

# Paired permutation helpers
.tsenat_permute_paired <- function(x, samples, control, method, randomizations, paired_method) {
    ncols <- ncol(x)
    if (ncols %% 2 != 0) stop("Paired permutation requires an even number of samples and paired column ordering", call. = FALSE)
    npairs <- ncols / 2
    if (paired_method == "swap") {
        perm_mat <- matrix(NA_real_, nrow = nrow(x), ncol = randomizations)
        for (r in seq_len(randomizations)) {
            swap <- sample(c(TRUE, FALSE), size = npairs, replace = TRUE)
            perm_samples <- samples
            for (p in seq_len(npairs)) {
                if (swap[p]) {
                    i1 <- (p - 1) * 2 + 1
                    i2 <- i1 + 1
                    perm_samples[c(i1, i2)] <- perm_samples[c(i2, i1)]
                }
            }
            df_perm <- calculate_fc(x, perm_samples, control, method)
            perm_mat[, r] <- as.numeric(df_perm[, 4])
        }
        return(perm_mat)
    } else if (paired_method == "signflip") {
        total_comb <- 2^npairs
        if (randomizations <= 0 || randomizations >= total_comb) {
            combos <- expand.grid(rep(list(c(0, 1)), npairs))
            nrep <- nrow(combos)
            perm_mat <- matrix(NA_real_, nrow = nrow(x), ncol = nrep)
            for (r in seq_len(nrep)) {
                swap <- as.logical(as.integer(combos[r, ]))
                perm_samples <- samples
                for (p in seq_len(npairs)) {
                    if (swap[p]) {
                        i1 <- (p - 1) * 2 + 1
                        i2 <- i1 + 1
                        perm_samples[c(i1, i2)] <- perm_samples[c(i2, i1)]
                    }
                }
                df_perm <- calculate_fc(x, perm_samples, control, method)
                perm_mat[, r] <- as.numeric(df_perm[, 4])
            }
            return(perm_mat)
        } else {
            perm_mat <- matrix(NA_real_, nrow = nrow(x), ncol = randomizations)
            for (r in seq_len(randomizations)) {
                swap <- sample(c(TRUE, FALSE), size = npairs, replace = TRUE)
                perm_samples <- samples
                for (p in seq_len(npairs)) {
                    if (swap[p]) {
                        i1 <- (p - 1) * 2 + 1
                        i2 <- i1 + 1
                        perm_samples[c(i1, i2)] <- perm_samples[c(i2, i1)]
                    }
                }
                df_perm <- calculate_fc(x, perm_samples, control, method)
                perm_mat[, r] <- as.numeric(df_perm[, 4])
            }
            return(perm_mat)
        }
    }
}
