#' Calculate splicing diversity changes between two conditions.
#'
#' @param x A \code{matrix} with the splicing diversity values.
#' #' @param samples Character vector with an equal length to the number of columns
#'   in the input dataset, specifying the category of each sample.
#' @param control Name of the control sample category, defined in the
#'   \code{samples} vector, e.g. \code{control = 'Normal'} or \code{control =
#'   'WT'}.
#' @param method Method to use for calculating the average splicing diversity
#'   value in a condition. Can be \code{'mean'} or \code{'median'}.
#' @return A \code{data.frame} with mean or median value of splicing diversity
#' #'   across sample categories, the difference between these values and the log2
#'   fold change values.
#' #' @details The function uses a matrix of splicing diversity values in order to
#' calculate mean or median differences and log2 fold changes between two
#' conditions.
#' @import stats
calculate_fc <- function(x, samples, control, method = "mean") {
  if (method == "mean") {
    value <- aggregate(t(x), by = list(samples), mean, na.rm = TRUE)
  }

  if (method == "median") {
    value <- aggregate(t(x), by = list(samples), median, na.rm = TRUE)
  }

  sorted <- value[value$Group.1 != control, ]
  sorted[2, ] <- value[value$Group.1 == control, ]
  value <- t(sorted[, -1])
  value[is.na(value[, 1]), c(1)] <- NA
  value[is.na(value[, 2]), c(2)] <- NA

  result <- data.frame(value, ifelse(as.matrix(is.na(value[, 1]) | is.na(value[, 2])), NA, value[, 1] - value[
    ,
    2
  ]), ifelse(as.matrix(is.na(value[, 1]) | is.na(value[, 2])), NA, as.matrix(log(value[, 1] / value[, 2],
    base = 2
  ))))
  colnames(result) <- c(paste0(sorted[1, 1], "_", method), paste0(sorted[2, 1], "_", method), paste0(
    method,
    "_difference"
  ), "log2_fold_change")
  return(result)
}

#' Calculate p-values using Wilcoxon rank sum test.
#'
#' @param x A \code{matrix} with the splicing diversity values.
#' #' @param samples Character vector with an equal length to the number of columns
#' in the input dataset, specifying the category of each sample.
#' #' @param pcorr P-value correction method applied to the results, as defined in
#' the \code{p.adjust} function.
#' @param paired If \code{TRUE}, the Wilcox-test will be paired, and therefore
#' it will be a signed rank test instead of the rank sum test.
#' @param exact If \code{TRUE}, an exact p-value will be computed.
#' @return Raw and corrected p-values in a matrix.
#' @import stats
wilcoxon <- function(x, samples, pcorr = "BH", paired = FALSE, exact = FALSE) {
  p_values <- vector("list", nrow(x))

  for (i in seq_len(nrow(x))) {
    p_values[i] <- wilcox.test(x[i, as.numeric(which(samples %in% unique(sort(samples))[1]))], x[i, as.numeric(which(samples %in%
      unique(sort(samples))[2]))], paired = paired, exact = exact)$p.value
  }

  raw_p_values <- ifelse(is.na(vapply(p_values, c, numeric(1))), 1, vapply(p_values, c, numeric(1)))
  adjusted_p_values <- p.adjust(raw_p_values, method = pcorr)
  return(cbind(raw_p_values, adjusted_p_values))
}

#' Calculate p-values using label shuffling.
#'
#' @param x A \code{matrix} with the splicing diversity values.
#' #' @param samples Character vector with an equal length to the number of columns
#' in the input dataset, specifying the category of each sample.
#' @param control Name of the control sample category, defined in the
#' #' \code{samples} vector, e.g. \code{control = 'Normal'} or \code{control =
#' #' 'WT'}.
#' @param method Method to use for calculating the average splicing diversity
#' value in a condition. Can be \code{'mean'} or \code{'median'}.
#' @param randomizations The number of random shuffles.
#' #' @param pcorr P-value correction method applied to the results, as defined in
#' the \code{p.adjust} function.
#' @return Raw and corrected p-values.
#' @details The permutation p-values are computed two-sided as the proportion
#'  of permuted log2 fold-changes at least as extreme as the observed value,
#'  with a pseudocount added: (count + 1) / (n_perm + 1).
#' @import stats
#' @note The permutation test returns two-sided empirical p-values using a
#' pseudocount to avoid zero p-values for small numbers of permutations. See
#' the function documentation for details.
label_shuffling <- function(x, samples, control, method, randomizations = 100,
                            pcorr = "BH") {
  # observed log2 fold changes
  log2_fc <- calculate_fc(x, samples, control, method)[, 4]

  # build permutation/null distribution of log2 fold changes
  permuted <- replicate(randomizations, calculate_fc(x, sample(samples), control, method), simplify = FALSE)
  # each element is a data.frame/matrix; extract log2_fold_change column (4th column)
  perm_mat <- vapply(permuted, function(z) as.numeric(z[, 4]), numeric(nrow(x)))

  # compute two-sided permutation p-value with pseudocount: (count >= |obs| + 1) / (n_perm + 1)
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
    # if all permuted values are identical, fallback to p = 1 (no evidence)
    if (length(unique(nulls_non_na)) == 1) {
      return(1)
    }
    cnt <- sum(abs(nulls_non_na) >= abs(obs))
    pval <- (cnt + 1) / (n_non_na + 1)
    return(pval)
  }, numeric(1))

  adjusted_p_values <- p.adjust(raw_p_values, method = pcorr)
  return(cbind(raw_p_values, adjusted_p_values))
}
