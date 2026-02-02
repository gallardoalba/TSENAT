#' Calculate splicing diversity changes between two conditions.
#'
#' @param x A \code{SummarizedExperiment} with splicing diversity values for
#'   each gene in each sample or a \code{data.frame} with gene names in the
#'   first column and splicing diversity values for each sample in additional
#'   columns.
#' @param samples A vector of length one, specifying the column name of the
#'   \code{colData} annotation column from the \code{SummarizedExperiment}
#'   object, that should be used as the category column or a character vector
#'   with an equal length to the number of columns in the input dataset,
#'   specifying the category of each sample in the case of a \code{data.frame}
#'   input.
#' @param control Name of the control sample category, defined in the
#'   \code{samples} vector, e.g. \code{control = 'Normal'} or \code{control =
#'   'WT'}.
#' @param method Method to use for calculating the average splicing diversity
#'   value in a condition. Can be \code{'mean'} or \code{'median'}.
#' #' @param test Method to use for p-value calculation: use \code{'wilcoxon'} for
#'   Wilcoxon rank sum test or \code{'shuffle'} for a label shuffling test.
#' #' @param randomizations Number of random shuffles, used for the label shuffling
#'   test (default = 100).
#' #' @param pcorr P-value correction method applied to the Wilcoxon rank sum test
#'   or label shuffling test results, as defined in the \code{p.adjust}
#'   function.
#' @param assayno An integer value. In case of multiple assays in a
#' #'    \code{SummarizedExperiment} input, the argument specifies the assay number
#'    to use for difference calculations.
#' #' @param verbose If \code{TRUE}, the function will print additional diagnostic
#'    messages.
#' @param ... Further arguments to be passed on for other methods.
#' @return A \code{data.frame} with the mean or median values of splicing
#'   diversity across sample categories and all samples, log2(fold change) of
#'   the two different conditions, raw and corrected p-values.
#' @import methods
#' @importFrom SummarizedExperiment SummarizedExperiment assays assay colData
#' @export
#' @details The function calculates diversity changes between two sample
#' conditions. It uses the output of the diversity calculation function, which
#' is a \code{SummarizedExperiment} object of splicing diversity values.
#' #' Additionally, it can use a \code{data.frame} as input, where the first column
#' contains gene names, and all additional columns contain splicing diversity
#' values for each sample. A vector of sample conditions also serves as input,
#' used for aggregating the samples by condition.
#'
#' It calculates the mean or median of the splicing diversity data per sample
#' #' condition, the difference of these values and the log2 fold change of the two
#' conditions. Furthermore, the user can select a statistical method to
#' #' calculate the significance of the changes. The p-values and adjusted p-values
#' are calculated using a Wilcoxon sum rank test or label shuffling test.
#'
#' The function will exclude genes of low sample size from the significance
#' calculation, depending on which statistical test is applied.
#'
#' @examples
#' # data.frame with splicing diversity values
#' x <- data.frame(Genes = letters[seq_len(10)], matrix(runif(80), ncol = 8))
#'
#' # sample categories
#' samples <- c(rep("Healthy", 4), rep("Pathogenic", 4))
#'
#' # To calculate the difference of splicing diversity changes between the
#' #' # 'Healthy' and 'Pathogenic' condition together with the significance values,
#' # using mean and Wilcoxon rank sum test, use:
#' #' calculate_difference(x, samples, control = "Healthy", method = "mean", test =
#' #' "wilcoxon")
calculate_difference <- function(x, samples, control, method = "mean",
                                 test = "wilcoxon", randomizations = 100,
                                 pcorr = "BH", assayno = 1, verbose = FALSE,
                                 ...) {
  # internal small helpers (kept here to avoid adding new files)
  .tsenat_prepare_df <- function(x, samples, assayno) {
    if (inherits(x, "RangedSummarizedExperiment") || inherits(x, "SummarizedExperiment")) {
      if (length(samples) != 1) stop("For SummarizedExperiment input, 'samples' must be the name of a single colData column.", call. = FALSE)
      samples <- SummarizedExperiment::colData(x)[[samples]]
      if (!is.numeric(assayno) || length(SummarizedExperiment::assays(x)) < assayno) stop("Invalid 'assayno'.", call. = FALSE)
      df <- as.data.frame(SummarizedExperiment::assays(x)[[assayno]])
      genes <- rownames(df)
      df <- cbind(genes = genes, df)
    } else {
      df <- as.data.frame(x)
    }
    list(df = df, samples = samples)
  }

  .tsenat_sample_matrix <- function(dfr) as.matrix(dfr[, -c(1, ncol(dfr) - 1, ncol(dfr)), drop = FALSE])

  # Validate input container
  # Reject matrices explicitly (tests expect this error for matrix input)
  if (is.matrix(x)) {
    stop("Input data type is not supported! Please use `?calculate_difference`
             to see the possible arguments and details.", call. = FALSE)
  }
  if (!(is.data.frame(x) || inherits(x, "RangedSummarizedExperiment") || inherits(x, "SummarizedExperiment"))) {
    stop("Input data type is not supported! Please use `?calculate_difference`
             to see the possible arguments and details.", call. = FALSE)
  }

  # prepare data.frame and sample vector (handles SummarizedExperiment)
  pd <- .tsenat_prepare_df(x, samples, assayno)
  df <- pd$df
  samples <- pd$samples

  # Basic consistency checks
  if (ncol(df) - 1 != length(samples)) stop("The number of columns in the data.frame is not equal to the number of          samples defined in the samples argument.", call. = FALSE)
  groups <- levels(as.factor(samples))
  # capture dots early and initialize 'case' to avoid unbound variable errors
  dots <- list(...)
  case <- if (!is.null(dots$case)) dots$case else NULL
  if (length(groups) > 2) {
    stop("The number of conditions are higher than two. Please use exactly two
             different sample conditions, e.g. healthy and pathogenic.", call. = FALSE)
  }
  if (length(groups) < 2) stop("The number of conditions are smaller than two. Please use exactly two
           different sample conditions, e.g. healthy and pathogenic.", call. = FALSE)
  if (!(control %in% samples)) stop("This control sample type cannot be found in your samples.", call. = FALSE)
  if (!(method %in% c("mean", "median"))) stop("Invalid method. Please use `?calculate_difference` to see the possible
           arguments and details.", call. = FALSE)
  if (!(test %in% c("wilcoxon", "shuffle"))) stop("Invalid test method. Please use `?calculate_difference` to see the
           possible arguments and details.", call. = FALSE)
  valid_pcorr <- c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
  if (!(pcorr %in% valid_pcorr)) stop("Invalid p-value correction method. Please use `?calculate_difference` to see the
           possible arguments and details.", call. = FALSE)

  # Informational warnings about sample size
  tab <- table(samples)
  if (test == "wilcoxon") {
    if (randomizations != 100 && verbose) {
      message("Note: The 'randomizations' argument is an option for label shuffling,
                    it won't have any effect on the Wilcoxon rank sum test.")
    }
    if (any(tab < 3) || sum(tab) < 8) {
      warning("Low sample size. Wilcoxon rank sum test requires at least
            three samples in a given category and at least 8 samples overall for a
                          theoretical p-value smaller than 0.05.", call. = FALSE)
    }
  }
  if (test == "shuffle") {
    if (sum(tab) <= 5) {
      warning("Low sample size, not enough samples for label shuffling!", call. = FALSE)
    }
    if (sum(tab) > 5 && sum(tab) < 10) {
      warning("Low sample size, label shuffling might not give informative and
                    correct results.", call. = FALSE)
    }
  }

  # Count non-NA observations per gene per group
  idx1 <- which(samples == groups[1])
  idx2 <- which(samples == groups[2])
  df$cond_1 <- rowSums(!is.na(df[, idx1 + 1, drop = FALSE]))
  df$cond_2 <- rowSums(!is.na(df[, idx2 + 1, drop = FALSE]))

  # Partition genes by sufficient observations for the chosen test
  if (test == "wilcoxon") {
    keep_mask <- (df$cond_1 >= 3 & df$cond_2 >= 3 & (df$cond_1 + df$cond_2) >= 8)
  } else {
    keep_mask <- (df$cond_1 + df$cond_2) >= 5
  }

  df_keep <- df[keep_mask, , drop = FALSE]
  df_small <- df[!keep_mask, , drop = FALSE]

  result_list <- list()

  # Helper to extract the numeric matrix of sample columns (keeps original order)
  sample_matrix <- .tsenat_sample_matrix

  if (nrow(df_keep) > 0) {
    if (nrow(df_small) > 0 && verbose) message(sprintf("Note: %d genes excluded from testing due to low sample counts.", nrow(df_small)))
    ymat <- sample_matrix(df_keep)
    # p-value calculation
    if (test == "wilcoxon") {
      ptab <- wilcoxon(ymat, samples, ...)
      # wilcoxon should return a data.frame of p-values named appropriately
    } else {
      ptab <- label_shuffling(ymat, samples, control, method, randomizations)
    }
    result_list$tested <- data.frame(genes = df_keep[, 1], calculate_fc(ymat, samples, control, method), ptab, stringsAsFactors = FALSE)
  }

  if (nrow(df_small) > 0) {
    small_mat <- sample_matrix(df_small)
    result_list$small <- data.frame(genes = df_small[, 1], calculate_fc(small_mat, samples, control, method), raw_p_values = NA, adjusted_p_values = NA, stringsAsFactors = FALSE)
  }

  # Combine results preserving tested rows first
  if (length(result_list) == 0) {
    return(data.frame())
  }
  res <- do.call(rbind, result_list)
  rownames(res) <- NULL
  return(res)
}


#' Linear-model interaction test for Tsallis entropy
#'
#' For each gene, fit a linear model of the form `entropy ~ q * group` and
#' extract the p-value for the interaction term (whether the effect of `q`
#' differs between groups). The function expects a `SummarizedExperiment`
#' produced by `calculate_diversity()` when multiple `q` values have been
#' computed (column names contain `_q=`).
#'
#' @param se A `SummarizedExperiment` containing a `diversity` assay produced
#'   by `calculate_diversity(..., q = <vector>)`.
#' @param sample_type_col Optional column name in `colData(se)` that contains
#'   a grouping factor for samples (character). If `NULL`, the function will
#'   attempt to infer group from column names (suffix `_N` interpreted as
#'   "Normal").
#' @param min_obs Minimum number of non-NA observations required to fit a
#'   model for a gene (default: 10).
#' @param method Modeling method to use for interaction testing: one of
#'   \code{c("linear", "gam", "fpca")}.
#' #' @param nthreads Number of threads (mc.cores) to use for parallel processing
#' #' (default: 1).
#' #' @param assay_name Name of the assay in the SummarizedExperiment to use
#' #' (default: "diversity").
#' @return A data.frame with columns `gene`, `p_interaction`, and
#'   `adj_p_interaction`, ordered by ascending `p_interaction`.
#' @export
#' @examples
#' data("tcga_brca_luma_dataset", package = "TSENAT")
#' rc <- as.matrix(tcga_brca_luma_dataset[1:20, -1, drop = FALSE])
#' gs <- tcga_brca_luma_dataset$genes[1:20]
#' se <- calculate_diversity(rc, gs, q = c(0.1, 1), norm = TRUE)
#' calculate_lm_interaction(se)
calculate_lm_interaction <- function(se, sample_type_col = NULL, min_obs = 10,
                                     method = c("linear", "gam", "fpca"), nthreads = 1,
                                     assay_name = "diversity") {
  method <- match.arg(method)
  message("[calculate_lm_interaction] method=", method)
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) stop("SummarizedExperiment required")

  mat <- SummarizedExperiment::assay(se, assay_name)
  if (is.null(mat)) stop(sprintf("`%s` assay not found in provided SummarizedExperiment", assay_name))

  sample_q <- colnames(mat)
  if (is.null(sample_q) || length(sample_q) == 0) stop("No column names found on diversity assay")

  # parse sample names and q values from column names like 'Sample_q=0.01'
  sample_names <- sub("_q=.*", "", sample_q)
  q_vals <- as.numeric(sub(".*_q=", "", sample_q))
  if (all(is.na(q_vals))) stop("Could not parse q values from column names; expected pattern '_q=' in names")

  # determine group for each sample
  if (!is.null(sample_type_col) && sample_type_col %in% colnames(SummarizedExperiment::colData(se))) {
    st <- as.character(SummarizedExperiment::colData(se)[, sample_type_col])
    names(st) <- SummarizedExperiment::colData(se)$samples %||% colnames(mat)
    # when user provides a sample_type_col, index by the sample names (strip q suffix)
    group_vec <- unname(st[sample_names])
  } else {
    # infer group from sample names using helper (supports TCGA barcodes and _N/_T suffixes)
    group_vec <- infer_sample_group(sample_names)
  }

  message("[calculate_lm_interaction] parsed samples and groups; starting per-gene fits")
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
      return(data.frame(gene = g, p_interaction = p_interaction, stringsAsFactors = FALSE))
    }

    if (method == "gam") {
      if (!requireNamespace("mgcv", quietly = TRUE)) stop("Package 'mgcv' is required for method = 'gam'")
      # choose smoothing basis dimension k based on number of unique q values
      uq_len <- length(unique(na.omit(q_vals)))
      k_q <- max(2, min(10, uq_len - 1))
      fit_null <- try(mgcv::gam(entropy ~ group + s(q, k = k_q), data = df), silent = TRUE)
      fit_alt <- try(mgcv::gam(entropy ~ group + s(q, by = group, k = k_q), data = df), silent = TRUE)
      if (inherits(fit_null, "try-error") || inherits(fit_alt, "try-error")) {
        return(NULL)
      }
      an <- try(mgcv::anova.gam(fit_null, fit_alt, test = "F"), silent = TRUE)
      if (inherits(an, "try-error")) {
        return(NULL)
      }
      # anova.gam returns a table; p-value typically in second row, column 'Pr(F)'
      p_interaction <- NA_real_
      if (nrow(an) >= 2) {
        if ("Pr(F)" %in% colnames(an)) {
          p_interaction <- an[2, "Pr(F)"]
        } else if ("Pr(>F)" %in% colnames(an)) {
          p_interaction <- an[2, "Pr(>F)"]
        } else if ("p-value" %in% colnames(an)) p_interaction <- an[2, "p-value"]
      }
      return(data.frame(gene = g, p_interaction = p_interaction, stringsAsFactors = FALSE))
    }

    if (method == "fpca") {
      message("[fpca] processing gene: ", g)
      # Build per-sample curves across q: rows = samples, cols = unique q values
      uq <- sort(unique(q_vals))
      samples_u <- unique(sample_names)
      message(sprintf("[fpca] uq=%s samples_u=%s", paste(uq, collapse = ","), paste(samples_u, collapse = ",")))
      curve_mat <- matrix(NA_real_, nrow = length(samples_u), ncol = length(uq))
      # assign rownames defensively; avoid assigning colnames to prevent dimname mismatch
      if (length(samples_u) > 0) rownames(curve_mat) <- samples_u
      message(sprintf("[fpca] created curve_mat with dims: %s", paste(dim(curve_mat), collapse = ",")))
      for (i in seq_along(sample_names)) {
        s <- sample_names[i]
        qv <- q_vals[i]
        qi <- match(qv, uq)
        message("[fpca] i=", i, " s=", s, " qv=", qv, " qi=", qi)
        if (is.na(qi)) next
        # protect in case sample name not present in rownames
        if (s %in% rownames(curve_mat)) curve_mat[s, qi] <- as.numeric(mat[g, i])
      }
      # keep samples with at least half of q points present
      good_rows <- which(rowSums(!is.na(curve_mat)) >= max(2, ceiling(ncol(curve_mat) / 2)))
      if (length(good_rows) < min_obs) {
        return(NULL)
      }
      mat_sub <- curve_mat[good_rows, , drop = FALSE]
      # simple imputation for remaining NAs using column means
      col_means <- apply(mat_sub, 2, function(col) mean(col, na.rm = TRUE))
      for (r in seq_len(nrow(mat_sub))) mat_sub[r, is.na(mat_sub[r, ])] <- col_means[is.na(mat_sub[r, ])]
      # perform PCA across q (observations = samples)
      pca <- try(stats::prcomp(mat_sub, center = TRUE, scale. = FALSE), silent = TRUE)
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
      return(data.frame(gene = g, p_interaction = pval, stringsAsFactors = FALSE))
    }
    return(NULL)
  }

  if (nthreads > 1 && .Platform$OS.type == "unix") {
    if (!requireNamespace("parallel", quietly = TRUE)) stop("parallel package required for multi-threading")
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
  res <- res[order(res$p_interaction), , drop = FALSE]
  rownames(res) <- NULL
  return(res)
}

# small helper (replacement for `%||%`) to provide default when NULL
`%||%` <- function(a, b) if (is.null(a)) b else a
