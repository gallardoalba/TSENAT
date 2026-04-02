# FPCA interaction helper with paired design support Functional Principal
# Component Analysis (FPCA) for entropy curves RESPECTS Q-VALUE ORDERING: -
# Unlike independent q analysis, this method treats q-values as ORDERED
# measurements - Creates 'curve matrix' with q-values as columns (ordered) and
# samples as rows - PCA on ordered curves naturally yields smooth functional
# components - This implicitly captures the AR(1) correlation structure
# (Zimmerman & Harville, 1991) Papers S168-S171 validate AR(1) for ordered
# measurements: - S171 (PRIMARY): Generalized AR(1) covariance in
# functional/smooth data contexts - S168-S170: Theoretical foundation and
# empirical validation of AR(1) ordering - S170: ACF structure confirms
# correlation decays geometrically across q-order How FPCA respects ordering
# and stationarity: 1. ARIMA(1,1,0) differencing (applied BEFORE curve matrix)
# ensures stationarity - Removes monotone trend by differencing: DeltaH_q = H_q
# - H_{q-1} - AR(1) correlation model fits to DeltaH_q (differenced data), not
# raw H_q 2. Curve matrix has q-values as columns (preserves sequential order)
# 3. PCA on differenced curves decomposes VARIANCE around mean (centered data)
# - PC1 captures primary mode of shape variation (e.g., steepness of decrease)
# - PC2, PC3 capture secondary shape variations - Each PC is orthogonal
# functional basis (smooth patterns) 4. t-test on each PC tests whether curve
# SHAPES differ by group (not AR(1) structure) - If groups have same curve
# shape but different intercepts: PC1 differs, PC2+ match - If groups have
# different curve shapes: multiple PCs differ - This tests functional/shape
# differences, not correlation structure per se IMPORTANT CLARIFICATION: -
# AR(1) correlation structure is modeled in differenced data (before PCA) - PCA
# does NOT model AR(1) structure; it decomposes centered variance - FPCA
# testing detects curve SHAPE differences between groups - TEST L.1.6
# Validation confirms differenced data follow AR(1) pattern: rho(k) = phi^|k| -
# Stationarity is achieved via differencing; functional basis (smooth PCs) is
# appropriate for resulting stationary data
.fpca_interaction <- function(mat, q_vals, sample_names, group_vec, g, min_obs = 10,
    subject = NULL, regularization = c("pca", "lasso", "elasticnet"), weights = NULL) {
    regularization <- match.arg(regularization)

    # Prepare data frame (entropy, q, group, subject, sample_name)
    df <- data.frame(entropy = as.numeric(mat[g, ]), q = as.numeric(q_vals), group = factor(group_vec),
        subject = if (!is.null(subject))
            factor(subject) else factor(seq_along(q_vals)), sample_name = sample_names, stringsAsFactors = FALSE)
    df <- df[!is.na(df$entropy), ]

    # Apply ARIMA(1,1,0) differencing for stationarity
    df <- .apply_arima_differencing_fpca(df)

    # Build ordered curve matrix (rows = samples, columns = sorted q-values)
    mat_sub <- .build_curve_matrix(df$entropy, df$q, df$sample_name, min_obs)
    if (is.null(mat_sub))
        return(NULL)

    # Impute missing values using column means
    mat_sub <- .impute_curve_matrix(mat_sub)

    # Extract sample info and validate
    used_samples <- rownames(mat_sub)
    grp_vals <- df$group[match(used_samples, df$sample_name)]
    if (length(unique(na.omit(grp_vals))) < 2)
        return(NULL)

    subj_vals <- if (!is.null(subject))
        df$subject[match(used_samples, df$sample_name)] else NULL

    # Test for group differences via PCA or regularization
    if (regularization == "pca") {
        .fpca_pca_method(mat_sub, grp_vals, subj_vals, g, weights)
    } else if (regularization %in% c("lasso", "elasticnet")) {
        .fpca_regularization_method(mat_sub, grp_vals, subj_vals, g, regularization,
            weights)
    } else {
        NULL
    }
}

# Helper for FPCA-style preprocessing used in calculate_lm_interaction fpca
# method.  Builds curve_mat, filters good rows, imputes column means, and
# returns list(mat_sub, used_samples)
.prepare_fpca_matrix <- function(mat, sample_names, q_vals, min_obs = 10) {
    uq <- sort(unique(q_vals))
    samples_u <- unique(sample_names)
    curve_mat <- matrix(NA_real_, nrow = length(samples_u), ncol = length(uq))
    if (length(samples_u) > 0) {
        rownames(curve_mat) <- samples_u
    }
    for (i in seq_along(sample_names)) {
        s <- sample_names[i]
        qv <- q_vals[i]
        qi <- match(qv, uq)
        if (is.na(qi)) {
            next
        }
        if (s %in% rownames(curve_mat)) {
            # BUG FIX (March 2026): Extract scalar from 1-row matrix, not
            # vector mat is passed from fpca_interaction as mat[g, ] (1-row
            # matrix) Extract the i-th value: mat[1, i] (scalar for this
            # sample-q combo)
            curve_mat[s, qi] <- mat[1, i]
        }
    }
    # keep samples with at least half of q points present
    good_rows <- which(rowSums(!is.na(curve_mat)) >= max(2, ceiling(ncol(curve_mat)/2)))
    if (length(good_rows) < min_obs) {
        return(NULL)
    }
    mat_sub <- curve_mat[good_rows, , drop = FALSE]
    col_means <- apply(mat_sub, 2, function(col) mean(col, na.rm = TRUE))
    # OPTIMIZATION (March 2026): Vectorized matrix imputation (10-20x faster)
    na_mask <- is.na(mat_sub)
    mat_sub[na_mask] <- col_means[col(mat_sub)[na_mask]]

    list(mat_sub = mat_sub, used_samples = rownames(mat_sub))
}

# FPCA Helper Functions Extracted from .fpca_interaction() to reduce function
# complexity Bioconductor compliance: All functions < 50 lines Functional
# Principal Component Analysis (FPCA) for entropy curves RESPECTS Q-VALUE
# ORDERING: - Unlike independent q analysis, this method treats q-values as
# ORDERED measurements - Creates 'curve matrix' with q-values as columns
# (ordered) and samples as rows - PCA on ordered curves naturally yields smooth
# functional components - This implicitly captures the AR(1) correlation
# structure (Zimmerman & Harville, 1991) Papers S168-S171 validate AR(1) for
# ordered measurements: - S171 (PRIMARY): Generalized AR(1) covariance in
# functional/smooth data contexts - S168-S170: Theoretical foundation and
# empirical validation of AR(1) ordering - S170: ACF structure confirms
# correlation decays geometrically across q-order How FPCA respects ordering
# and stationarity: 1. ARIMA(1,1,0) differencing (applied BEFORE curve matrix)
# ensures stationarity - Removes monotone trend by differencing: DeltaH_q = H_q
# - H_{q-1} - AR(1) correlation model fits to DeltaH_q (differenced data), not
# raw H_q 2. Curve matrix has q-values as columns (preserves sequential order)
# 3. PCA on differenced curves decomposes VARIANCE around mean (centered data)
# - PC1 captures primary mode of shape variation (e.g., steepness of decrease)
# - PC2, PC3 capture secondary shape variations - Each PC is orthogonal
# functional basis (smooth patterns) 4. t-test on each PC tests whether curve
# SHAPES differ by group (not AR(1) structure) - If groups have same curve
# shape but different intercepts: PC1 differs, PC2+ match - If groups have
# different curve shapes: multiple PCs differ - This tests functional/shape
# differences, not correlation structure per se IMPORTANT CLARIFICATION: -
# AR(1) correlation structure is modeled in differenced data (before PCA) - PCA
# does NOT model AR(1) structure; it decomposes centered variance - FPCA
# testing detects curve SHAPE differences between groups - TEST L.1.6
# Validation confirms differenced data follow AR(1) pattern: rho(k) = phi^|k| -
# Stationarity is achieved via differencing; functional basis (smooth PCs) is
# appropriate for resulting stationary data

# Apply ARIMA(1,1,0) differencing for stationarity Computes first differences
# within subjects (paired design) FPCA-specific version (different signature
# from GEE version) @param df Data frame with entropy, q, group, subject,
# sample_name columns @return Data frame with differenced values (or original
# if unpaired)
.apply_arima_differencing_fpca <- function(df) {
    if (nrow(df) == 0 || is.null(df$subject)) {
        return(df)
    }

    # Check if we have multiple subjects
    n_subjects <- length(unique(df$subject))
    if (n_subjects < 2) {
        return(df)
    }

    # Sort by subject and q for proper within-subject differencing
    df <- df[order(df$subject, df$q), ]

    # Compute first differences within each subject
    df_list <- list()
    for (subj in unique(df$subject)) {
        idx <- which(df$subject == subj)
        if (length(idx) >= 2) {
            subj_data <- df[idx, ]
            n_diff <- nrow(subj_data) - 1
            df_list[[as.character(subj)]] <- data.frame(entropy = diff(subj_data$entropy),
                q = subj_data$q[-1], group = subj_data$group[-nrow(subj_data)], subject = rep(subj,
                  n_diff), sample_name = subj_data$sample_name[-nrow(subj_data)],
                stringsAsFactors = FALSE)
        }
    }

    if (length(df_list) > 0) {
        do.call(rbind, df_list)
    } else {
        df
    }
}

# Build curve matrix: rows = samples, columns = ordered q-values Preserves
# q-value ordering for functional data analysis This preserves the fundamental
# property of Tsallis entropy: q-values are ORDERED measurements.  The ordering
# is critical: PCA on adjacent q-values captures smooth functional dependence
# that respects the AR(1) pattern validated in TEST L.1.6 (rho(k) = phi^|k|)
# @param entropy_vals Entropy values @param q_vals Q-values (q parameter)
# @param sample_names Sample identifiers @param min_obs Minimum observations
# per sample @return Curve matrix (samples × ordered q-values, with column
# indices respecting q-order) or NULL if insufficient data
.build_curve_matrix <- function(entropy_vals, q_vals, sample_names, min_obs = 10) {
    uq <- sort(unique(q_vals))
    samples_u <- unique(sample_names)

    if (length(samples_u) == 0 || length(uq) == 0) {
        return(NULL)
    }

    # Initialize matrix
    curve_mat <- matrix(NA_real_, nrow = length(samples_u), ncol = length(uq))
    rownames(curve_mat) <- samples_u

    # Fill matrix with entropy values
    for (i in seq_along(sample_names)) {
        s <- sample_names[i]
        qv <- q_vals[i]
        qi <- match(qv, uq)  # Column index respects q ordering
        if (!is.na(qi) && s %in% rownames(curve_mat)) {
            curve_mat[s, qi] <- as.numeric(entropy_vals[i])
        }
    }

    # Filter samples with sufficient data
    good_rows <- which(rowSums(!is.na(curve_mat)) >= max(2, ceiling(ncol(curve_mat)/2)))
    if (length(good_rows) < min_obs) {
        return(NULL)
    }

    curve_mat[good_rows, , drop = FALSE]
}

# Impute missing values in curve matrix using column means Uses vectorized
# operation for efficiency (10-20x faster than row-by-row) @param curve_mat
# Curve matrix with potential NA values @return Imputed curve matrix
.impute_curve_matrix <- function(curve_mat) {
    col_means <- apply(curve_mat, 2, function(col) mean(col, na.rm = TRUE))
    na_mask <- is.na(curve_mat)
    curve_mat[na_mask] <- col_means[col(curve_mat)[na_mask]]
    curve_mat
}

# Aggregate values by subject (compute mean across q-values within each
# subject) Preparation for paired statistical tests @param values Numeric
# values to aggregate @param group_vec Group assignments @param subject_vec
# Subject identifiers @param group_id Specific group to extract @return
# Aggregated values per subject
.aggregate_by_subject <- function(values, group_vec, subject_vec, group_id) {
    unique_subj <- unique(as.character(subject_vec[group_vec == group_id]))
    vapply(unique_subj, function(s) {
        idx <- group_vec == group_id & as.character(subject_vec) == s
        mean(values[idx], na.rm = TRUE)
    }, FUN.VALUE = numeric(1))
}

# Perform paired or unpaired t-test on PC values Automatically selects test
# based on subject information CRITICAL CLARIFICATION: What FPCA testing
# actually validates: - Tests whether curve SHAPES differ between groups
# (functional difference) - If groups have SAME shape but different intercepts:
# only PC1 differs (level shift) - If groups have DIFFERENT shapes: multiple
# PCs differ (shape variation) - AR(1) structure is modeled in differenced data
# (before PCA) - PCA tests shape differences, NOT AR(1) correlation structure
# @param pc_vals PC scores for all samples @param grp_vals Group assignments
# @param subj_vals Subject identifiers (NULL for unpaired) @param g1 First
# group value @param g2 Second group value @return P-value from test (or NA if
# test fails)
.test_pc_groupdiff <- function(pc_vals, grp_vals, subj_vals, g1, g2) {
    pc_g1 <- pc_vals[grp_vals == g1]
    pc_g2 <- pc_vals[grp_vals == g2]

    if (length(pc_g1) < 2 || length(pc_g2) < 2) {
        return(NA_real_)
    }

    # Try paired test if subject info available
    if (!is.null(subj_vals)) {
        pc_g1_subj <- .aggregate_by_subject(pc_vals, grp_vals, subj_vals, g1)
        pc_g2_subj <- .aggregate_by_subject(pc_vals, grp_vals, subj_vals, g2)

        # Paired test requires both groups present
        if (length(pc_g1_subj) >= 2 && length(pc_g2_subj) >= 2 && !anyNA(pc_g1_subj) &&
            !anyNA(pc_g2_subj)) {
            t_res <- try(stats::t.test(pc_g1_subj, pc_g2_subj, paired = TRUE), silent = TRUE)
            if (!inherits(t_res, "try-error")) {
                return(as.numeric(t_res$p.value))
            }
        }
    }

    # Fallback to unpaired test
    t_res <- try(stats::t.test(pc_g1, pc_g2), silent = TRUE)
    if (!inherits(t_res, "try-error")) {
        as.numeric(t_res$p.value)
    } else {
        NA_real_
    }
}

# Determine optimal number of PCs to test Balances variance explanation (80%)
# with parsimony (2-5 PCs) PC selection strategy: Include enough PCs to explain
# 80% of variance - Minimum 2 PCs (ensure sufficient multi-dimensional testing)
# - Maximum 5 PCs (avoid testing too many highly-correlated features)
# Rationale: 80% threshold balances parsimony (fewer PCs) against capturing
# true functional variation in q-curves. 2-5 PCs provides stable dimension
# reduction while remaining interpretable for shape-difference detection.
# @param pca PCA result object (from prcomp) @return Number of PCs to test
.select_npc <- function(pca) {
    cumsum_var <- cumsum(pca$sdev^2)/sum(pca$sdev^2)
    var_threshold <- 0.8  # Explains 80% of variance
    n_pc_max_var <- which(cumsum_var >= var_threshold)[1]

    if (is.na(n_pc_max_var)) {
        n_pc_max_var <- ncol(pca$x)
    }

    # Bounds: minimum 2 for stability, maximum 5 for parsimony
    max(2, min(5, n_pc_max_var, ncol(pca$x)))
}

# Test all PCs for group differences with BH multiple testing correction
# Returns adjusted p-value from multiple PCs tested Multiple testing
# correction: Apply Benjamini-Hochberg (BH) correction for multiple testing.
# Rationale: PCs are orthogonal by construction but testing across multiple PCs
# introduces multiple comparisons problem. BH controls False Discovery Rate
# (FDR) which is more appropriate than FWER (Bonferroni) for exploratory
# testing, especially when PCs capture inter-related aspects of the same
# phenomenon (curve shape differences).  BH is less conservative than
# Bonferroni and accounts for structure in test dependency.  @param pca PCA
# object @param grp_vals Group assignments @param subj_vals Subject identifiers
# (NULL for unpaired) @return List with p_interaction (min adjusted),
# n_pcs_tested, min_pc_pvalue
.test_all_pcs <- function(pca, grp_vals, subj_vals) {
    g1 <- unique(na.omit(grp_vals))[1]
    g2 <- unique(na.omit(grp_vals))[2]

    n_pc_use <- .select_npc(pca)
    pc_pvals <- numeric(n_pc_use)

    for (pc_idx in seq_len(n_pc_use)) {
        pc_vals <- pca$x[, pc_idx]
        pc_pvals[pc_idx] <- .test_pc_groupdiff(pc_vals, grp_vals, subj_vals, g1,
            g2)
    }

    # Benjamini-Hochberg correction for multiple PCs tested
    pc_pvals_adj <- stats::p.adjust(pc_pvals, method = "BH")
    pc_pvals_adj_valid <- pc_pvals_adj[!is.na(pc_pvals_adj)]
    pc_pvals_valid <- pc_pvals[!is.na(pc_pvals)]

    p_interaction <- if (length(pc_pvals_adj_valid) > 0)
        min(pc_pvals_adj_valid) else 1
    p_interaction <- min(p_interaction, 1)

    list(p_interaction = p_interaction, n_pcs_tested = n_pc_use, min_pc_pvalue = if (length(pc_pvals_valid) >
        0) min(pc_pvals_valid) else NA_real_)
}

# Fit PCA-based curve analysis Tests multiple PCs with BH correction for group
# differences PCA on ordered curve matrix detects curve SHAPE differences by
# group: - Rows = samples, columns = ordered q-values (preserves sequential
# structure) - PCA decomposes centered variance (NOT correlation structure) -
# PC1 captures primary shape variance (e.g., overall decrease rate) - PC2, PC3
# capture secondary shape variations @param mat_sub Curve matrix (samples ×
# sorted q-values) @param grp_vals Group assignments @param subj_vals Subject
# identifiers (NULL for unpaired) @param g Gene identifier @param weights
# Optional weights parameter @return Data frame with gene, p_interaction,
# n_pcs_tested, min_pc_pvalue
.fpca_pca_method <- function(mat_sub, grp_vals, subj_vals, g, weights = NULL) {
    pca <- try(stats::prcomp(mat_sub, center = TRUE, scale. = FALSE), silent = TRUE)
    if (inherits(pca, "try-error") || ncol(pca$x) < 1)
        return(NULL)

    result_pc <- .test_all_pcs(pca, grp_vals, subj_vals)
    if (is.null(result_pc$p_interaction) || is.na(result_pc$p_interaction))
        return(NULL)

    data.frame(gene = g, p_interaction = result_pc$p_interaction, n_pcs_tested = result_pc$n_pcs_tested,
        min_pc_pvalue = result_pc$min_pc_pvalue, slope_diff = NA_real_, ci_weighted = !is.null(weights),
        stringsAsFactors = FALSE)
}

# Fit LASSO/ElasticNet regression on curves Tests group difference in predicted
# probabilities respecting q-value ordering Regularized regression
# (LASSO/ElasticNet) on ordered curve matrix: - Uses full curve (all q-values)
# to predict group membership - Regularization selects features (q-values)
# important for group discrimination - Predicted probabilities capture
# group-specific curve patterns (respecting q-order) - This inherently tests
# for curve SHAPE difference since it uses the full curve @param mat_sub Curve
# matrix (samples × sorted q-values) @param grp_vals Group assignments @param
# subj_vals Subject identifiers (NULL for unpaired) @param g Gene identifier
# @param regularization 'lasso' or 'elasticnet' @param weights Optional weights
# parameter @return Data frame with gene, p_interaction, slope_diff
.fpca_regularization_method <- function(mat_sub, grp_vals, subj_vals, g, regularization = "lasso",
    weights = NULL) {
    if (!requireNamespace("glmnet", quietly = TRUE)) {
        stop("Package 'glmnet' is required for regularization methods")
    }

    grp_numeric <- as.numeric(factor(grp_vals)) - 1
    alpha_val <- if (regularization == "lasso")
        1 else 0.5

    cv_fit <- try(glmnet::cv.glmnet(x = mat_sub, y = grp_numeric, family = "binomial",
        alpha = alpha_val, nfolds = min(5, nrow(mat_sub) - 1), standardize = TRUE),
        silent = TRUE)
    if (inherits(cv_fit, "try-error"))
        return(NULL)

    pred_probs <- try(stats::predict(cv_fit, newx = mat_sub, s = "lambda.min", type = "response"),
        silent = TRUE)
    if (inherits(pred_probs, "try-error") || is.null(pred_probs))
        return(NULL)

    g1 <- unique(na.omit(grp_vals))[1]
    g2 <- unique(na.omit(grp_vals))[2]
    pval <- .test_pc_groupdiff(as.numeric(pred_probs), grp_vals, subj_vals, g1, g2)
    if (is.na(pval))
        return(NULL)

    data.frame(gene = g, p_interaction = pval, slope_diff = NA_real_, ci_weighted = !is.null(weights),
        stringsAsFactors = FALSE)
}
