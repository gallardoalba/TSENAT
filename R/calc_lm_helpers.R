# Summary reporting helper
.tsenat_report_fit_summary <- function(res, verbose = TRUE) {
    if (verbose && "fit_method" %in% colnames(res)) {
        total_genes <- nrow(res)
        fallback_mask <- !is.na(res$fit_method) & res$fit_method != "lmer"
        n_fallback <- sum(fallback_mask)
        if (n_fallback > 0) {
            tab <- table(res$fit_method[fallback_mask])
            tab_str <- paste(sprintf("%s=%d", names(tab), as.integer(tab)), collapse = ", ")
            message(sprintf("[calculate_lm_interaction] fallback fits used: %d/%d genes (%s)",
                n_fallback, total_genes, tab_str))
        }
        if ("singular" %in% colnames(res)) {
            n_sing <- sum(as.logical(res$singular), na.rm = TRUE)
            if (n_sing > 0)
                message(sprintf("[calculate_lm_interaction] singular fits detected: %d/%d genes",
                  n_sing, total_genes))
        }
    }
}
# GAM interaction helper
.tsenat_gam_interaction <- function(df, q_vals, g, min_obs = 10) {
    if (!requireNamespace("mgcv", quietly = TRUE)) {
        stop("Package 'mgcv' is required for method = 'gam'")
    }
    uq_len <- length(unique(na.omit(q_vals)))
    k_q <- max(2, min(10, uq_len - 1))
    fit_null <- try(mgcv::gam(entropy ~ group + s(q, k = k_q), data = df), silent = TRUE)
    fit_alt <- try(mgcv::gam(entropy ~ group + s(q, by = group, k = k_q), data = df),
        silent = TRUE)
    if (inherits(fit_null, "try-error") || inherits(fit_alt, "try-error"))
        return(NULL)
    an <- try(mgcv::anova.gam(fit_null, fit_alt, test = "F"), silent = TRUE)
    if (inherits(an, "try-error"))
        return(NULL)
    p_interaction <- NA_real_
    if (nrow(an) >= 2) {
        if ("Pr(F)" %in% colnames(an))
            p_interaction <- an[2, "Pr(F)"] else if ("Pr(>F)" %in% colnames(an))
            p_interaction <- an[2, "Pr(>F)"] else if ("p-value" %in% colnames(an))
            p_interaction <- an[2, "p-value"]
    }
    return(data.frame(gene = g, p_interaction = p_interaction, stringsAsFactors = FALSE))
}

# FPCA interaction helper
.tsenat_fpca_interaction <- function(mat, q_vals, sample_names, group_vec, g, min_obs = 10) {
    uq <- sort(unique(q_vals))
    samples_u <- unique(sample_names)
    curve_mat <- matrix(NA_real_, nrow = length(samples_u), ncol = length(uq))
    if (length(samples_u) > 0)
        rownames(curve_mat) <- samples_u
    for (i in seq_along(sample_names)) {
        s <- sample_names[i]
        qv <- q_vals[i]
        qi <- match(qv, uq)
        if (is.na(qi))
            next
        if (s %in% rownames(curve_mat))
            curve_mat[s, qi] <- as.numeric(mat[g, i])
    }
    good_rows <- which(rowSums(!is.na(curve_mat)) >= max(2, ceiling(ncol(curve_mat)/2)))
    if (length(good_rows) < min_obs)
        return(NULL)
    mat_sub <- curve_mat[good_rows, , drop = FALSE]
    col_means <- apply(mat_sub, 2, function(col) mean(col, na.rm = TRUE))
    for (r in seq_len(nrow(mat_sub))) mat_sub[r, is.na(mat_sub[r, ])] <- col_means[is.na(mat_sub[r,
        ])]
    pca <- try(stats::prcomp(mat_sub, center = TRUE, scale. = FALSE), silent = TRUE)
    if (inherits(pca, "try-error"))
        return(NULL)
    if (ncol(pca$x) < 1)
        return(NULL)
    pc1 <- pca$x[, 1]
    used_samples <- rownames(mat_sub)
    grp_vals <- group_vec[match(used_samples, sample_names)]
    if (length(unique(na.omit(grp_vals))) < 2)
        return(NULL)
    g1 <- unique(na.omit(grp_vals))[1]
    g2 <- unique(na.omit(grp_vals))[2]
    x1 <- pc1[grp_vals == g1]
    x2 <- pc1[grp_vals == g2]
    if (length(x1) < 2 || length(x2) < 2)
        return(NULL)
    t_res <- try(stats::t.test(x1, x2), silent = TRUE)
    if (inherits(t_res, "try-error"))
        return(NULL)
    pval <- as.numeric(t_res$p.value)
    return(data.frame(gene = g, p_interaction = pval, stringsAsFactors = FALSE))
}
## All helpers for calculate_lm_interaction

# Try lme4::lmer with multiple optimizers and controlled warnings.
.tsenat_try_lmer <- function(formula, data, suppress_lme4_warnings = TRUE, verbose = FALSE,
    mm_suppress_pattern = "boundary \\(singular\\) fit|Computed variance-covariance matrix problem|not a positive definite matrix") {
    if (!requireNamespace("lme4", quietly = TRUE)) {
        stop("Package 'lme4' is required for mixed-model fitting")
    }
    opts <- list(list(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+05)), list(optimizer = "nloptwrap",
        optCtrl = list(maxfun = 5e+05)))
    for (o in opts) {
        ctrl <- lme4::lmerControl(optimizer = o$optimizer, optCtrl = o$optCtrl)
        muffle_cond <- suppress_lme4_warnings || (!verbose)
        fit_try <- withCallingHandlers(try(lme4::lmer(formula, data = data, REML = FALSE,
            control = ctrl), silent = TRUE), warning = function(w) {
            if (muffle_cond && grepl(mm_suppress_pattern, conditionMessage(w), ignore.case = TRUE))
                invokeRestart("muffleWarning")
        }, message = function(m) {
            if (muffle_cond && grepl(mm_suppress_pattern, conditionMessage(m), ignore.case = TRUE))
                invokeRestart("muffleMessage")
        })
        if (!inherits(fit_try, "try-error")) {
            # check singularity if function available
            is_sing <- FALSE
            if (exists("isSingular", where = asNamespace("lme4"), inherits = FALSE)) {
                is_sing <- tryCatch(lme4::isSingular(fit_try, tol = 1e-04), error = function(e) FALSE)
            }
            attr(fit_try, "singular") <- is_sing
            return(fit_try)
        }
    }
    # all attempts failed
    return(structure("error", class = "try-error"))
}

# LM fallback helpers
.tsenat_try_lm_fallbacks <- function(df, verbose = FALSE) {
    fit0_lm <- try(stats::lm(entropy ~ q + group + subject, data = df), silent = TRUE)
    fit1_lm <- try(stats::lm(entropy ~ q * group + subject, data = df), silent = TRUE)
    if (!inherits(fit0_lm, "try-error") && !inherits(fit1_lm, "try-error")) {
        return(list(fit0 = fit0_lm, fit1 = fit1_lm, method = "lm_subject"))
    }
    fit0_lm2 <- try(stats::lm(entropy ~ q + group, data = df), silent = TRUE)
    fit1_lm2 <- try(stats::lm(entropy ~ q * group, data = df), silent = TRUE)
    if (!inherits(fit0_lm2, "try-error") && !inherits(fit1_lm2, "try-error")) {
        return(list(fit0 = fit0_lm2, fit1 = fit1_lm2, method = "lm_nosubject"))
    }
    return(NULL)
}

# LRT p-value extraction
.tsenat_extract_lrt_p <- function(fit0, fit1) {
    an <- try(stats::anova(fit0, fit1), silent = TRUE)
    if (!inherits(an, "try-error") && nrow(an) >= 2) {
        pcol <- grep("Pr\\(>F\\)|Pr\\(>Chisq\\)|Pr\\(>Chi\\)", colnames(an), value = TRUE)
        if (length(pcol) == 0)
            return(as.numeric(an[2, ncol(an)])) else return(as.numeric(an[2, pcol[1]]))
    }
    return(NA_real_)
}

# Satterthwaite p-value extraction
.tsenat_extract_satterthwaite_p <- function(fit1, fallback_lm = NULL) {
    if (!is.null(fallback_lm)) {
        coefs <- try(summary(fallback_lm$fit1)$coefficients, silent = TRUE)
        if (!inherits(coefs, "try-error")) {
            ia_idx <- grep("^q:group", rownames(coefs))
            if (length(ia_idx) > 0)
                return(coefs[ia_idx[1], "Pr(>|t|)"])
        }
        return(NA_real_)
    }
    if (requireNamespace("lmerTest", quietly = TRUE) && inherits(fit1, "lmerMod") &&
        !isTRUE(attr(fit1, "singular"))) {
        fit_lt <- try(lmerTest::lmer(stats::formula(fit1), data = stats::model.frame(fit1),
            REML = FALSE), silent = TRUE)
        if (!inherits(fit_lt, "try-error")) {
            coefs <- summary(fit_lt)$coefficients
            ia_idx <- grep("^q:group", rownames(coefs))
            if (length(ia_idx) > 0)
                return(coefs[ia_idx[1], "Pr(>|t|)"])
        }
    }
    return(NA_real_)
}

# FPCA matrix preparation
.tsenat_prepare_fpca_matrix <- function(mat, min_frac = 0.01) {
    if (!is.matrix(mat))
        mat <- as.matrix(mat)
    row_vars <- apply(mat, 1, stats::var, na.rm = TRUE)
    keep <- row_vars > (min_frac * max(row_vars, na.rm = TRUE))
    if (sum(keep) == 0)
        keep <- rep(TRUE, nrow(mat))
    m2 <- mat[keep, , drop = FALSE]
    m2 <- t(scale(t(m2)))
    return(list(mat = m2, keep = keep))
}

## Helpers for calculate_lm_interaction fallbacks, LRT and Satterthwaite
## p-values
.tsenat_try_lm_fallbacks <- function(df, verbose = FALSE) {
    # try lm with subject as fixed effect
    fit0_lm <- try(stats::lm(entropy ~ q + group + subject, data = df), silent = TRUE)
    fit1_lm <- try(stats::lm(entropy ~ q * group + subject, data = df), silent = TRUE)
    if (!inherits(fit0_lm, "try-error") && !inherits(fit1_lm, "try-error")) {
        return(list(fit0 = fit0_lm, fit1 = fit1_lm, method = "lm_subject"))
    }
    # last resort: drop subject, plain lm
    fit0_lm2 <- try(stats::lm(entropy ~ q + group, data = df), silent = TRUE)
    fit1_lm2 <- try(stats::lm(entropy ~ q * group, data = df), silent = TRUE)
    if (!inherits(fit0_lm2, "try-error") && !inherits(fit1_lm2, "try-error")) {
        return(list(fit0 = fit0_lm2, fit1 = fit1_lm2, method = "lm_nosubject"))
    }
    return(NULL)
}

.tsenat_extract_lrt_p <- function(fit0, fit1) {
    an <- try(stats::anova(fit0, fit1), silent = TRUE)
    if (!inherits(an, "try-error") && nrow(an) >= 2) {
        pcol <- grep("Pr\\(>F\\)|Pr\\(>Chisq\\)|Pr\\(>Chi\\)", colnames(an), value = TRUE)
        if (length(pcol) == 0)
            return(as.numeric(an[2, ncol(an)])) else return(as.numeric(an[2, pcol[1]]))
    }
    return(NA_real_)
}

.tsenat_extract_satterthwaite_p <- function(fit1, fallback_lm = NULL) {
    # if lmerTest available and fit1 is not singular use that
    if (!is.null(fallback_lm)) {
        coefs <- try(summary(fallback_lm$fit1)$coefficients, silent = TRUE)
        if (!inherits(coefs, "try-error")) {
            ia_idx <- grep("^q:group", rownames(coefs))
            if (length(ia_idx) > 0)
                return(coefs[ia_idx[1], "Pr(>|t|)"])
        }
        return(NA_real_)
    }
    if (requireNamespace("lmerTest", quietly = TRUE) && inherits(fit1, "lmerMod") &&
        !isTRUE(attr(fit1, "singular"))) {
        fit_lt <- try(lmerTest::lmer(stats::formula(fit1), data = stats::model.frame(fit1),
            REML = FALSE), silent = TRUE)
        if (!inherits(fit_lt, "try-error")) {
            coefs <- summary(fit_lt)$coefficients
            ia_idx <- grep("^q:group", rownames(coefs))
            if (length(ia_idx) > 0)
                return(coefs[ia_idx[1], "Pr(>|t|)"])
        }
    }
    return(NA_real_)
}

.tsenat_prepare_fpca_matrix <- function(mat, min_frac = 0.01) {
    # prepare matrix for FPCA: center, scale, and drop near-constant rows
    if (!is.matrix(mat))
        mat <- as.matrix(mat)
    row_vars <- apply(mat, 1, stats::var, na.rm = TRUE)
    keep <- row_vars > (min_frac * max(row_vars, na.rm = TRUE))
    if (sum(keep) == 0)
        keep <- rep(TRUE, nrow(mat))
    m2 <- mat[keep, , drop = FALSE]
    m2 <- t(scale(t(m2)))
    return(list(mat = m2, keep = keep))
}
## Additional helpers for calculate_lm_interaction

# Try simple lm fallbacks when lmer fails. Returns list(fit0, fit1, method)
.tsenat_try_lm_fallbacks <- function(df, verbose = FALSE) {
    # try lm with subject as fixed effect
    fit0_lm <- try(stats::lm(entropy ~ q + group + subject, data = df), silent = TRUE)
    fit1_lm <- try(stats::lm(entropy ~ q * group + subject, data = df), silent = TRUE)
    if (!inherits(fit0_lm, "try-error") && !inherits(fit1_lm, "try-error")) {
        return(list(fit0 = fit0_lm, fit1 = fit1_lm, method = "lm_subject"))
    }
    # last resort: drop subject
    fit0_lm2 <- try(stats::lm(entropy ~ q + group, data = df), silent = TRUE)
    fit1_lm2 <- try(stats::lm(entropy ~ q * group, data = df), silent = TRUE)
    if (!inherits(fit0_lm2, "try-error") && !inherits(fit1_lm2, "try-error")) {
        return(list(fit0 = fit0_lm2, fit1 = fit1_lm2, method = "lm_nosubject"))
    }
    # indicate failure
    return(NULL)
}

# Extract LRT p-value from two fitted models (lm or lmer). Returns NA_real_ on
# failure.
.tsenat_extract_lrt_p <- function(fit0, fit1) {
    an <- try(stats::anova(fit0, fit1), silent = TRUE)
    if (inherits(an, "try-error") || nrow(an) < 2)
        return(NA_real_)
    pcol <- grep("Pr\\(>F\\)|Pr\\(>Chisq\\)|Pr\\(>Chi\\)", colnames(an), value = TRUE)
    if (length(pcol) == 0) {
        return(as.numeric(an[2, ncol(an)]))
    }
    return(as.numeric(an[2, pcol[1]]))
}

# Extract Satterthwaite p-value for q:group interaction from lmer/lmerTest or
# fallback lm.
.tsenat_extract_satterthwaite_p <- function(fit1, fallback_lm = NULL) {
    # prefer lmerTest when available
    if (!is.null(fallback_lm)) {
        coefs <- summary(fallback_lm$fit1)$coefficients
        ia_idx <- grep("^q:group", rownames(coefs))
        if (length(ia_idx) > 0)
            return(coefs[ia_idx[1], "Pr(>|t|)"])
        return(NA_real_)
    }
    if (requireNamespace("lmerTest", quietly = TRUE) && inherits(fit1, "lmerMod")) {
        fit_lt <- try(lmerTest::lmer(as.formula(formula(fit1)), data = model.frame(fit1),
            REML = FALSE), silent = TRUE)
        if (!inherits(fit_lt, "try-error")) {
            coefs <- summary(fit_lt)$coefficients
            ia_idx <- grep("^q:group", rownames(coefs))
            if (length(ia_idx) > 0)
                return(coefs[ia_idx[1], "Pr(>|t|)"])
        }
    }
    return(NA_real_)
}

# Helper for FPCA-style preprocessing used in calculate_lm_interaction fpca
# method.  Builds curve_mat, filters good rows, imputes column means, and
# returns list(mat_sub, used_samples)
.tsenat_prepare_fpca_matrix <- function(mat, sample_names, q_vals, min_obs = 10) {
    uq <- sort(unique(q_vals))
    samples_u <- unique(sample_names)
    curve_mat <- matrix(NA_real_, nrow = length(samples_u), ncol = length(uq))
    if (length(samples_u) > 0)
        rownames(curve_mat) <- samples_u
    for (i in seq_along(sample_names)) {
        s <- sample_names[i]
        qv <- q_vals[i]
        qi <- match(qv, uq)
        if (is.na(qi))
            next
        if (s %in% rownames(curve_mat))
            curve_mat[s, qi] <- as.numeric(mat[, i])
    }
    # keep samples with at least half of q points present
    good_rows <- which(rowSums(!is.na(curve_mat)) >= max(2, ceiling(ncol(curve_mat)/2)))
    if (length(good_rows) < min_obs)
        return(NULL)
    mat_sub <- curve_mat[good_rows, , drop = FALSE]
    col_means <- apply(mat_sub, 2, function(col) mean(col, na.rm = TRUE))
    for (r in seq_len(nrow(mat_sub))) mat_sub[r, is.na(mat_sub[r, ])] <- col_means[is.na(mat_sub[r,
        ])]
    list(mat_sub = mat_sub, used_samples = rownames(mat_sub))
}
