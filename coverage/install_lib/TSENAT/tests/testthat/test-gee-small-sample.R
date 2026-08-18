# =============================================================================
# test-gee-small-sample.R — small-sample GEE validation:
# compares the joint Wald p-value under chi-squared vs. the
# F(k, n_clusters - k) correction for different numbers of clusters, under H0.
# Criterion: the type I CI must include 0.05.
# =============================================================================
`%||%` <- function(a, b) if (is.null(a)) b else a
script_dir <- get0("TSENAT_VALIDATION_DIR", envir = environment(),
    ifnotfound = get0("TSENAT_VALIDATION_DIR", envir = .GlobalEnv,
        ifnotfound = {
            .args <- commandArgs(trailingOnly = FALSE)
            .fa <- .args[startsWith(.args, "--file=")]
            if (length(.fa)) {
                .sd <- dirname(sub("^--file=", "", .fa[[1]]))
                if (file.exists(file.path(.sd, "helpers.R"))) .sd else getwd()
            } else getwd()
        }))
sys.source(file.path(script_dir, "helpers.R"), envir = environment())

reps <- as.integer(Sys.getenv("TSENAT_GEE_REPS", "150"))
seed <- get0("TSENAT_VALIDATION_SEED", envir = environment(),
    ifnotfound = get0("TSENAT_VALIDATION_SEED", envir = .GlobalEnv,
        ifnotfound = as.integer(Sys.getenv("TSENAT_VALIDATION_SEED", "42"))))
set.seed(seed)

# testthat integration: the Monte Carlo suite is heavy and must NOT run under
# R CMD check (it made CI exceed its 30-minute no-output timeout). Skipped by
# default under testthat (Bioconductor builds always; elsewhere unless
# explicitly requested). Run the full suite locally with:
#   Rscript tests/testthat/run-validation.R
# or opt in under testthat with env TSENAT_RUN_VALIDATION=true
# (reps/seed via TSENAT_VALIDATION_REPS / TSENAT_VALIDATION_SEED).
if (requireNamespace("testthat", quietly = TRUE) && testthat::is_testing()) {
    testthat::skip_on_bioc()
    if (!isTRUE(as.logical(Sys.getenv("TSENAT_RUN_VALIDATION", "false")))) {
        testthat::skip("Monte Carlo validation suite skipped under testthat; run tests/testthat/run-validation.R or set TSENAT_RUN_VALIDATION=true")
    }
}

#' Fit the interaction GEE and return the joint Wald with both
#' reference distributions: chi-squared and small-cluster F with
#' df2 = n_clusters - p (p = estimated coefficients; the package's
#' mitigation, sait_gee.R)
run_gee_both <- function(df) {
    df$subject <- factor(df$subject)
    df$group <- factor(df$condition)
    df <- df[order(df$subject, df$group, df$q), ]
    fit <- geepack::geeglm(entropy ~ q * group, id = subject, data = df,
        family = gaussian(), corstr = "ar1")
    coefs <- stats::coef(fit)
    ia <- which(grepl("^q:group", names(coefs)))
    if (length(ia) == 0) return(c(p_chisq = NA_real_, p_f = NA_real_))
    V <- stats::vcov(fit)[ia, ia, drop = FALSE]
    w <- as.numeric(t(coefs[ia]) %*% solve(V) %*% coefs[ia])
    n_clusters <- length(unique(df$subject))
    p_model <- length(coefs)
    c(p_chisq = stats::pchisq(w, df = length(ia), lower.tail = FALSE),
        p_f = stats::pf(w/length(ia), df1 = length(ia),
            df2 = max(1, n_clusters - p_model), lower.tail = FALSE))
}

cells <- c(6L, 10L, 20L)
rows <- do.call(rbind, lapply(cells, function(ncl) {
    out <- lapply(seq_len(reps), function(i) {
        df <- sim_paired_h0(n_subjects = ncl, n_q = 8, rho = 0.5,
            seed = seed * 1000 + ncl * reps + i)
        tryCatch(run_gee_both(df), error = function(e) c(NA_real_, NA_real_))
    })
    mat <- do.call(rbind, out)
    do.call(rbind, lapply(c("p_chisq", "p_f"), function(col) {
        pv <- mat[, col]
        k <- sum(pv < 0.05, na.rm = TRUE)
        n <- sum(!is.na(pv))
        ci <- binom_ci(k, n)
        is_f <- col == "p_f"
        data.frame(statistic = if (col == "p_chisq") "Wald chi2 (evidence)" else "Wald F small-cluster",
            n_clusters = ncl, n_valid = n, type1 = round(k/n, 4),
            ci_low = round(ci["lower"], 4), ci_high = round(ci["upper"], 4),
            includes_0.05 = contains_value(ci, 0.05),
            # The nominal criterion is kept for the F version (the one chosen by
            # the package): the real finding is that at n=6 clusters it is still
            # anti-conservative.
            criterion_met = if (is_f) contains_value(ci, 0.05) else TRUE)
    }))
}))

validation_summary <- rows
print(validation_summary)

# testthat integration: assert Monte Carlo criteria when running under testthat
if (requireNamespace("testthat", quietly = TRUE) && testthat::is_testing()) {
    testthat::test_that("test-gee-small-sample: Monte Carlo criteria met", {
        testthat::expect_true(all(validation_summary$criterion_met,
            na.rm = TRUE))
    })
}
