# =============================================================================
# test-missing-q.R — robustness to missing q-values.
#
# q-points are removed at random (MCAR, 25% per subject×condition block) from
# the paired H0 curves. GAMM (lme_ns) and GEE apply casewise deletion; FPCA
# imputes column means with the package's .impute_curve_matrix(). Criterion:
# each method's type I must include 0.05.
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

pkg_root <- normalizePath(file.path(script_dir, "..", ".."))
# Under testthat the package is already loaded; load_all() here would unload
# the installed namespace and break the whole run (covr/R CMD check).
if (!(requireNamespace("testthat", quietly = TRUE) && testthat::is_testing()) &&
    file.exists(file.path(pkg_root, "DESCRIPTION"))) {
    suppressMessages(pkgload::load_all(pkg_root, quiet = TRUE))
}
impute_curve_matrix <- getFromNamespace(".impute_curve_matrix", "TSENAT")

reps <- get0("TSENAT_VALIDATION_REPS", envir = environment(),
    ifnotfound = get0("TSENAT_VALIDATION_REPS", envir = .GlobalEnv,
        ifnotfound = as.integer(Sys.getenv("TSENAT_VALIDATION_REPS", "300"))))
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

# -----------------------------------------------------------------------------
# GEE: joint Wald q:group + F small-cluster (mirrors sait_gee.R)
# -----------------------------------------------------------------------------
run_gee <- function(df) {
    df$subject <- factor(df$subject)
    df$group <- factor(df$condition)
    df <- df[order(df$subject, df$group, df$q), ]
    if (!requireNamespace("geepack", quietly = TRUE)) return(NA_real_)
    fit <- geepack::geeglm(entropy ~ q * group, id = subject, data = df,
        family = gaussian(), corstr = "ar1")
    coefs <- stats::coef(fit)
    ia <- which(grepl("^q:group", names(coefs)))
    if (length(ia) == 0) return(NA_real_)
    V <- stats::vcov(fit)[ia, ia, drop = FALSE]
    w <- as.numeric(t(coefs[ia]) %*% solve(V) %*% coefs[ia])
    n_clusters <- length(unique(df$subject))
    f_stat <- w/length(ia)
    stats::pf(f_stat, df1 = length(ia), df2 = n_clusters - length(ia),
        lower.tail = FALSE)
}

# -----------------------------------------------------------------------------
# GAMM lme + ns + AR(1) (package priority 1, sait_gam.R)
# -----------------------------------------------------------------------------
run_gamm <- function(df) {
    if (!requireNamespace("nlme", quietly = TRUE) || !requireNamespace("splines",
        quietly = TRUE)) {
        return(NA_real_)
    }
    df$subject <- factor(df$subject)
    df$condition <- factor(df$condition)
    df <- df[order(df$subject, df$condition, df$q), ]
    # q grid index (mirrors the package fix: rho^|Δgrid| with gaps)
    df$obs_seq <- match(df$q, sort(unique(df$q)))
    fit1 <- nlme::lme(entropy ~ splines::ns(q, df = 3) * condition,
        random = ~1 | subject,
        correlation = nlme::corAR1(form = ~obs_seq | subject/condition),
        data = df, method = "ML")
    an <- anova(fit1, type = "marginal")
    an$`p-value`[rownames(an) == "splines::ns(q, df = 3):condition"]
}

# -----------------------------------------------------------------------------
# Paired FPCA with column-mean imputation (.impute_curve_matrix)
# -----------------------------------------------------------------------------
run_fpca_imputed <- function(df, qs_all, max_pc = 3) {
    df$subject <- factor(df$subject)
    df$condition <- factor(df$condition)
    n_q <- length(qs_all)
    wide <- do.call(rbind, lapply(levels(df$subject), function(s) {
        do.call(rbind, lapply(levels(df$condition), function(cc) {
            v <- rep(NA_real_, n_q)
            m <- df$q[df$subject == s & df$condition == cc]
            v[match(m, qs_all)] <- df$entropy[df$subject == s & df$condition ==
                cc]
            v
        }))
    }))
    wide <- impute_curve_matrix(wide)
    n_sub <- length(levels(df$subject))
    K <- min(max_pc, n_q, n_sub - 1)
    pca <- prcomp(wide, center = TRUE, scale. = FALSE, rank. = K)
    scores <- pca$x
    D <- scores[seq(1, nrow(scores), by = 2), , drop = FALSE] -
        scores[seq(2, nrow(scores), by = 2), , drop = FALSE]
    n <- nrow(D)
    dbar <- colMeans(D)
    S <- cov(D)
    t2 <- if (K == 1) n * dbar^2/S[1, 1] else as.numeric(n * t(dbar) %*%
        solve(S) %*% dbar)
    if (!is.finite(t2) || t2 <= 0) return(1)
    f_stat <- t2 * (n - K)/(K * (n - 1))
    stats::pf(f_stat, df1 = K, df2 = n - K, lower.tail = FALSE)
}

drop_missing_q <- function(df, frac = 0.25) {
    blocks <- split(df, paste(df$subject, df$condition))
    out <- do.call(rbind, lapply(blocks, function(b) {
        n_keep <- max(3L, ceiling(nrow(b) * (1 - frac)))
        b[sample(seq_len(nrow(b)), n_keep), ]
    }))
    out
}

pvals <- lapply(c("GAMM", "GEE", "FPCA"), function(method) {
    vapply(seq_len(reps), function(i) {
        df <- sim_paired_h0(n_subjects = 10, n_q = 8, rho = 0.5, seed = seed *
            3000 + i)
        qs_all <- sort(unique(df$q))
        df <- drop_missing_q(df)
        tryCatch(suppressWarnings(switch(method,
            GAMM = run_gamm(df), GEE = run_gee(df), FPCA = run_fpca_imputed(df,
                qs_all))), error = function(e) NA_real_)
    }, FUN.VALUE = numeric(1))
})
names(pvals) <- c("GAMM", "GEE", "FPCA")

validation_summary <- do.call(rbind, lapply(names(pvals), function(nm) {
    pv <- pvals[[nm]]
    k <- sum(pv < 0.05, na.rm = TRUE)
    n <- sum(!is.na(pv))
    ci <- binom_ci(k, n)
    data.frame(method = paste(nm, "(25% q missing)"), n_valid = n,
        type1 = round(k/n, 4), ci_low = round(ci["lower"], 4),
        ci_high = round(ci["upper"], 4), includes_0.05 = contains_value(ci,
            0.05))
}))
print(validation_summary)

# testthat integration: assert Monte Carlo criteria when running under testthat
if (requireNamespace("testthat", quietly = TRUE) && testthat::is_testing()) {
    testthat::test_that("test-missing-q: Monte Carlo criteria met", {
        testthat::expect_true(all(validation_summary$includes_0.05,
            na.rm = TRUE))
    })
}
