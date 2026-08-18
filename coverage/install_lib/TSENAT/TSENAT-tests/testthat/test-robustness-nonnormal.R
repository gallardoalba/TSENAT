# =============================================================================
# test-robustness-nonnormal.R — type I robustness under non-normality
# (robustness; the validation appendix reports normality p=3.6e-27 and
# heterogeneity p=3.1e-128).
#
# Under H0 with (a) increasing heteroscedasticity in q and (b) 10% outliers,
# the primary GAMM type I (lme_ns + corAR1, mirror of the package) must not
# exclude 0.05; GEE is reported as evidence.
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

n_sims <- as.integer(Sys.getenv("TSENAT_ROBUST_SIMS", "200"))
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

# Paired H0 with increasing heteroscedasticity in q + 10% outliers
sim_h0_robust <- function(seed_i, n_subjects = 10, n_q = 8, rho = 0.5,
    sd_scale = 3, outlier_frac = 0.1, outlier_sd = 3) {
    set.seed(seed_i)
    q_seq <- seq(0.1, 2, length.out = n_q)
    fq <- log1p(q_seq)
    sd_q <- 0.15 + 0.15 * (q_seq - min(q_seq)) * sd_scale  # sd increasing in q
    out <- do.call(rbind, lapply(seq_len(n_subjects), function(s) {
        u <- rnorm(1, sd = 0.7)
        do.call(rbind, lapply(c("A", "B"), function(cond) {
            e <- as.numeric(arima.sim(list(ar = rho), n = n_q, sd = 1)) * sd_q
            y <- u + fq + e
            out_idx <- sample(seq_len(n_q), max(1, round(outlier_frac * n_q)))
            y[out_idx] <- y[out_idx] + rnorm(length(out_idx), sd = outlier_sd)
            data.frame(subject = paste0("S", s), condition = cond, q = q_seq,
                entropy = y)
        }))
    }))
    out
}

run_gamm <- function(df) {
    if (!requireNamespace("nlme", quietly = TRUE) || !requireNamespace("splines",
        quietly = TRUE)) {
        return(NA_real_)
    }
    df$subject <- factor(df$subject)
    df$condition <- factor(df$condition)
    df <- df[order(df$subject, df$condition, df$q), ]
    df$obs_seq <- match(df$q, sort(unique(df$q)))
    fit1 <- nlme::lme(entropy ~ splines::ns(q, df = 3) * condition,
        random = ~1 | subject,
        correlation = nlme::corAR1(form = ~obs_seq | subject/condition),
        data = df, method = "ML")
    an <- anova(fit1, type = "marginal")
    an$`p-value`[rownames(an) == "splines::ns(q, df = 3):condition"]
}

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

p_gamm <- p_gee <- numeric(n_sims)
for (i in seq_len(n_sims)) {
    df <- sim_h0_robust(seed * 8000 + i)
    p_gamm[i] <- tryCatch(suppressWarnings(run_gamm(df)), error = function(e) NA_real_)
    p_gee[i] <- tryCatch(suppressWarnings(run_gee(df)), error = function(e) NA_real_)
}

mk <- function(nm, pv) {
    k <- sum(pv < 0.05, na.rm = TRUE)
    n <- sum(!is.na(pv))
    ci <- binom_ci(k, n)
    data.frame(method = nm, n_valid = n, type1 = round(k/n, 4),
        ci_low = round(ci["lower"], 4), ci_high = round(ci["upper"], 4),
        includes_0.05 = contains_value(ci, 0.05), stringsAsFactors = FALSE)
}

validation_summary <- rbind(
    mk("GAMM heteroscedastic + outliers H0 — CI must not exclude 0.05", p_gamm),
    mk("GEE heteroscedastic + outliers H0 — evidence: anti-conservative (0.085), documented limitation", p_gee))
print(validation_summary)

# testthat integration: assert Monte Carlo criteria when running under testthat
if (requireNamespace("testthat", quietly = TRUE) && testthat::is_testing()) {
    testthat::test_that("test-robustness-nonnormal: Monte Carlo criteria met", {
        testthat::expect_true(isTRUE(validation_summary$includes_0.05[1]))
    })
}
