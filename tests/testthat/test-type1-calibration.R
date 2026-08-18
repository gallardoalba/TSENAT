# =============================================================================
# test-type1-calibration.R — size calibration under H0
# For each method: simulate R datasets under H0 and measure P(p < .05) with a
# Monte Carlo CI. Criterion: the CI must include 0.05.
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

if (!requireNamespace("TSENAT", quietly = TRUE)) {
    # Development use: load R/ into the environment
    for (f in list.files(file.path(script_dir, "..", "..", "R"),
        full.names = TRUE, pattern = "\\.R$")) {
        sys.source(f, envir = .GlobalEnv)
    }
}

run_type1 <- function(fit_fn, n_sim = reps, alpha = 0.05, seed = 1) {
    set.seed(seed)
    pvals <- vapply(seq_len(n_sim), function(i) {
        df <- sim_paired_h0(n_subjects = 10, n_q = 8, rho = 0.5, seed = seed * 1000 + i)
        p <- tryCatch(suppressWarnings(fit_fn(df)), error = function(e) NA_real_)
        as.numeric(p)
    }, FUN.VALUE = numeric(1))
    list(pvals = pvals, rejection_rate = mean(pvals < alpha, na.rm = TRUE),
        n_valid = sum(!is.na(pvals)))
}

# -----------------------------------------------------------------------------
# GEE (geepack) sobre H(q) original con corstr preespecificada
# -----------------------------------------------------------------------------
run_gee <- function(df) {
    df$subject <- factor(df$subject)
    df$group <- factor(df$condition)
    df <- df[order(df$subject, df$group, df$q), ]
    if (!requireNamespace("geepack", quietly = TRUE)) return(NA_real_)
    # Joint Wald on the interaction coefficients (mirrors
    # sait_gee.R .compute_joint_wald_pvalue with small-cluster F correction)
    fit <- geepack::geeglm(entropy ~ q * group, id = subject, data = df,
        family = gaussian(), corstr = "ar1")
    coefs <- stats::coef(fit)
    ia <- which(grepl("^q:group", names(coefs)))
    if (length(ia) == 0) return(NA_real_)
    V <- stats::vcov(fit)[ia, ia, drop = FALSE]
    w <- as.numeric(t(coefs[ia]) %*% solve(V) %*% coefs[ia])
    n_clusters <- length(unique(df$subject))
    f_stat <- w / length(ia)
    stats::pf(f_stat, df1 = length(ia), df2 = n_clusters - length(ia),
        lower.tail = FALSE)
}

# -----------------------------------------------------------------------------
# LMM (nlme) with AR(1) within condition (ADR-001)
# -----------------------------------------------------------------------------
run_lmm <- function(df) {
    if (!requireNamespace("nlme", quietly = TRUE)) return(NA_real_)
    df$subject <- factor(df$subject)
    df$condition <- factor(df$condition)
    df <- df[order(df$subject, df$condition, df$q), ]
    df$time_idx <- sequence(rle(paste(df$subject, df$condition))$lengths)
    fit0 <- nlme::lme(entropy ~ q + condition, random = ~1 | subject,
        correlation = nlme::corAR1(form = ~time_idx | subject/condition),
        data = df, method = "ML")
    fit1 <- nlme::lme(entropy ~ q * condition, random = ~1 | subject,
        correlation = nlme::corAR1(form = ~time_idx | subject/condition),
        data = df, method = "ML")
    anova(fit0, fit1)$`p-value`[2]
}

# -----------------------------------------------------------------------------
# GAMM (mgcv) with AR(1) within condition (ADR-001)
# -----------------------------------------------------------------------------
run_gamm <- function(df) {
    if (!requireNamespace("nlme", quietly = TRUE) || !requireNamespace("splines", quietly = TRUE)) {
        return(NA_real_)
    }
    df$subject <- factor(df$subject)
    df$condition <- factor(df$condition)
    df$group <- df$condition
    df <- df[order(df$subject, df$condition, df$q), ]
    df$obs_seq <- sequence(rle(paste(df$subject, df$condition))$lengths)
    # Corrected AR(1) GAMM: regression splines (ns) + subject random
    # effect + AR(1) within subject×condition (ADR-001). Marginal F-test
    # on the interaction (better calibrated than the LRT with few subjects).
    # The mgcv::gamm formulation with s(q)+s(q,by=group) is singular in lme
    # (null-space collinearity).
    fit1 <- nlme::lme(entropy ~ splines::ns(q, df = 3) * condition,
        random = ~1 | subject,
        correlation = nlme::corAR1(form = ~obs_seq | subject/condition),
        data = df, method = "ML")
    an <- anova(fit1, type = "marginal")
    an$`p-value`[rownames(an) == "splines::ns(q, df = 3):condition"]
}

# -----------------------------------------------------------------------------
# Subject-blocked permutation (nonparametric reference test)
# -----------------------------------------------------------------------------
run_perm <- function(df, n_perm = 199) {
    df$subject <- factor(df$subject)
    df$condition <- factor(df$condition)
    df$q_f <- factor(df$q)
    stat <- function(d) {
        m <- suppressWarnings(lm(rank(entropy) ~ q_f * condition + subject, data = d))
        anova(m)$`F value`[nrow(anova(m)) - 1]
    }
    t_obs <- stat(df)
    subjects <- unique(as.character(df$subject))
    # Paired permutation: swap labels within each (subject, q) pair,
    # which is the exchangeable unit of the paired design.
    t_perm <- replicate(n_perm, {
        d2 <- df
        for (s in subjects) {
            for (qq in unique(df$q)) {
                idx <- which(as.character(d2$subject) == s & d2$q == qq)
                d2$condition[idx] <- sample(d2$condition[idx])
            }
        }
        stat(d2)
    })
    (1 + sum(t_perm >= t_obs, na.rm = TRUE))/(n_perm + 1)
}

results <- list(
    GEE = run_type1(run_gee),
    LMM = run_type1(run_lmm),
    GAMM = run_type1(run_gamm),
    blocked_permutation = run_type1(run_perm)
)

summary_df <- do.call(rbind, lapply(names(results), function(nm) {
    r <- results[[nm]]
    ci <- binom_ci(round(r$rejection_rate * r$n_valid), r$n_valid)
    data.frame(method = nm, n_valid = r$n_valid,
        type1 = round(r$rejection_rate, 4), ci_low = round(ci["lower"], 4),
        ci_high = round(ci["upper"], 4),
        includes_0.05 = contains_value(ci, 0.05))
}))

validation_summary <- summary_df
print(validation_summary)

# testthat integration: assert Monte Carlo criteria when running under testthat
if (requireNamespace("testthat", quietly = TRUE) && testthat::is_testing()) {
    testthat::test_that("test-type1-calibration: Monte Carlo criteria met", {
        testthat::expect_true(all(validation_summary$includes_0.05,
            na.rm = TRUE))
    })
}
