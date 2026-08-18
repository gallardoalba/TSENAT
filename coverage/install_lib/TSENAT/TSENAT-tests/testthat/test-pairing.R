# =============================================================================
# test-pairing.R — pairing-specific validation.
# Y_{s,c,q} = u_s + beta_c*q + eps under H0: beta_A = beta_B.
# Size should be ≈ 0.05 for each method.
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
        ifnotfound = as.integer(Sys.getenv("TSENAT_VALIDATION_REPS", "100"))))
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

sim_paired_slopes_h0 <- function(n_subjects = 20, n_q = 8, sd_subject = 1.5,
    sd_eps = 0.2, seed = NULL) {
    if (!is.null(seed)) set.seed(seed)
    q_seq <- seq(0.1, 2, length.out = n_q)
    do.call(rbind, lapply(seq_len(n_subjects), function(s) {
        u <- rnorm(1, sd = sd_subject)
        do.call(rbind, lapply(c("A", "B"), function(cond) {
            e <- rnorm(n_q, sd = sd_eps)
            data.frame(subject = paste0("S", s), condition = cond, q = q_seq,
                entropy = u + 0.5 * q_seq + e)  # same slope in both conditions
        }))
    }))
}

run_lmm_paired <- function(df) {
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

run_perm_paired <- function(df, n_perm = 199) {
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

p_lmm <- vapply(seq_len(reps), function(i) {
    df <- sim_paired_slopes_h0(seed = i)
    tryCatch(run_lmm_paired(df), error = function(e) NA_real_)
}, FUN.VALUE = numeric(1))
p_perm <- vapply(seq_len(reps), function(i) {
    df <- sim_paired_slopes_h0(seed = i)
    tryCatch(run_perm_paired(df), error = function(e) NA_real_)
}, FUN.VALUE = numeric(1))

make_row <- function(pv, nm) {
    k <- sum(pv < 0.05, na.rm = TRUE)
    n <- sum(!is.na(pv))
    ci <- binom_ci(k, n)
    data.frame(method = nm, n_valid = n, type1 = round(k/n, 4),
        ci_low = round(ci["lower"], 4), ci_high = round(ci["upper"], 4),
        includes_0.05 = contains_value(ci, 0.05))
}

validation_summary <- rbind(
    make_row(p_lmm, "LMM AR(1) within condition"),
    make_row(p_perm, "subject-blocked permutation")
)
print(validation_summary)

# testthat integration: assert Monte Carlo criteria when running under testthat
if (requireNamespace("testthat", quietly = TRUE) && testthat::is_testing()) {
    testthat::test_that("test-pairing: Monte Carlo criteria met", {
        testthat::expect_true(all(validation_summary$includes_0.05,
            na.rm = TRUE))
    })
}
