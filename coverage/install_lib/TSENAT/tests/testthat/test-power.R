# =============================================================================
# test-power.R — power curves under the alternative.
# Evaluates the power of GAMM (lme_ns AR(1)), LMM AR(1) and GEE (Wald F)
# for different effects and sample sizes.
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

reps <- as.integer(Sys.getenv("TSENAT_POWER_REPS", "50"))
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

#' Simulate paired curves under the alternative: different slope per condition
sim_paired_h1 <- function(n_subjects = 10, n_q = 8, effect = 0.5, rho = 0.5,
    sd_subject = 1, sd_eps = 0.25, seed = NULL) {
    if (!is.null(seed)) set.seed(seed)
    q_seq <- seq(0.1, 2, length.out = n_q)
    do.call(rbind, lapply(seq_len(n_subjects), function(s) {
        u <- rnorm(1, sd = sd_subject)
        do.call(rbind, lapply(c("A", "B"), function(cond) {
            e <- as.numeric(arima.sim(list(ar = rho), n = n_q, sd = sd_eps))
            slope <- if (cond == "B") 0.3 + effect else 0.3
            data.frame(subject = paste0("S", s), condition = cond, q = q_seq,
                entropy = u + slope * q_seq + e)
        }))
    }))
}

# GEE joint Wald with small-cluster F correction (package implementation)
run_gee <- function(df) {
    df$subject <- factor(df$subject)
    df$group <- factor(df$condition)
    df <- df[order(df$subject, df$group, df$q), ]
    fit <- geepack::geeglm(entropy ~ q * group, id = subject, data = df,
        family = gaussian(), corstr = "ar1")
    coefs <- stats::coef(fit)
    ia <- which(grepl("^q:group", names(coefs)))
    if (length(ia) == 0) return(NA_real_)
    V <- stats::vcov(fit)[ia, ia, drop = FALSE]
    w <- as.numeric(t(coefs[ia]) %*% solve(V) %*% coefs[ia])
    n_clusters <- length(unique(df$subject))
    stats::pf(w/length(ia), df1 = length(ia), df2 = n_clusters - length(ia),
        lower.tail = FALSE)
}

# LMM AR(1) within condition (nested LRT)
run_lmm <- function(df) {
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

# GAMM lme_ns AR(1) with marginal F (package implementation)
run_gamm <- function(df) {
    df$subject <- factor(df$subject)
    df$condition <- factor(df$condition)
    df <- df[order(df$subject, df$condition, df$q), ]
    df$obs_seq <- sequence(rle(paste(df$subject, df$condition))$lengths)
    fit <- nlme::lme(entropy ~ splines::ns(q, df = 3) * condition,
        random = ~1 | subject,
        correlation = nlme::corAR1(form = ~obs_seq | subject/condition),
        data = df, method = "ML")
    an <- anova(fit, type = "marginal")
    an$`p-value`[rownames(an) == "splines::ns(q, df = 3):condition"]
}

methods <- list(GAMM = run_gamm, LMM = run_lmm, GEE = run_gee)
cells <- expand.grid(n_subjects = c(8L, 12L, 20L), effect = c(0.2, 0.5, 1))

rows <- lapply(names(methods), function(mn) {
    do.call(rbind, lapply(seq_len(nrow(cells)), function(ci) {
        n_sub <- cells$n_subjects[ci]
        eff <- cells$effect[ci]
        pv <- vapply(seq_len(reps), function(i) {
            df <- sim_paired_h1(n_subjects = n_sub, n_q = 8, effect = eff,
                seed = seed * 1000 + ci * reps + i)
            tryCatch(methods[[mn]](df), error = function(e) NA_real_)
        }, FUN.VALUE = numeric(1))
        data.frame(method = mn, n_subjects = n_sub, effect = eff,
            power = mean(pv < 0.05, na.rm = TRUE),
            n_valid = sum(!is.na(pv)))
    }))
})

power_table <- do.call(rbind, rows)

# Criterion: power monotonically increasing with effect (within each method and n).
# Not strict: power may saturate at 1.0 with large effects.
monotonic_ok <- vapply(split(power_table, interaction(power_table$method,
    power_table$n_subjects)), function(x) {
    x <- x[order(x$effect), ]
    if (nrow(x) < 3 || any(is.na(x$power))) return(FALSE)
    x$power[3] >= x$power[2] && x$power[2] >= x$power[1] && x$power[3] >
        x$power[1]
}, FUN.VALUE = logical(1))

validation_summary <- power_table
validation_summary$criterion_met <- unname(monotonic_ok[interaction(power_table$method,
    power_table$n_subjects)])
print(validation_summary)

# testthat integration: assert Monte Carlo criteria when running under testthat
if (requireNamespace("testthat", quietly = TRUE) && testthat::is_testing()) {
    testthat::test_that("test-power: Monte Carlo criteria met", {
        testthat::expect_true(all(validation_summary$criterion_met,
            na.rm = TRUE))
    })
}
