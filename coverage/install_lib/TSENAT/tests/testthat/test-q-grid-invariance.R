# =============================================================================
# test-q-grid-invariance.R — p-value stability under q discretization.
# If q is treated as a numeric covariate, refining the grid must not
# qualitatively change the p-value ranking. Evaluated under the alternative
# (different slope per condition) with the GAMM lme_ns AR(1) (the package
# implementation).
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

seed <- get0("TSENAT_VALIDATION_SEED", envir = environment(),
    ifnotfound = get0("TSENAT_VALIDATION_SEED", envir = .GlobalEnv,
        ifnotfound = as.integer(Sys.getenv("TSENAT_VALIDATION_SEED", "42"))))
n_genes <- as.integer(Sys.getenv("TSENAT_INVARIANCE_GENES", "60"))
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
sim_paired_h1 <- function(n_subjects = 10, n_q = 8, effect = 0.8, rho = 0.5,
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

#' GAMM lme_ns AR(1) with marginal F (same implementation as the package,
#' including the guard against duplicated q per block)
run_gamm_lme_ns <- function(df) {
    df$subject <- factor(df$subject)
    df$condition <- factor(df$condition)
    df <- df[order(df$subject, df$condition, df$q), ]
    df$obs_seq <- sequence(rle(paste(df$subject, df$condition))$lengths)
    q_dup <- tapply(df$q, paste(df$subject, df$condition), anyDuplicated)
    if (any(q_dup > 0, na.rm = TRUE)) {
        stop("Duplicated q within subject x condition block")
    }
    fit <- nlme::lme(entropy ~ splines::ns(q, df = 3) * condition,
        random = ~1 | subject,
        correlation = nlme::corAR1(form = ~obs_seq | subject/condition),
        data = df, method = "ML")
    an <- anova(fit, type = "marginal")
    an$`p-value`[rownames(an) == "splines::ns(q, df = 3):condition"]
}

p_grid <- function(n_q, eff) {
    vapply(seq_len(n_genes), function(g) {
        df <- sim_paired_h1(n_q = n_q, effect = eff[g], seed = seed * 1000 + g)
        tryCatch(run_gamm_lme_ns(df), error = function(e) NA_real_)
    }, FUN.VALUE = numeric(1))
}

nq_A <- 8L
nq_B <- 16L
# Per-gene varying effects to cover a range of p-values so that the
# ranking correlation is informative (all p-values ~0 would give noise)
set.seed(seed)
eff_g <- stats::runif(n_genes, 0.1, 1.2)
p_A <- p_grid(nq_A, eff_g)
p_B <- p_grid(nq_B, eff_g)
valid <- !is.na(p_A) & !is.na(p_B)
spearman <- suppressWarnings(cor(p_A[valid], p_B[valid], method = "spearman"))
sig_A <- p.adjust(p_A, "BH") < 0.05
sig_B <- p.adjust(p_B, "BH") < 0.05
inter <- sum(sig_A & sig_B, na.rm = TRUE)
union_ <- sum(sig_A | sig_B, na.rm = TRUE)
jaccard <- if (union_ > 0) inter/union_ else NA_real_

# Row reordering invariance (10 genes, H1). The method
# sorts internally by subject/condition/q, so the p-value must be
# identical.
order_diff <- tryCatch({
    p_ord <- sapply(1:10, function(g) {
        df <- sim_paired_h1(n_q = 8, effect = 0.8, seed = seed * 100 + g)
        p1 <- run_gamm_lme_ns(df)
        p2 <- run_gamm_lme_ns(df[sample(nrow(df)), ])
        abs(p1 - p2)
    })
    max(p_ord, na.rm = TRUE)
}, error = function(e) NA_real_)

# Entropy scale invariance (10 genes, H1).
# H' = c*H must not change the interaction test.
scale_diff <- tryCatch({
    p_scale <- sapply(1:10, function(g) {
        df <- sim_paired_h1(n_q = 8, effect = 0.8, seed = seed * 200 + g)
        p1 <- run_gamm_lme_ns(df)
        df2 <- df
        df2$entropy <- 2 * df2$entropy
        p2 <- run_gamm_lme_ns(df2)
        df3 <- df
        df3$entropy <- df3$entropy/log(2)
        p3 <- run_gamm_lme_ns(df3)
        max(abs(p1 - p2), abs(p1 - p3))
    })
    max(p_scale, na.rm = TRUE)
}, error = function(e) NA_real_)

# Duplicated q grid (functional pseudoreplication) must be rejected:
# corAR1 requires a unique covariate per group, so the method must fail.
p_dup <- sapply(1:10, function(g) {
    df <- sim_paired_h1(n_q = 8, effect = 0.8, seed = seed * 300 + g)
    df <- rbind(df, df)  # artificially duplicate the grid
    tryCatch(run_gamm_lme_ns(df), error = function(e) NA_real_)
})
n_dup_valid <- sum(!is.na(p_dup))

validation_summary <- data.frame(
    metric = c("Spearman rank cor (p-values, grid 8q vs 16q)",
        "Jaccard (BH<0.05 significant genes)",
        "Max |p_diff| reordering rows",
        "Max |p_diff| entropy scale",
        "n_valid duplicated q grid (expect 0)",
        "n_genes", "n_valid"),
    estimate = c(round(spearman, 4), round(jaccard, 4), round(order_diff, 8),
        round(scale_diff, 8), n_dup_valid, n_genes, sum(valid)),
    ci_low = NA_real_, ci_high = NA_real_,
    criterion_met = c(!is.na(spearman) && spearman > 0.8,
        !is.na(jaccard) && jaccard > 0.4,
        !is.na(order_diff) && order_diff < 1e-06,
        !is.na(scale_diff) && scale_diff < 1e-06,
        n_dup_valid == 0, TRUE, TRUE))
print(validation_summary)

# testthat integration: assert Monte Carlo criteria when running under testthat
if (requireNamespace("testthat", quietly = TRUE) && testthat::is_testing()) {
    testthat::test_that("test-q-grid-invariance: Monte Carlo criteria met", {
        testthat::expect_true(all(validation_summary$criterion_met,
            na.rm = TRUE))
    })
}
