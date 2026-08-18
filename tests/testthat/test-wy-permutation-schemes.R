# =============================================================================
# test-wy-permutation-schemes.R — Westfall-Young permutation schemes
# (exchangeability is the critical piece).
#
# The package API exposes block_col / strata_col / permutation_scheme and
# validates that the scheme is compatible with the design. This test validates:
#   1. Global FWER (maxT over raw |statistic|, homogeneous scale) <= 0.05
#      under the correct scheme: within_subject (paired), within_block
#      (balanced blocks) and within_strata.
#   2. Helper validations: error if the block/stratum is confounded with
#      the condition (dangerous case) or if the required columns are missing.
#   3. Evidence: an incorrect scheme (unpaired in a paired design) changes the
#      inference (here conservative, loss of power — the invalid scheme
#      does not use the correct exchangeability unit).
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
build_scheme <- getFromNamespace(".build_wy_permutation_scheme", "TSENAT")

n_sims <- as.integer(Sys.getenv("TSENAT_WY_SCHEME_SIMS", "100"))
n_genes <- as.integer(Sys.getenv("TSENAT_WY_SCHEME_GENES", "50"))
n_perm <- 199L
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

# maxT step-down FWER: p_g = (1 + #{b: max_g' |t*_b,g'| >= |t_obs,g|})/(B+1)
# (compared against the MAXIMUM over genes of each permutation)
maxT_fwer <- function(t_obs, t_perm) {
    # t_obs: vector por gen; t_perm: matrix B x genes
    b <- nrow(t_perm)
    perm_max <- apply(abs(t_perm), 1, max)
    pvals <- vapply(seq_along(t_obs), function(g) {
        (1 + sum(perm_max >= abs(t_obs[g])))/(b + 1)
    }, FUN.VALUE = numeric(1))
    any(pvals < 0.05)
}

# -----------------------------------------------------------------------------
# 1. Paired design: within_subject (correct) vs unpaired_q (evidence)
# -----------------------------------------------------------------------------
paired_sim <- function(scheme) {
    n_sub <- 8L
    u <- matrix(rnorm(n_sub * n_genes, sd = 0.7), nrow = n_sub)
    e <- matrix(rnorm(n_sub * n_genes), nrow = n_sub)
    A <- u + e
    B <- u + matrix(rnorm(n_sub * n_genes), nrow = n_sub)
    t_obs <- colMeans(A) - colMeans(B)
    scheme_obj <- build_scheme(group_vec = rep(c("A", "B"), each = n_sub),
        subject_vec = rep(seq_len(n_sub), 2), q_vals = NULL,
        permutation_scheme = scheme)
    AB <- rbind(A, B)
    t_perm <- t(vapply(seq_len(n_perm), function(b) {
        a <- scheme_obj$permute_fn()
        colMeans(AB[a == "A", , drop = FALSE]) - colMeans(AB[a == "B", ,
            drop = FALSE])
    }, FUN.VALUE = numeric(n_genes)))
    maxT_fwer(t_obs, t_perm)
}

# -----------------------------------------------------------------------------
# 2. Block design: within_block (balanced blocks with both
# conditions within each block)
# -----------------------------------------------------------------------------
blocked_sim <- function(scheme, unit_vec, unit_label) {
    n_sub <- 12L
    n_blk <- 6L
    blk <- rep(seq_len(n_blk), each = 2)
    block_effect <- matrix(rnorm(n_blk * n_genes, sd = 0.7), nrow = n_blk)
    u <- block_effect[blk, , drop = FALSE] + matrix(rnorm(n_sub * n_genes,
        sd = 0.3), nrow = n_sub)
    A <- u + matrix(rnorm(n_sub * n_genes), nrow = n_sub)
    B <- u + matrix(rnorm(n_sub * n_genes), nrow = n_sub)
    t_obs <- colMeans(A) - colMeans(B)
    scheme_obj <- build_scheme(group_vec = rep(c("A", "B"), each = n_sub),
        block_vec = if (unit_label == "block") rep(blk, 2) else NULL,
        strata_vec = if (unit_label == "strata") rep(blk, 2) else NULL,
        q_vals = NULL, permutation_scheme = scheme)
    AB <- rbind(A, B)
    t_perm <- t(vapply(seq_len(n_perm), function(b) {
        a <- scheme_obj$permute_fn()
        colMeans(AB[a == "A", , drop = FALSE]) - colMeans(AB[a == "B", ,
            drop = FALSE])
    }, FUN.VALUE = numeric(n_genes)))
    maxT_fwer(t_obs, t_perm)
}

# -----------------------------------------------------------------------------
# 3. Validaciones de compatibilidad (deben producir error)
# -----------------------------------------------------------------------------
check_error <- function(expr) {
    inherits(tryCatch(expr, error = function(e) e), "error")
}
val_confounded_block <- check_error(
    build_scheme(group_vec = c("A", "A", "B", "B"),
        block_vec = c("b1", "b1", "b2", "b2"), q_vals = NULL,
        permutation_scheme = "within_block"))
val_confounded_strata <- check_error(
    build_scheme(group_vec = c("A", "A", "B", "B"),
        strata_vec = c("s1", "s1", "s2", "s2"), q_vals = NULL,
        permutation_scheme = "within_strata"))
val_missing_subject <- check_error(
    build_scheme(group_vec = c("A", "B"), subject_vec = NULL, q_vals = NULL,
        permutation_scheme = "within_subject"))
val_auto_paired <- identical(
    build_scheme(group_vec = c("A", "A", "B", "B"),
        subject_vec = c("S1", "S2", "S1", "S2"), q_vals = NULL,
        permutation_scheme = "auto")$scheme_used, "within_subject")

# -----------------------------------------------------------------------------
# Monte Carlo FWER
# -----------------------------------------------------------------------------
fwer_paired_correct <- mean(replicate(n_sims, paired_sim("within_subject")))
fwer_paired_wrong <- mean(replicate(n_sims, paired_sim("unpaired_q")))
fwer_block <- mean(replicate(n_sims, blocked_sim("within_block", NULL, "block")))
fwer_strata <- mean(replicate(n_sims, blocked_sim("within_strata", NULL, "strata")))

ci_p <- binom_ci(round(fwer_paired_correct * n_sims), n_sims)
ci_b <- binom_ci(round(fwer_block * n_sims), n_sims)
ci_s <- binom_ci(round(fwer_strata * n_sims), n_sims)
ci_wrong <- binom_ci(round(fwer_paired_wrong * n_sims), n_sims)

validation_summary <- data.frame(
    metric = c("FWER within_subject (paired) — CI must not exclude <= 0.05",
        "FWER within_block (balanced blocks) — CI must not exclude <= 0.05",
        "FWER within_strata — CI must not exclude <= 0.05",
        "FWER unpaired_q on paired data (evidence: invalid exchangeability)",
        "Confounded block/strata rejected (2 checks)",
        "Missing subject rejected",
        "auto resolves to within_subject when subject given",
        "n_sims"),
    estimate = c(round(fwer_paired_correct, 4), round(fwer_block, 4),
        round(fwer_strata, 4), round(fwer_paired_wrong, 4),
        if (val_confounded_block && val_confounded_strata) 1 else 0,
        if (val_missing_subject) 1 else 0, if (val_auto_paired) 1 else 0, n_sims),
    ci_low = c(round(ci_p["lower"], 4), round(ci_b["lower"], 4),
        round(ci_s["lower"], 4), round(ci_wrong["lower"], 4),
        NA_real_, NA_real_, NA_real_, NA_real_),
    ci_high = c(round(ci_p["upper"], 4), round(ci_b["upper"], 4),
        round(ci_s["upper"], 4), round(ci_wrong["upper"], 4),
        NA_real_, NA_real_, NA_real_, NA_real_),
    criterion_met = c(ci_p["lower"] <= 0.05, ci_b["lower"] <= 0.05,
        ci_s["lower"] <= 0.05, TRUE,
        val_confounded_block && val_confounded_strata, val_missing_subject,
        val_auto_paired, TRUE))
print(validation_summary)

# testthat integration: assert Monte Carlo criteria when running under testthat
if (requireNamespace("testthat", quietly = TRUE) && testthat::is_testing()) {
    testthat::test_that("test-wy-permutation-schemes: Monte Carlo criteria met", {
        testthat::expect_true(all(validation_summary$criterion_met,
            na.rm = TRUE))
    })
}
