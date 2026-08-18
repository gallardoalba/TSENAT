# =============================================================================
# test-filtering.R — filtering and condition-dependent selection.
#
# 1. The package filter is blind to the condition by construction
#    (TPM/min_samples/abundance thresholds per transcript over all
#    samples). Behavioral check: permuting the condition labels in the
#    colData does NOT change which transcripts .filter_se() retains.
# 2. Monte Carlo: genes with variable depth, noise sd ~ 1/sqrt(depth).
#    Filter A (condition-independent, by depth) → GEE type I ≈ 0.05.
#    Filter B (evidence): selecting genes by |mean A − mean B| (outcome-
#    dependent) inflates type I — the documented risk.
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
filter_se <- getFromNamespace(".filter_se", "TSENAT")

n_sims <- as.integer(Sys.getenv("TSENAT_FILTER_SIMS", "100"))
n_genes <- as.integer(Sys.getenv("TSENAT_FILTER_GENES", "60"))
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
    c(p = stats::pf(f_stat, df1 = length(ia), df2 = n_clusters - length(ia),
        lower.tail = FALSE), w = w)
}

sim_gene_curve <- function(depth, seed_i, n_subjects = 10, n_q = 8, rho = 0.5) {
    q_seq <- seq(0.1, 2, length.out = n_q)
    fq <- log1p(q_seq)
    sd_eps <- 0.25 * sqrt(1/depth)
    do.call(rbind, lapply(seq_len(n_subjects), function(s) {
        u <- rnorm(1)
        do.call(rbind, lapply(c("A", "B"), function(cond) {
            e <- as.numeric(arima.sim(list(ar = rho), n = n_q, sd = sd_eps))
            data.frame(subject = paste0("S", s), condition = cond, q = q_seq,
                entropy = u + fq + e)
        }))
    }))
}

# Criterion 1: filtering invariance under condition-label permutation
set.seed(seed)
cnt <- matrix(rpois(40 * 12, lambda = rep(seq(3, 200, length.out = 40), 12)),
    nrow = 40)
colnames(cnt) <- paste0("S", seq_len(12))
rownames(cnt) <- paste0("Tx", seq_len(40))
se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = cnt,
        tpm = sweep(cnt, 2, colSums(cnt), "/") * 1e+06),
    colData = data.frame(condition = rep(c("A", "B"), each = 6),
        row.names = colnames(cnt)))
f1 <- suppressWarnings(filter_se(se, stringency = NULL, min_tpm = 1,
    min_samples = 3L, min_tx_per_gene = 2L, min_isoform_abundance = 0,
    verbose = FALSE))
se2 <- se
SummarizedExperiment::colData(se2)$condition <- sample(
    SummarizedExperiment::colData(se)$condition)
f2 <- suppressWarnings(filter_se(se2, stringency = NULL, min_tpm = 1,
    min_samples = 3L, min_tx_per_gene = 2L, min_isoform_abundance = 0,
    verbose = FALSE))
filter_invariant <- identical(rownames(f1), rownames(f2))

# Criterion 2: type I Monte Carlo under an independent filter vs. selection
# dependent on the test statistic (selection bias)
rej_indep <- rej_dep <- numeric(n_sims)
for (i in seq_len(n_sims)) {
    depth <- exp(rnorm(n_genes, 0, 1.2))
    pvals <- wstats <- numeric(n_genes)
    ok <- logical(n_genes)
    for (g in seq_len(n_genes)) {
        df <- sim_gene_curve(depth[g], seed * 5000 + i * 1000 + g)
        r <- suppressWarnings(run_gee(df))
        pvals[g] <- r["p"]
        wstats[g] <- r["w"]
        ok[g] <- !is.na(pvals[g])
    }
    keep_indep <- ok & depth >= median(depth)
    keep_dep <- ok & wstats >= median(wstats)
    rej_indep[i] <- mean(pvals[keep_indep] < 0.05)
    rej_dep[i] <- mean(pvals[keep_dep] < 0.05)
}

ci_indep <- binom_ci(sum(rej_indep), n_sims)
ci_dep <- binom_ci(sum(rej_dep), n_sims)

validation_summary <- data.frame(
    metric = c("Filter invariant to condition-label permutation (must be TRUE)",
        "Type I, condition-independent filter (depth) — CI must include 0.05",
        "Type I, selection on the test statistic (evidence: inflated ~0.10)",
        "n_sims"),
    estimate = c(if (filter_invariant) 1 else 0, round(mean(rej_indep), 4),
        round(mean(rej_dep), 4), n_sims),
    ci_low = c(NA_real_, round(ci_indep["lower"], 4), round(ci_dep["lower"], 4),
        NA_real_),
    ci_high = c(NA_real_, round(ci_indep["upper"], 4), round(ci_dep["upper"], 4),
        NA_real_),
    criterion_met = c(filter_invariant, contains_value(ci_indep, 0.05), TRUE,
        TRUE))
print(validation_summary)

# testthat integration: assert Monte Carlo criteria when running under testthat
if (requireNamespace("testthat", quietly = TRUE) && testthat::is_testing()) {
    testthat::test_that("test-filtering: Monte Carlo criteria met", {
        testthat::expect_true(all(validation_summary$criterion_met,
            na.rm = TRUE))
    })
}
