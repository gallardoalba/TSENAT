# =============================================================================
# test-estimation-sensitivity.R — estimation validations:
# analytical behavior of the Tsallis entropy at q=0, q=1, q=2,
# normalization, and pseudocount sensitivity.
# Uses the package's real implementation (pkgload::load_all) with fallback to
# local formulas if the package is not available.
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

# Load the package from source (to validate ITS implementation);
# in the BBS environment (copied tests, no DESCRIPTION) use the installed package.
# Under testthat the package is already loaded; load_all() here would unload
# the installed namespace and break the whole run (covr/R CMD check).
pkg_loaded <- if (requireNamespace("testthat", quietly = TRUE) && testthat::is_testing()) {
    requireNamespace("TSENAT", quietly = TRUE)
} else {
    tryCatch({
        pkgload::load_all(file.path(script_dir, "..", ".."), quiet = TRUE)
        TRUE
    }, error = function(e) requireNamespace("TSENAT", quietly = TRUE))
}

entropy_fn <- if (pkg_loaded && exists(".entropy_core", where = asNamespace("TSENAT"))) {
    get(".entropy_core", envir = asNamespace("TSENAT"))
} else {
    # Local fallback: standard Tsallis formula
    function(proportions, q = 1, norm = FALSE, log_base = exp(1), q_tol = 1e-06) {
        p <- proportions[proportions > 0]
        if (abs(q - 1) < q_tol) return(-sum(p * log(p)))
        if (abs(q) < q_tol) return(length(p) - 1)
        (1 - sum(p^q))/(q - 1)
    }
}

tol <- 1e-06

# Test distributions (with and without zeros)
p1 <- c(0.4, 0.3, 0.2, 0.1)
p2 <- c(0.5, 0.5, 0, 0)

# 1) q=0 → richness n_present - 1
q0_ok <- abs(entropy_fn(p1, q = 0) - 3) < tol && abs(entropy_fn(p2, q = 0) - 1) <
    tol

# 2) q=1 → Shannon
sh1 <- -sum(p1 * log(p1))
sh2 <- -sum(p2[p2 > 0] * log(p2[p2 > 0]))
q1_ok <- abs(entropy_fn(p1, q = 1) - sh1) < tol && abs(entropy_fn(p2, q = 1) -
    sh2) < tol

# 3) q=2 → Gini-Simpson (1 - sum p^2)
q2_ok <- abs(entropy_fn(p1, q = 2) - (1 - sum(p1^2))) < tol && abs(entropy_fn(p2,
    q = 2) - (1 - sum(p2^2))) < tol

# 4) Normalization: normalized entropy == entropy / maximum (q=0.5)
q_half <- 0.5
max_half <- (1 - 4^(1 - q_half))/(q_half - 1)
norm_ok <- abs(entropy_fn(p1, q = q_half, norm = TRUE) - entropy_fn(p1, q = q_half)/max_half) <
    tol

# 5) Pseudocount sensitivity (evidence): relative entropy change for
#    counts with zeros when adding a pseudocount
counts <- c(100, 80, 50, 30, 0, 0, 0, 0)
pcs <- c(0, 0.1, 0.5, 1, 2)
ent_pc <- sapply(pcs, function(pc) {
    p <- (counts + pc)/sum(counts + pc)
    entropy_fn(p, q = 1)
})
rel_change <- abs(ent_pc - ent_pc[1])/ent_pc[1]
# Criterion: typical pseudocounts (<= 0.5, usual range in quantification)
# must change the entropy < 5%
pc_ok <- max(rel_change[pcs <= 0.5]) < 0.05

# 6) Tsallis divergence with zero support: the package implementation
#    must return finite values for distributions with zeros
div_fn <- if (pkg_loaded && exists(".compute_tsallis_divergence", where = asNamespace("TSENAT"))) {
    get(".compute_tsallis_divergence", envir = asNamespace("TSENAT"))
} else {
    function(p, r, q) NA_real_
}
d_zero <- c(1, 0)
d_half <- c(0.5, 0.5)
d05 <- div_fn(d_zero, d_half, 0.5)
d1 <- div_fn(d_zero, d_half, 1)
d2 <- div_fn(d_zero, d_half, 2)
d_self <- div_fn(d_zero, d_zero, 1)
div_finite_ok <- all(is.finite(c(d05, d1, d2))) && !is.na(d_self) &&
    abs(d_self) < tol

validation_summary <- data.frame(
    metric = c("q=0 equals richness n_present-1", "q=1 equals Shannon entropy",
        "q=2 equals Gini-Simpson", "Normalization = entropy / max",
        "Pseudocount (<=0.5) relative change", "Divergence finite with zero support",
        "Package loaded from source"),
    estimate = c(as.numeric(q0_ok), as.numeric(q1_ok), as.numeric(q2_ok),
        as.numeric(norm_ok), round(max(rel_change[pcs <= 0.5]), 4),
        as.numeric(div_finite_ok), as.numeric(pkg_loaded)),
    ci_low = NA_real_, ci_high = NA_real_,
    criterion_met = c(q0_ok, q1_ok, q2_ok, norm_ok, pc_ok, div_finite_ok,
        pkg_loaded))
print(validation_summary)

# testthat integration: assert Monte Carlo criteria when running under testthat
if (requireNamespace("testthat", quietly = TRUE) && testthat::is_testing()) {
    testthat::test_that("test-estimation-sensitivity: Monte Carlo criteria met", {
        testthat::expect_true(all(validation_summary$criterion_met,
            na.rm = TRUE))
    })
}
