#!/usr/bin/env Rscript
# =============================================================================
# run-validation.R — Monte Carlo validation suite orchestrator.
# Run outside the fast suite:
#   Rscript tests/testthat/run-validation.R
# On Bioconductor the individual test files are skipped via skip_on_bioc();
# this orchestrator runs the whole suite with a fixed seed for manual use.
# =============================================================================

reps <- as.integer(Sys.getenv("TSENAT_VALIDATION_REPS", "300"))
seed <- as.integer(Sys.getenv("TSENAT_VALIDATION_SEED", "42"))
set.seed(seed)

`%||%` <- function(a, b) if (is.null(a)) b else a

.args <- commandArgs(trailingOnly = FALSE)
.file_arg <- .args[startsWith(.args, "--file=")]
script_dir <- if (length(.file_arg)) {
    dirname(normalizePath(sub("^--file=", "", .file_arg[[1]])))
} else {
    getwd()
}
scripts <- file.path(script_dir,
    c("test-type1-calibration.R", "test-type1-rank.R", "test-type1-fpca.R",
        "test-pairing.R", "test-fwer-fdr.R", "test-fwer-westfall-young.R",
        "test-q-grid-invariance.R", "test-estimation-sensitivity.R",
        "test-m-estimator-coverage.R", "test-bootstrap-coverage.R",
        "test-power.R", "test-gee-small-sample.R", "test-ar1-structure.R",
        "test-missing-q.R", "test-filtering.R", "test-effect-sizes.R",
        "test-wy-permutation-schemes.R", "test-postselection-lasso.R",
        "test-lmm-functional.R", "test-type1-art.R",
        "test-qic-auto-selection.R", "test-robustness-nonnormal.R"))

cat(sprintf("== TSENAT statistical validation suite ==\n  reps = %d, seed = %d\n\n", reps, seed))

summaries <- list()
for (s in scripts) {
    if (!file.exists(s)) {
        cat(sprintf("  [SKIP] %s (not found)\n", basename(s)))
        next
    }
    cat(sprintf("  [RUN ] %s\n", basename(s)))
    env <- new.env(parent = globalenv())
    env$TSENAT_VALIDATION_REPS <- reps
    env$TSENAT_VALIDATION_SEED <- seed
    env$TSENAT_VALIDATION_DIR <- script_dir
    out <- try(sys.source(s, envir = env), silent = TRUE)
    if (inherits(out, "try-error")) {
        cat(sprintf("  [FAIL] %s: %s\n", basename(s), conditionMessage(attr(out, "condition"))))
        next
    }
    if (exists("validation_summary", envir = env)) {
        summaries[[basename(s)]] <- get("validation_summary", envir = env)
    }
}

cat("\n== Summary ==\n")
for (nm in names(summaries)) {
    cat(sprintf("-- %s --\n", nm))
    print(summaries[[nm]])
    cat("\n")
}
cat("End of suite. Results should be compared against the acceptance criteria\n")
cat("(Type I≈0.05, FWER≤0.05, FDR≤0.05, coverage≈0.95).\n")
