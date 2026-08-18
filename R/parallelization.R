# Helper functions for parallel processing using BiocParallel and parallel
# package

# Validate nthreads input — ensure numeric, finite, non-negative
.validate_nthreads <- function(nthreads) {
    if (is.null(nthreads)) return(invisible(NULL))
    if (!is.numeric(nthreads) || length(nthreads) != 1) {
        stop("nthreads must be a single numeric value", call. = FALSE)
    }
    if (!is.finite(nthreads) || nthreads < 0) {
        stop("nthreads must be a finite, non-negative number", call. = FALSE)
    }
    invisible(NULL)
}

# Get effective number of threads, respecting environment constraints
.get_effective_nthreads <- function(nthreads = 1) {
    .validate_nthreads(nthreads)

    if (nthreads <= 1) {
        return(1)
    }

    # Check for _R_CHECK_LIMIT_CORES_ environment variable (used by R CMD
    # check)
    core_limit <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
    if (nchar(core_limit) > 0) {
        # Safely coerce environment variable (will return NA if not valid
        # integer)
        core_limit <- tryCatch(as.integer(core_limit), warning = function(w) NA_integer_,
            error = function(e) NA_integer_)
        if (!is.na(core_limit) && is.finite(core_limit) && core_limit > 0) {
            nthreads <- min(nthreads, core_limit)
        }
    }

    # Respect R's standard 'mc.cores' option (set by users or system admins)
    mc_cores <- getOption("mc.cores")
    if (is.numeric(mc_cores) && length(mc_cores) == 1 && is.finite(mc_cores) && mc_cores > 0) {
        nthreads <- min(nthreads, as.integer(mc_cores))
    }

    return(max(1L, as.integer(nthreads)))
}

# ==============================================================================
# Select and Initialize Parallel Backend
# ==============================================================================
# PARAMETERS:
#   nthreads: Number of threads to use (default: 1)
#   seed:     Integer seed for reproducible parallel execution
#
# RETURNS:
#   BiocParallel BPPARAM object
.get_bpparam <- function(nthreads = 1, seed = NULL) {
    # Apply environment variable constraints
    nthreads <- .get_effective_nthreads(nthreads)

    if (nthreads <= 1) {
        return(BiocParallel::SerialParam())
    }

    # RNGseed for reproducible parallel execution.
    # Note: set.seed() in the main R process does NOT propagate to forked
    # workers on Unix (mclapply). The caller should derive a seed from the
    # current RNG state (e.g., sample.int()) and pass it explicitly.
    #
    # BiocParallel >= 1.16.0 (MulticoreParam) / >= 1.14.0 (SnowParam) supports
    # RNGseed. Fall back to no-seed for older versions.
    has_rngseed <- "RNGseed" %in% names(formals(BiocParallel::MulticoreParam))

    if (.Platform$OS.type == "unix") {
        if (has_rngseed && !is.null(seed)) {
            return(BiocParallel::MulticoreParam(workers = nthreads, RNGseed = seed))
        } else {
            return(BiocParallel::MulticoreParam(workers = nthreads))
        }
    } else {
        if (has_rngseed && !is.null(seed)) {
            return(BiocParallel::SnowParam(workers = nthreads, RNGseed = seed))
        } else {
            return(BiocParallel::SnowParam(workers = nthreads))
        }
    }
}

# ==============================================================================
# Apply Function in Parallel Using Best Available Backend
# ==============================================================================
# PARAMETERS:
#   X:         Vector or list to iterate over
#   FUN:       Function to apply
#   nthreads:  Number of threads (default: 1)
#
# NOTE: Always returns a list. Callers that need simplified output should
#   use vapply() / unlist() on the result themselves.
.bplapply <- function(X, FUN, nthreads = 1) {
    .validate_nthreads(nthreads)

    if (nthreads <= 1) {
        return(lapply(X, FUN))
    }

    # Parallel execution
    # Derive seed from current RNG state so that outer set.seed() calls
    # propagate deterministically to workers.
    seed <- sample.int(.Machine$integer.max, 1)
    bpparam <- .get_bpparam(nthreads, seed = seed)

    # BLAS/OpenMP oversubscription guard: the per-gene workers run SMALL
    # linear/mixed-model fits (ART/lm/lme). If every worker uses a
    # multithreaded BLAS, the workers contend for cores and parallel can be
    # SLOWER than serial (measured: ART 42s parallel vs 23s serial). Pin each
    # worker process to a single BLAS/OpenMP thread for the duration of this
    # call; withr::local_envvar() restores (or unsets) the variables
    # automatically on exit, even on error.
    withr::local_envvar(c(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1",
        MKL_NUM_THREADS = "1"))

    return(BiocParallel::bplapply(X, FUN, BPPARAM = bpparam))
}

# ==============================================================================
# Auto-Detect and Validate Number of Threads
# ==============================================================================
# DESCRIPTION:
#   Auto-detects available cores (minus 1) if nthreads is NULL or < 1.
#   Otherwise uses provided value. Always applies environment constraints.
#
# PARAMETERS:
#   nthreads: Integer or NULL; number of threads (default: NULL for auto-detect)
#
# RETURNS:
#   Validated number of threads respecting environment limits
.get_nthreads_auto_detect <- function(nthreads = NULL) {
    if (is.null(nthreads) || (is.numeric(nthreads) && nthreads < 1)) {
        # Auto-detect available cores, leaving one free for system
        # detectCores() can return NA in container/restricted environments;
        # fall back to 1 in that case (na.rm = TRUE handles NA gracefully)
        detected <- parallel::detectCores()
        if (is.na(detected) || !is.finite(detected) || detected < 1) {
            nthreads <- 1L
        } else {
            nthreads <- max(1L, as.integer(detected - 1L))
        }
    } else {
        nthreads <- as.integer(nthreads)
    }

    # Apply environment variable constraints (R CMD check limits, mc.cores)
    return(.get_effective_nthreads(nthreads))
}
