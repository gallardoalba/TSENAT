# Helper functions for parallel processing using BiocParallel and parallel
# package

# Get effective number of threads, respecting environment constraints
.get_effective_nthreads <- function(nthreads = 1) {
    if (nthreads <= 1) {
        return(1)
    }

    # Check for _R_CHECK_LIMIT_CORES_ environment variable (used by R CMD check)
    core_limit <- Sys.getenv("_R_CHECK_LIMIT_CORES_", NA)
    if (!is.na(core_limit)) {
        core_limit <- as.integer(core_limit)
        if (is.finite(core_limit) && core_limit > 0) {
            nthreads <- min(nthreads, core_limit)
        }
    }

    return(max(1, nthreads))
}

# Select and initialize parallel backend @param nthreads Number of threads to
# use (default: 1) @return BiocParallel BPPARAM object
.get_bpparam <- function(nthreads = 1) {
    # Apply environment variable constraints
    nthreads <- .get_effective_nthreads(nthreads)

    if (nthreads <= 1) {
        return(BiocParallel::SerialParam())
    }

    # On Unix-like systems, use MulticoreParam; on Windows use SnowParam
    if (.Platform$OS.type == "unix") {
        return(BiocParallel::MulticoreParam(workers = nthreads))
    } else {
        return(BiocParallel::SnowParam(workers = nthreads))
    }
}

# Apply function in parallel using the best available backend @param X Vector
# or list to iterate over @param FUN Function to apply @param nthreads Number
# of threads (default: 1) @param SIMPLIFY Whether to simplify results (default:
# TRUE) @param FUN.VALUE Template for vapply (optional)
.bplapply <- function(X, FUN, nthreads = 1, SIMPLIFY = TRUE, FUN.VALUE = NULL) {
    if (nthreads <= 1) {
        # Serial execution
        if (is.null(FUN.VALUE)) {
            return(lapply(X, FUN))
        } else {
            return(unname(vapply(X, FUN, FUN.VALUE = FUN.VALUE)))
        }
    }

    # Parallel execution
    bpparam <- .get_bpparam(nthreads)

    if (is.null(FUN.VALUE)) {
        return(BiocParallel::bplapply(X, FUN, BPPARAM = bpparam))
    } else {
        # Use bplapply and then simplify with vapply
        result_list <- BiocParallel::bplapply(X, FUN, BPPARAM = bpparam)
        return(vapply(result_list, identity, FUN.VALUE = FUN.VALUE))
    }
}

# Apply function over two vectors in parallel @param X First vector or list
# @param Y Second vector or list @param FUN Function to apply (takes two
# arguments) @param nthreads Number of threads (default: 1)
.bpmapply <- function(X, Y, FUN, nthreads = 1) {
    if (nthreads <= 1) {
        return(unname(mapply(FUN, X, Y, SIMPLIFY = FALSE)))
    }

    bpparam <- .get_bpparam(nthreads)
    return(unname(BiocParallel::bpmapply(FUN, X, Y, BPPARAM = bpparam, SIMPLIFY = FALSE)))
}
