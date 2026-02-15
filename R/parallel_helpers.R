# Helper functions for parallel processing using BiocParallel and parallel
# package

# Select and initialize parallel backend @param nthreads Number of threads to
# use (default: 1) @return BiocParallel BPPARAM object
.tsenat_get_bpparam <- function(nthreads = 1) {
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
.tsenat_bplapply <- function(X, FUN, nthreads = 1, SIMPLIFY = TRUE, FUN.VALUE = NULL) {
    if (nthreads <= 1) {
        # Serial execution
        if (is.null(FUN.VALUE)) {
            return(lapply(X, FUN))
        } else {
            return(unname(vapply(X, FUN, FUN.VALUE = FUN.VALUE)))
        }
    }

    # Parallel execution
    bpparam <- .tsenat_get_bpparam(nthreads)

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
.tsenat_bpmapply <- function(X, Y, FUN, nthreads = 1) {
    if (nthreads <= 1) {
        return(unname(mapply(FUN, X, Y, SIMPLIFY = FALSE)))
    }

    bpparam <- .tsenat_get_bpparam(nthreads)
    return(unname(BiocParallel::bpmapply(FUN, X, Y, BPPARAM = bpparam, SIMPLIFY = FALSE)))
}
