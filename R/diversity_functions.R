## Minimalized diversity functions: only Tsallis entropy retained

#' Calculate Tsallis entropy for a vector of transcript-level
#' expression values of one gene.
#'
#' @param x Vector of (non-negative) expression values.
#' @param q Tsallis entropy parameter (q > 0). Can be a single value or a numeric vector. Default is 2.
#' @param norm If \code{TRUE}, the entropy values are normalized by their theoretical maximum for the number of transcripts (so values lie in [0,1]).
#' @param what Which quantity to return: \code{"S"} for Tsallis entropy (S_q), \code{"D"} for Hill numbers (D_q), or \code{"both"} for a list with both.
#' @param log_base Base of the logarithm used for Shannon limits and normalization (default: \code{exp(1)}).
#' @export
#' @return A numeric vector (named when length(q) > 1) for \code{what = "S"} or \code{"D"}, or a list with components \code{$S} and \code{$D} when \code{what = "both"}.
#' @details
#' Implements S_q = (1 - sum p^q) / (q - 1) and D_q = (sum p^q)^(1/(1-q)) with the q->1 limits given by Shannon entropy and the exponential of Shannon respectively. Natural logarithms are used for the q->1 limit and normalization.
 #' @examples
 #' Basic usage with a small numeric vector
 #' x <- c(10, 5, 0)
 #' calculate_tsallis_entropy(x, q = c(0.5, 1, 2), norm = TRUE)
calculate_tsallis_entropy <- function(x, q = 2, norm = TRUE, what = c("S", "D", "both"), log_base = exp(1)) {
  what <- match.arg(what)
  if (!is.numeric(q)) stop("q must be numeric.")
  if (any(q <= 0)) stop("q must be greater than 0.")
  if (!is.numeric(x)) stop("x must be numeric")

  n <- length(x)
  # If all counts sum to zero, return NA; allow single-element vectors to proceed
  if (sum(x, na.rm = TRUE) <= 0) {
    if (what == "both") {
      return(list(S = rep(NA_real_, length(q)), D = rep(NA_real_, length(q))))
    }
    return(rep(NA_real_, length(q)))
  }

  p <- x / sum(x)

  tol <- sqrt(.Machine$double.eps)
  S_vec <- vapply(q, function(qi) {
    if (abs(qi - 1) < tol) {
      sh <- -sum(ifelse(p > 0, p * log(p, base = log_base), 0))
      if (norm) sh <- sh / log(n, base = log_base)
      return(sh)
    } else {
      ts <- (1 - sum(p^qi)) / (qi - 1)
      if (norm) {
        max_ts <- (1 - n^(1 - qi)) / (qi - 1)
        ts <- ts / max_ts
      }
      return(ts)
    }
  }, numeric(1))

  D_vec <- vapply(q, function(qi) {
    if (abs(qi - 1) < tol) {
      sh <- -sum(ifelse(p > 0, p * log(p, base = log_base), 0))
      D1 <- (log_base)^(sh)
      return(D1)
    } else {
      spq <- sum(p^qi)
      Dq <- spq^(1 / (1 - qi))
      return(Dq)
    }
  }, numeric(1))

  if (what == "S") {
    out <- S_vec
    if (length(q) > 1) names(out) <- paste0("q=", q)
    if (length(q) == 1) {
      return(unname(out))
    }
    return(out)
  }
  if (what == "D") {
    out <- D_vec
    if (length(q) > 1) names(out) <- paste0("q=", q)
    if (length(q) == 1) {
      return(unname(out))
    }
    return(out)
  }
  # both
  names(S_vec) <- paste0("q=", q)
  names(D_vec) <- paste0("q=", q)
  return(list(S = S_vec, D = D_vec))
}
