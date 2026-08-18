# =============================================================================
# helpers.R — common utilities for the Monte Carlo validation suite
# =============================================================================

`%||%` <- function(a, b) if (is.null(a)) b else a

#' Simulate paired entropy curves under H0
#'
#' Y_{s,c,q} = u_s + f(q) + eps, no condition×q interaction.
#' @param n_subjects number of subjects
#' @param n_q number of q-values
#' @param rho AR(1) within condition
#' @param sd_subject SD of the subject random effect
#' @param seed optional seed
sim_paired_h0 <- function(n_subjects = 12, n_q = 8, rho = 0.5, sd_subject = 1,
    sd_eps = 0.25, seed = NULL) {
    if (!is.null(seed)) set.seed(seed)
    q_seq <- seq(0.1, 2, length.out = n_q)
    fq <- log1p(q_seq)  # common mean trend
    out <- do.call(rbind, lapply(seq_len(n_subjects), function(s) {
        u <- rnorm(1, sd = sd_subject)
        do.call(rbind, lapply(c("A", "B"), function(cond) {
            e <- as.numeric(arima.sim(list(ar = rho), n = n_q, sd = sd_eps))
            data.frame(subject = paste0("S", s), condition = cond,
                q = q_seq, entropy = u + fq + e)
        }))
    }))
    out
}

#' Monte Carlo proportion with binomial CI (Wilson)
binom_ci <- function(k, n) {
    p <- k/n
    z <- qnorm(0.975)
    denom <- 1 + z^2/n
    center <- (p + z^2/(2 * n))/denom
    half <- z * sqrt(p * (1 - p)/n + z^2/(4 * n^2))/denom
    c(lower = center - half, estimate = p, upper = center + half)
}

#' Check whether a proportion CI contains the target value
contains_value <- function(ci, value) ci["lower"] <= value && value <= ci["upper"]
