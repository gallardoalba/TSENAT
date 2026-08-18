context("statistical_outputs_golden: Golden Reference Tests")

library(TSENAT)

# Helper: access internal entropy/divergence functions
.ent <- getFromNamespace(".entropy_core", "TSENAT")
.ent_vec <- getFromNamespace(".entropy_vectorized", "TSENAT")
.ent_max <- getFromNamespace(".entropy_max", "TSENAT")

# ---------------------------------------------------------------------------
# SUITE 1: Tsallis entropy — analytical reference values
# ---------------------------------------------------------------------------

test_that("entropy of uniform distribution matches analytical formula at q=0", {
    # q=0 species richness: S_0 = n - 1 for n species with p>0
    p <- c(0.25, 0.25, 0.25, 0.25)  # 4 equally abundant species
    expect_equal(.ent(p, q = 0), 3)  # 4 - 1 = 3
})

test_that("entropy of uniform distribution matches analytical formula at q=1 (Shannon, base e)", {
    # Shannon entropy: -sum(p * ln(p)) = ln(n) for uniform distribution
    p <- rep(1/5, 5)
    expected <- log(5)  # natural log
    expect_equal(.ent(p, q = 1, log_base = exp(1)), expected, tolerance = 1e-10)
})

test_that("entropy of uniform distribution matches analytical formula at q=1 (Shannon, base 2)", {
    p <- rep(1/8, 8)
    expected <- log2(8)  # = 3
    expect_equal(.ent(p, q = 1, log_base = 2), expected, tolerance = 1e-10)
})

test_that("entropy of uniform distribution matches analytical formula at q=2", {
    # Tsallis q=2: S_2 = (1 - sum(p^2)) / (2-1) = 1 - sum(p^2)
    # For uniform: sum(p^2) = n * (1/n)^2 = 1/n, so S_2 = 1 - 1/n
    n <- 6
    p <- rep(1/n, n)
    expected <- 1 - 1/n
    expect_equal(.ent(p, q = 2), expected, tolerance = 1e-10)
})

test_that("entropy of single-species distribution is zero for all q", {
    # One species with 100%: entropy should be 0 regardless of q
    for (q_val in c(0, 0.5, 1, 1.5, 2, 3)) {
        expect_equal(.ent(c(1), q = q_val), 0, tolerance = 1e-10,
            info = paste("q =", q_val))
    }
})

test_that("entropy of (1,0,0,...) distribution is zero for all q", {
    # One dominant + zeros: only 1 species with p>0
    p <- c(1, 0, 0, 0, 0)
    for (q_val in c(0, 0.5, 1, 1.5, 2, 3)) {
        expect_equal(.ent(p, q = q_val), 0, tolerance = 1e-10,
            info = paste("q =", q_val))
    }
})

test_that("entropy of maximally diverse distribution is 1 when normalized", {
    # Uniform distribution normalized: should be exactly 1
    for (n in c(3, 5, 10, 20)) {
        p <- rep(1/n, n)
        for (q_val in c(0.5, 1, 1.5, 2)) {
            val <- .ent(p, q = q_val, norm = TRUE)
            expect_equal(val, 1, tolerance = 1e-8,
                info = paste("n =", n, ", q =", q_val))
        }
    }
})

test_that("entropy increases monotonically with number of equally abundant species", {
    # For uniform distribution with same q, more species → more entropy
    p3 <- rep(1/3, 3)
    p5 <- rep(1/5, 5)
    p10 <- rep(1/10, 10)
    for (q_val in c(0, 0.5, 1, 1.5, 2)) {
        s3 <- .ent(p3, q = q_val)
        s5 <- .ent(p5, q = q_val)
        s10 <- .ent(p10, q = q_val)
        expect_true(s3 < s5, info = paste("q =", q_val, ": 3 vs 5"))
        expect_true(s5 < s10, info = paste("q =", q_val, ": 5 vs 10"))
    }
})

# ---------------------------------------------------------------------------
# SUITE 2: Hand-computed reference values for specific q
# ---------------------------------------------------------------------------

test_that("entropy of (0.7, 0.2, 0.1) at q=0 matches species richness", {
    # 3 species with p>0 → S_0 = 3 - 1 = 2
    expect_equal(.ent(c(0.7, 0.2, 0.1), q = 0), 2)
})

test_that("entropy of (0.7, 0.2, 0.1) at q=1 matches hand-computed Shannon (base e)", {
    p <- c(0.7, 0.2, 0.1)
    # Shannon: -0.7*ln(0.7) - 0.2*ln(0.2) - 0.1*ln(0.1)
    expected <- -(0.7*log(0.7) + 0.2*log(0.2) + 0.1*log(0.1))
    expect_equal(.ent(p, q = 1, log_base = exp(1)), expected, tolerance = 1e-10)
})

test_that("entropy of (0.7, 0.2, 0.1) at q=2 matches hand-computed Tsallis", {
    p <- c(0.7, 0.2, 0.1)
    # S_2 = (1 - (0.7^2 + 0.2^2 + 0.1^2)) / (2-1) = 1 - (0.49 + 0.04 + 0.01) = 0.46
    expected <- 1 - (0.7^2 + 0.2^2 + 0.1^2)  # = 0.46
    expect_equal(.ent(p, q = 2), expected, tolerance = 1e-10)
})

test_that("entropy of (0.5, 0.3, 0.2) at q=0.5 matches hand-computed Tsallis", {
    p <- c(0.5, 0.3, 0.2)
    # S_{0.5} = (1 - sum(p^0.5)) / (0.5 - 1) = (1 - sum(sqrt(p))) / (-0.5)
    expected <- (1 - (sqrt(0.5) + sqrt(0.3) + sqrt(0.2))) / (0.5 - 1)
    expect_equal(.ent(p, q = 0.5), expected, tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# SUITE 3: Maximum entropy reference values
# ---------------------------------------------------------------------------

test_that(".entropy_max at q=0 returns n-1", {
    expect_equal(.ent_max(5, q = 0), 4)
    expect_equal(.ent_max(1, q = 0), 0)
    expect_equal(.ent_max(10, q = 0), 9)
})

test_that(".entropy_max at q=1 returns log(n) / log(log_base)", {
    expect_equal(.ent_max(8, q = 1, log_base = 2), 3, tolerance = 1e-10)  # log2(8) = 3
    expect_equal(.ent_max(10, q = 1, log_base = exp(1)), log(10), tolerance = 1e-10)
})

test_that(".entropy_max at q=2 returns 1 - 1/n", {
    # H_max = (1 - n^(1-2)) / (2-1) = (1 - n^(-1)) = 1 - 1/n
    expect_equal(.ent_max(4, q = 2), 1 - 1/4, tolerance = 1e-10)
    expect_equal(.ent_max(10, q = 2), 1 - 1/10, tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# SUITE 4: Vectorized entropy consistency
# ---------------------------------------------------------------------------

test_that("vectorized entropy handles uniform proportions correctly", {
    # Equal-length rows: 4 species vs 8 species, both uniform
    counts <- rbind(
        c(1, 1, 1, 1, 0, 0, 0, 0),   # 4 species with data + 4 zeros
        c(1, 1, 1, 1, 1, 1, 1, 1)    # 8 species
    )
    res <- .ent_vec(counts, q = 2, pseudocount = 0)
    # Row 1: 4 uniform p>0 → S_2 = 1 - 1/4 = 0.75
    # Row 2: 8 uniform p>0 → S_2 = 1 - 1/8 = 0.875
    expect_equal(res[1], 1 - 1/4, tolerance = 1e-10)
    expect_equal(res[2], 1 - 1/8, tolerance = 1e-10)
    expect_true(res[2] > res[1])
})

# ---------------------------------------------------------------------------
# SUITE 5: Tsallis divergence — identical distributions → 0
# ---------------------------------------------------------------------------

test_that("divergence between identical distributions is zero for all q", {
    # Access internal divergence function
    .ts_div <- getFromNamespace(".tsallis_divergence_scalar", "TSENAT")
    
    x <- c(10, 20, 30, 40)
    y <- c(10, 20, 30, 40)  # identical
    
    for (q_val in c(0.5, 1, 1.5, 2, 3)) {
        d <- .ts_div(x, y, q_val, pseudocount = 0.5, log_base = exp(1))
        expect_equal(d, 0, tolerance = 1e-8,
            info = paste("q =", q_val))
    }
})

test_that("divergence between identical count vectors is zero", {
    .ts_div <- getFromNamespace(".tsallis_divergence_scalar", "TSENAT")
    
    x <- c(10, 20, 30, 40)
    y <- c(10, 20, 30, 40)  # exact same
    
    for (q_val in c(0.5, 1, 1.5, 2)) {
        d <- .ts_div(x, y, q_val, pseudocount = 0.5, log_base = exp(1))
        expect_equal(d, 0, tolerance = 1e-8,
            info = paste("q =", q_val))
    }
})

test_that("divergence between different distributions is positive", {
    .ts_div <- getFromNamespace(".tsallis_divergence_scalar", "TSENAT")
    
    x <- c(100, 10, 10)
    y <- c(10, 100, 10)  # different
    
    for (q_val in c(0.5, 1, 1.5, 2)) {
        d <- .ts_div(x, y, q_val, pseudocount = 0.5, log_base = exp(1))
        expect_true(d > 0, info = paste("q =", q_val))
    }
})

# ---------------------------------------------------------------------------
# SUITE 6: Edge case robustness
# ---------------------------------------------------------------------------

test_that("entropy of empty vector returns NA", {
    expect_true(is.na(.ent(numeric(0))))
})

test_that("entropy of all-NA returns NA", {
    expect_true(is.na(.ent(c(NA, NA, NA))))
})

test_that("entropy rejects negative proportions instead of silently filtering", {
    # Negative values are invalid input: silently discarding them would
    # compute entropy on a different abundance vector than the user supplied
    # (audit 2026-08-17).
    expect_error(.ent(c(-0.1, 0.5, 0.6)), "non-negative")
})

test_that("entropy rejects negative q", {
    expect_error(.ent(c(0.5, 0.5), q = -1), "must be non-negative")
})

test_that("entropy rejects vector q", {
    expect_error(.ent(c(0.5, 0.5), q = c(1, 2)), "must be a scalar")
})

# ---------------------------------------------------------------------------
# SUITE 7: log_base consistency
# ---------------------------------------------------------------------------

test_that("Shannon entropy with base 2 produces consistent values", {
    p <- c(0.5, 0.25, 0.25)
    # Shannon base 2: -0.5*log2(0.5) - 0.25*log2(0.25) - 0.25*log2(0.25)
    # = 0.5*1 + 0.25*2 + 0.25*2 = 0.5 + 0.5 + 0.5 = 1.5
    expected <- 1.5
    expect_equal(.ent(p, q = 1, log_base = 2), expected, tolerance = 1e-10)
})

test_that("Shannon entropy with base 10 produces consistent values", {
    p <- c(1/10, 9/10)
    # Shannon base 10: -0.1*log10(0.1) - 0.9*log10(0.9)
    expected <- -(0.1*log10(0.1) + 0.9*log10(0.9))
    expect_equal(.ent(p, q = 1, log_base = 10), expected, tolerance = 1e-10)
})

test_that("Tsallis entropy is independent of log_base for q≠1", {
    # Tsallis formula is scale-invariant — log_base should not affect result
    p <- c(0.3, 0.3, 0.4)
    s_e <- .ent(p, q = 2, log_base = exp(1))
    s_2 <- .ent(p, q = 2, log_base = 2)
    s_10 <- .ent(p, q = 2, log_base = 10)
    expect_equal(s_e, s_2, tolerance = 1e-10)
    expect_equal(s_e, s_10, tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# SUITE 8: C++ consistency
# ---------------------------------------------------------------------------

test_that("R and C++ entropy implementations agree", {
    set.seed(123)
    for (q_val in c(0, 0.5, 1, 1.5, 2, 3)) {
        p <- runif(8)
        p <- p / sum(p)
        
        r_val <- .ent(p, q = q_val)
        cpp_val <- entropy_cpp(p, q = q_val, normalize = FALSE)
        
        expect_equal(r_val, cpp_val, tolerance = 1e-8,
            info = paste("q =", q_val))
    }
})
