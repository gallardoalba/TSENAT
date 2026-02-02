context("Wilcoxon paired behavior")

test_that("wilcoxon paired matches per-row wilcox.test on ordered pairs", {
  # create a small matrix with two genes and two pairs (columns: N,T,N,T)
  # ensure paired differences are not identical to avoid "ties" in differences
  x <- matrix(c(
    1, 2, 3, 6, # gene1 across 4 samples (diffs: -1, -3)
    5, 6, 7, 9 # gene2 across 4 samples (diffs: -1, -2)
  ), nrow = 2, byrow = TRUE)
  samples <- c("Normal", "Tumor", "Normal", "Tumor")

  res <- wilcoxon(x, samples, paired = TRUE, exact = TRUE)

  # compute expected raw p-values by calling wilcox.test per row with
  # paired=TRUE
  expected_raw <- sapply(seq_len(nrow(x)), function(i) {
    wilcox.test(
      x[
        i,
        which(samples == "Normal")
      ],
      x[
        i,
        which(samples == "Tumor")
      ],
      paired = TRUE, exact = TRUE
    )$p.value
  })

  expected_adj <- p.adjust(expected_raw, method = "BH")

  expect_equal(as.numeric(res[, 1]), as.numeric(expected_raw))
  expect_equal(as.numeric(res[, 2]), as.numeric(expected_adj))
})
