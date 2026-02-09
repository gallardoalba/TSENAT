context("plot_top_transcripts")

library(SummarizedExperiment)

test_that("plot_top_transcripts works on simple matrix input", {
  tx_counts <- matrix(sample(1:100, 24, replace = TRUE), nrow = 6)
  rownames(tx_counts) <- paste0("tx", seq_len(nrow(tx_counts)))
  colnames(tx_counts) <- paste0("S", seq_len(ncol(tx_counts)))

  tx2gene <- data.frame(Transcript = rownames(tx_counts), Gen = rep(paste0("G", seq_len(3)), each = 2), stringsAsFactors = FALSE)
  samples <- rep(c("Normal", "Tumor"), length.out = ncol(tx_counts))

  p <- plot_top_transcripts(tx_counts, gene = c("G1", "G2"), samples = samples, tx2gene = tx2gene, top_n = 2)
  expect_true(!is.null(p))
  # expect ggplot object or patchwork
  expect_true(inherits(p, "ggplot") || inherits(p, "patchwork") || inherits(p, "gtable") || inherits(p, "ggarrange"))
})

test_that("plot_top_transcripts selects genes from res when gene is NULL", {
  tx_counts <- matrix(sample(1:100, 36, replace = TRUE), nrow = 9)
  rownames(tx_counts) <- paste0("tx", seq_len(nrow(tx_counts)))
  colnames(tx_counts) <- paste0("S", seq_len(ncol(tx_counts)))

  # Make tx2gene mapping: three transcripts per gene
  tx2gene <- data.frame(Transcript = rownames(tx_counts), Gen = rep(paste0("G", seq_len(3)), each = 3), stringsAsFactors = FALSE)
  samples <- rep(c("Normal", "Tumor"), length.out = ncol(tx_counts))

  # create fake res with genes and p-values
  res <- data.frame(genes = paste0("G", seq_len(3)), adjusted_p_values = c(0.01, 0.05, 0.2), stringsAsFactors = FALSE)

  p <- plot_top_transcripts(tx_counts, res = res, tx2gene = tx2gene, samples = samples, top_n = 2)
  expect_true(!is.null(p))
})

test_that("plot_top_transcripts errors when counts lack rownames", {
  mat <- matrix(1:6, nrow = 2)
  expect_error(plot_top_transcripts(mat, gene = "G1", tx2gene = data.frame(Transcript = c("a","b"), Gen = c("G1","G1"))), "counts.*rownames")
})
