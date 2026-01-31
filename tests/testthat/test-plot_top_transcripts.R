testthat::test_that("plot_top_transcripts returns ggplot for synthetic data", {
  skip_on_cran()
  set.seed(42)
  counts <- matrix(rpois(3 * 8, lambda = 20), nrow = 3)
  rownames(counts) <- paste0("tx", 1:3)
  colnames(counts) <- paste0("S", 1:8)
  samples <- c(rep("Normal", 4), rep("Tumor", 4))
  # create simple tx2gene mapping
  tx2 <- data.frame(Transcript = rownames(counts), Gen = rep("GENE1", 3), stringsAsFactors = FALSE)

  p <- plot_top_transcripts(counts, gene = "GENE1", samples = samples, tx2gene = tx2, top_n = 2)
  testthat::expect_s3_class(p, "ggplot")
})
