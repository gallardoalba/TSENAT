context("generate_plots more tests")

library(SummarizedExperiment)

skip_on_bioc()

test_that(".ptt_prepare_inputs errors when tx2gene missing", {
  skip_if_not_installed("ggplot2")
  counts <- matrix(1:6, nrow = 3)
  rownames(counts) <- paste0("tx", seq_len(nrow(counts)))
  colnames(counts) <- c("S1","S2")
  expect_error(TSENAT:::.ptt_prepare_inputs(counts = counts, readcounts = NULL, samples = c("S1","S2"), coldata = NULL, sample_type_col = "sample_type", tx2gene = NULL, res = NULL, top_n = 2, pseudocount = 1e-6, output_file = NULL), "tx2gene")
})

test_that(".ptt_prepare_inputs returns list with mapping when provided", {
  skip_if_not_installed("ggplot2")
  counts <- matrix(1:6, nrow = 3)
  rownames(counts) <- paste0("tx",1:3)
  colnames(counts) <- c("S1","S2")
  tx2 <- data.frame(Transcript = rownames(counts), Gen = c("G1","G1","G2"), stringsAsFactors = FALSE)
  prep <- TSENAT:::.ptt_prepare_inputs(counts = counts, readcounts = NULL, samples = c("S1","S2"), coldata = NULL, sample_type_col = "sample_type", tx2gene = tx2, res = NULL, top_n = 2, pseudocount = 1e-6, output_file = NULL)
  expect_type(prep, "list")
  expect_true(all(c("counts","samples","mapping","agg_fun") %in% names(prep)))
})

test_that(".ptt_make_plot_for_gene returns ggplot object", {
  skip_if_not_installed("ggplot2")
  counts <- matrix(rpois(6,10), nrow = 3)
  rownames(counts) <- paste0("tx",1:3)
  colnames(counts) <- c("S1","S2")
  mapping <- data.frame(Transcript = rownames(counts), Gen = c("G1","G1","G2"), stringsAsFactors = FALSE)
  agg_fun <- function(x) median(x, na.rm = TRUE)
  p <- TSENAT:::.ptt_make_plot_for_gene("G1", mapping = mapping, counts = counts, samples = c("Normal","Tumor"), top_n = 2, agg_fun = agg_fun, pseudocount = 1e-6, agg_label_unique = "label")
  expect_s3_class(p, "ggplot")
})

test_that(".ptt_combine_plots returns a plot-like object", {
  skip_if_not_installed("ggplot2")
  p1 <- ggplot2::ggplot() + ggplot2::geom_point(mapping = ggplot2::aes(x = 1:3, y = 3:1))
  p2 <- ggplot2::ggplot() + ggplot2::geom_point(mapping = ggplot2::aes(x = 1:3, y = c(1,2,3)))
  out <- TSENAT:::.ptt_combine_plots(list(p1,p2), output_file = NULL, agg_label_unique = "agg")
  expect_true(!is.null(out))
})
