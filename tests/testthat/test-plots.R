# Comprehensive testing of all plotting functions
# Tests plot_diversity_density, plot_mean_violin, plot_ma, 
# plot_top_transcripts, plot_volcano, plot_tsallis_q_curve,
# plot_tsallis_gene_profile, plot_tsallis_density_multq, plot_tsallis_violin_multq

library(TSENAT)
skip_on_bioc()

test_that("plot_diversity_density returns ggplot object with valid data", {
    skip_if_not_installed("SummarizedExperiment")
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("tidyr")
    skip_if_not_installed("dplyr")

    library(SummarizedExperiment)
    library(ggplot2)

    # construct minimal SummarizedExperiment
    mat <- matrix(runif(20), nrow = 5, ncol = 4)
    colnames(mat) <- c("S1_N", "S2_T", "S3_N", "S4_T")
    rownames(mat) <- paste0("G", 1:5)
    rowData_df <- S4Vectors::DataFrame(genes = rownames(mat))
    colData_df <- S4Vectors::DataFrame(samples = colnames(mat))
    se <- SummarizedExperiment(
        assays = list(diversity = mat),
        rowData = rowData_df,
        colData = colData_df
    )

    p <- plot_diversity_density(se)
    expect_s3_class(p, "gg")
    expect_s3_class(p, "ggplot")
})

test_that("plot_mean_violin returns ggplot object", {
    skip_if_not_installed("SummarizedExperiment")
    skip_if_not_installed("ggplot2")

    library(SummarizedExperiment)
    library(ggplot2)

    mat <- matrix(runif(20), nrow = 5, ncol = 4)
    colnames(mat) <- c("S1_N", "S2_T", "S3_N", "S4_T")
    rownames(mat) <- paste0("G", 1:5)
    rowData_df <- S4Vectors::DataFrame(genes = rownames(mat))
    colData_df <- S4Vectors::DataFrame(samples = colnames(mat))
    se <- SummarizedExperiment(
        assays = list(diversity = mat),
        rowData = rowData_df,
        colData = colData_df
    )

    p <- plot_mean_violin(se)
    expect_s3_class(p, "gg")
    expect_s3_class(p, "ggplot")
})

test_that("plot_ma returns ggplot object with mean columns", {
    skip_if_not_installed("ggplot2")
    library(ggplot2)

    df <- data.frame(
        Gene = paste0("G", 1:10),
        A_mean = runif(10),
        B_mean = runif(10),
        log2_fold_change = rnorm(10),
        adjusted_p_values = runif(10)
    )
    p <- plot_ma_tsallis(df)
    expect_s3_class(p, "gg")
    expect_s3_class(p, "ggplot")
})

test_that("plot_ma returns ggplot object with median columns", {
    skip_if_not_installed("ggplot2")
    library(ggplot2)

    df <- data.frame(
        Gene = paste0("G", 1:10),
        A_median = runif(10),
        B_median = runif(10),
        log2_fold_change = rnorm(10),
        adjusted_p_values = runif(10)
    )
    p <- plot_ma_tsallis(df)
    expect_s3_class(p, "gg")
})

test_that("plot_ma errors on mixed mean/median columns", {
    skip_if_not_installed("ggplot2")
    library(ggplot2)

    df <- data.frame(
        Gene = paste0("G", 1:10),
        A_mean = runif(10),
        B_median = runif(10),
        log2_fold_change = rnorm(10),
        adjusted_p_values = runif(10)
    )
    expect_error(plot_ma_tsallis(df), "Could not find two mean or two median columns")
})

test_that("plot_top_transcripts returns ggplot for synthetic data", {
    skip_if_not_installed("ggplot2")
    library(ggplot2)

    set.seed(42)
    counts <- matrix(rpois(3 * 8, lambda = 20), nrow = 3)
    rownames(counts) <- paste0("tx", 1:3)
    colnames(counts) <- paste0("S", 1:8)
    samples <- c(rep("Normal", 4), rep("Tumor", 4))
    # create simple tx2gene mapping
    tx2 <- data.frame(
        Transcript = rownames(counts),
        Gen = rep("GENE1", 3),
        stringsAsFactors = FALSE
    )

    p <- plot_top_transcripts(counts,
        gene = "GENE1",
        samples = samples,
        tx2gene = tx2,
        top_n = 2
    )
    expect_s3_class(p, "ggplot")
})

test_that("plot_top_transcripts selects genes from res when gene is NULL", {
    skip_if_not_installed("ggplot2")
    library(ggplot2)

    set.seed(42)
    counts <- matrix(rpois(9 * 4, lambda = 20), nrow = 9)
    rownames(counts) <- paste0("tx", 1:9)
    colnames(counts) <- paste0("S", 1:4)
    samples <- c(rep("Normal", 2), rep("Tumor", 2))

    tx2 <- data.frame(
        Transcript = rownames(counts),
        Gen = rep(paste0("G", 1:3), each = 3),
        stringsAsFactors = FALSE
    )

    res <- data.frame(genes = paste0("G", 1:3), adjusted_p_values = c(0.01, 0.05, 0.2), stringsAsFactors = FALSE)

    p <- plot_top_transcripts(counts, res = res, tx2gene = tx2, samples = samples, top_n = 2)
    expect_s3_class(p, "ggplot")
})

test_that("plot_volcano returns a ggplot and annotates top genes", {
    skip_if_not_installed("ggplot2")
    library(ggplot2)

    set.seed(42)
    n <- 20
    df <- data.frame(
        genes = paste0("gene", seq_len(n)),
        mean_difference = rnorm(n),
        adjusted_p_values = p.adjust(runif(n))
    )

    p <- plot_volcano(df,
        x_col = "mean_difference",
        padj_col = "adjusted_p_values",
        top_n = 3
    )
    expect_s3_class(p, "ggplot")
    # building the plot should not error
    ggplot2::ggplot_build(p)
})

test_that("plot_volcano with custom columns", {
    skip_if_not_installed("ggplot2")
    library(ggplot2)

    set.seed(42)
    n <- 15
    df <- data.frame(
        genes = paste0("gene", seq_len(n)),
        logFC = rnorm(n),
        pval = p.adjust(runif(n))
    )

    p <- plot_volcano(df,
        x_col = "logFC",
        padj_col = "pval",
        top_n = 2
    )
    expect_s3_class(p, "ggplot")
})

test_that("plot_tsallis_q_curve returns ggplot with valid SE", {
    skip_if_not_installed("SummarizedExperiment")
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("tidyr")
    skip_if_not_installed("dplyr")

    library(SummarizedExperiment)
    library(ggplot2)

    set.seed(1)
    readcounts <- matrix(rpois(30 * 3, lambda = 10), nrow = 30, ncol = 3)
    colnames(readcounts) <- c("S1_N", "S2_T", "S3_N")
    genes <- rep(paste0("G", 1:10), length.out = nrow(readcounts))

    qvals <- seq(0.01, 0.05, by = 0.01)
    ts_se <- calculate_diversity(readcounts, genes, q = qvals, norm = TRUE)

    coldata_df <- data.frame(
        Sample = c("S1_N", "S2_T", "S3_N"),
        Condition = c("Normal", "Tumor", "Normal"),
        stringsAsFactors = FALSE
    )

    ts_se <- map_metadata(ts_se, coldata_df)

    p <- plot_tsallis_q_curve(ts_se)
    expect_true(inherits(p, "ggplot"))
})

test_that("plot_tsallis_gene_profile returns ggplot for single gene", {
    skip_if_not_installed("SummarizedExperiment")
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("tidyr")
    skip_if_not_installed("dplyr")

    library(SummarizedExperiment)
    library(ggplot2)

    set.seed(123)
    readcounts <- matrix(rpois(30 * 4, lambda = 15), nrow = 30, ncol = 4)
    colnames(readcounts) <- c("S1_N", "S2_N", "S3_T", "S4_T")
    genes <- rep(paste0("G", 1:10), length.out = nrow(readcounts))

    qvals <- seq(0.01, 0.05, by = 0.02)
    ts_se <- calculate_diversity(readcounts, genes, q = qvals, norm = TRUE)

    coldata_df <- data.frame(
        Sample = c("S1_N", "S2_N", "S3_T", "S4_T"),
        Condition = c("Normal", "Normal", "Tumor", "Tumor"),
        stringsAsFactors = FALSE
    )

    ts_se <- map_metadata(ts_se, coldata_df)

    p <- plot_tsallis_gene_profile(ts_se, gene = "G1")
    expect_s3_class(p, "ggplot")
})

test_that("plot_tsallis_density_multq returns ggplot", {
    skip_if_not_installed("SummarizedExperiment")
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("tidyr")
    skip_if_not_installed("dplyr")

    library(SummarizedExperiment)
    library(ggplot2)

    set.seed(456)
    readcounts <- matrix(rpois(25 * 4, lambda = 12), nrow = 25, ncol = 4)
    colnames(readcounts) <- c("S1_N", "S2_N", "S3_T", "S4_T")
    genes <- rep(paste0("G", 1:5), length.out = nrow(readcounts))

    qvals <- seq(0.01, 0.1, by = 0.03)
    ts_se <- calculate_diversity(readcounts, genes, q = qvals, norm = TRUE)

    coldata_df <- data.frame(
        Sample = c("S1_N", "S2_N", "S3_T", "S4_T"),
        Condition = c("Normal", "Normal", "Tumor", "Tumor"),
        stringsAsFactors = FALSE
    )

    ts_se <- map_metadata(ts_se, coldata_df)

    p <- plot_tsallis_density_multq(ts_se)
    expect_s3_class(p, "ggplot")
})

test_that("plot_tsallis_violin_multq returns ggplot", {
    skip_if_not_installed("SummarizedExperiment")
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("tidyr")
    skip_if_not_installed("dplyr")

    library(SummarizedExperiment)
    library(ggplot2)

    set.seed(789)
    readcounts <- matrix(rpois(20 * 4, lambda = 18), nrow = 20, ncol = 4)
    colnames(readcounts) <- c("S1_N", "S2_N", "S3_T", "S4_T")
    genes <- rep(paste0("G", 1:4), length.out = nrow(readcounts))

    qvals <- c(0.01, 0.05, 0.1)
    ts_se <- calculate_diversity(readcounts, genes, q = qvals, norm = TRUE)

    coldata_df <- data.frame(
        Sample = c("S1_N", "S2_N", "S3_T", "S4_T"),
        Condition = c("Normal", "Normal", "Tumor", "Tumor"),
        stringsAsFactors = FALSE
    )

    ts_se <- map_metadata(ts_se, coldata_df)

    p <- plot_tsallis_violin_multq(ts_se)
    expect_s3_class(p, "ggplot")
})

test_that("All plot functions produce buildable ggplot objects", {
    skip_if_not_installed("SummarizedExperiment")
    skip_if_not_installed("ggplot2")

    library(SummarizedExperiment)
    library(ggplot2)

    # Create test data
    mat <- matrix(runif(15), nrow = 5, ncol = 3)
    colnames(mat) <- c("S1_N", "S2_T", "S3_N")
    rownames(mat) <- paste0("G", 1:5)
    rowData_df <- S4Vectors::DataFrame(genes = rownames(mat))
    colData_df <- S4Vectors::DataFrame(samples = colnames(mat))
    se <- SummarizedExperiment(
        assays = list(diversity = mat),
        rowData = rowData_df,
        colData = colData_df
    )

    # Test that plots can be built without error
    expect_no_error(ggplot_build(plot_diversity_density(se)))
    expect_no_error(ggplot_build(plot_mean_violin(se)))
})

library(SummarizedExperiment)

test_that("infer_samples_from_se finds sample_type column and falls back", {
  mat <- matrix(runif(6), nrow = 3, ncol = 2)
  colnames(mat) <- c("a","b")
  se <- SummarizedExperiment(assays = list(diversity = mat), colData = S4Vectors::DataFrame(sample_type = c("X","Y")))
  samples <- infer_samples_from_se(se)
  expect_equal(samples, c("X","Y"))

  # If no sample_type, but a binary column exists
  se2 <- SummarizedExperiment(assays = list(diversity = mat), colData = S4Vectors::DataFrame(cond = c("A","B")))
  samples2 <- infer_samples_from_se(se2)
  expect_equal(samples2, c("A","B"))
})

test_that("get_readcounts_from_se accepts readcounts in metadata, assays and file", {
  mat <- matrix(1:6, nrow = 3)
  rownames(mat) <- paste0("tx", 1:3)
  se <- SummarizedExperiment(assays = list(dummy = matrix(0, nrow = 3, ncol = 2)))
  S4Vectors::metadata(se)$readcounts <- mat
  rc <- get_readcounts_from_se(se)
  expect_true(is.matrix(rc))
  expect_equal(rownames(rc), rownames(mat))

  # if first assay used
  se2 <- SummarizedExperiment(assays = list(readcounts = mat))
  rc2 <- get_readcounts_from_se(se2)
  expect_true(is.matrix(rc2))

  # file input: write temporary table
  tmpf <- tempfile(fileext = ".tsv")
  df <- data.frame(tx = rownames(mat), mat, stringsAsFactors = FALSE)
  write.table(df, file = tmpf, sep = "\t", row.names = FALSE, quote = FALSE)
  rcf <- get_readcounts_from_se(se, readcounts_arg = tmpf)
  expect_true(is.matrix(rcf))
})

test_that("get_tx2gene_from_se returns mapping from metadata, rowData or rownames", {
  mat <- matrix(1:6, nrow = 3)
  rownames(mat) <- paste0("tx", 1:3)
  se <- SummarizedExperiment(assays = list(diversity = mat))
  md <- list(tx2gene = data.frame(Transcript = rownames(mat), Gen = c("g1","g1","g2"), stringsAsFactors = FALSE))
  S4Vectors::metadata(se) <- md
  out <- get_tx2gene_from_se(se, readcounts_mat = mat)
  expect_equal(out$type, "vector")
  expect_equal(length(out$mapping), nrow(mat))

  # rowData case
  se2 <- SummarizedExperiment(assays = list(diversity = mat), rowData = S4Vectors::DataFrame(genes = c("g1","g1","g2")))
  out2 <- get_tx2gene_from_se(se2, readcounts_mat = mat)
  expect_equal(out2$type, "vector")
  expect_equal(length(out2$mapping), nrow(mat))

  # fallback to rownames
  se3 <- SummarizedExperiment(assays = list(diversity = mat))
  out3 <- get_tx2gene_from_se(se3, readcounts_mat = mat)
  expect_equal(out3$type, "vector")
})

test_that("validate_control_in_samples picks 'Normal' when present or first level", {
  samples <- c("Tumor","Normal","Tumor")
  expect_equal(validate_control_in_samples(NULL, samples), "Normal")
  samples2 <- c("A","B")
  expect_message(chosen <- validate_control_in_samples(NULL, samples2))
  expect_true(chosen %in% samples2)
  expect_equal(validate_control_in_samples("B", samples2), "B")
})

test_that("plot_ma_expression_impl works with precomputed fc df and SummarizedExperiment counts", {
  skip_if_not_installed("ggplot2")
  # Prepare x (diff results)
  x <- data.frame(genes = paste0("g", 1:5), mean = runif(5), log2_fold_change = rnorm(5), adjusted_p_values = runif(5))
  # precomputed fc as matrix/data.frame
  fc <- data.frame(genes = paste0("g", 1:5), log2_fold_change = rnorm(5), stringsAsFactors = FALSE)
  p1 <- plot_ma_expression_impl(x, se = fc)
  expect_s3_class(p1, "ggplot")

  # Now using SummarizedExperiment counts
  counts <- matrix(rpois(15, 10), nrow = 5)
  rownames(counts) <- paste0("g", 1:5)
  colnames(counts) <- paste0("S", 1:3)
  colData_df <- S4Vectors::DataFrame(sample = colnames(counts), sample_type = c("Normal","Tumor","Normal"))
  se <- SummarizedExperiment(assays = list(readcounts = counts), rowData = S4Vectors::DataFrame(genes = rownames(counts)), colData = colData_df)
  x2 <- data.frame(genes = rownames(counts), mean = runif(5), adjusted_p_values = runif(5))
  p2 <- plot_ma_expression_impl(x2, se = se, samples = c("Normal","Tumor","Normal"))
  expect_s3_class(p2, "ggplot")
})

test_that(".plot_ma_core errors when fold-change column missing or x axis missing", {
  skip_if_not_installed("ggplot2")
  df <- data.frame(genes = paste0("g", 1:4), val = runif(4))
  expect_error(.plot_ma_core(df), "Could not find a fold-change column")
  df2 <- data.frame(genes = paste0("g", 1:4), log2_fold_change = rnorm(4))
  expect_s3_class(.plot_ma_core(df2), "ggplot")
})

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

context("generate_plots extra tests")

library(SummarizedExperiment)

skip_on_bioc()

test_that("plot_top_transcripts works on simple matrix input", {
  skip_if_not_installed("ggplot2")
  tx_counts <- matrix(sample(1:100, 24, replace = TRUE), nrow = 6)
  rownames(tx_counts) <- paste0("tx", seq_len(nrow(tx_counts)))
  colnames(tx_counts) <- paste0("S", seq_len(ncol(tx_counts)))
  tx2gene <- data.frame(Transcript = rownames(tx_counts), Gen = rep(paste0("G", seq_len(3)), each = 2), stringsAsFactors = FALSE)
  samples <- rep(c("Normal", "Tumor"), length.out = ncol(tx_counts))
  p <- TSENAT:::plot_top_transcripts(tx_counts, gene = c("G1", "G2"), samples = samples, tx2gene = tx2gene, top_n = 2)
  expect_true(!is.null(p))
  expect_true(inherits(p, "ggplot") || inherits(p, "patchwork") || inherits(p, "gtable") || inherits(p, "ggarrange"))
})

test_that("plot_top_transcripts selects genes from res when gene is NULL", {
  skip_if_not_installed("ggplot2")
  tx_counts <- matrix(sample(1:100, 36, replace = TRUE), nrow = 9)
  rownames(tx_counts) <- paste0("tx", seq_len(nrow(tx_counts)))
  colnames(tx_counts) <- paste0("S", seq_len(ncol(tx_counts)))
  tx2gene <- data.frame(Transcript = rownames(tx_counts), Gen = rep(paste0("G", seq_len(3)), each = 3), stringsAsFactors = FALSE)
  samples <- rep(c("Normal", "Tumor"), length.out = ncol(tx_counts))
  res <- data.frame(genes = paste0("G", 1:3), adjusted_p_values = c(0.01, 0.05, 0.2), stringsAsFactors = FALSE)
  p <- TSENAT:::plot_top_transcripts(tx_counts, res = res, tx2gene = tx2gene, samples = samples, top_n = 2)
  expect_true(!is.null(p))
})

test_that("plot_top_transcripts errors when counts lack rownames", {
  skip_if_not_installed("ggplot2")
  mat <- matrix(1:6, nrow = 2)
  expect_error(TSENAT:::plot_top_transcripts(mat, gene = "G1", tx2gene = data.frame(Transcript = c("a","b"), Gen = c("G1","G1"))), "counts.*rownames")
})

test_that("plot_tsallis_gene_profile returns ggplot for simple SE", {
  skip_if_not_installed(c("ggplot2", "SummarizedExperiment"))
  mat <- matrix(runif(6), nrow = 2)
  rownames(mat) <- c("g1", "g2")
  colnames(mat) <- c("S1_q=0.1", "S1_q=1", "S2_q=0.1")
  se <- SummarizedExperiment::SummarizedExperiment(assays = list(diversity = mat))
  # align colData to assay columns (length must match)
  cd <- S4Vectors::DataFrame(sample_type = c("A", "A", "B"))
  rownames(cd) <- colnames(mat)
  SummarizedExperiment::colData(se) <- cd
  p <- TSENAT::plot_tsallis_gene_profile(se, gene = "g1")
  expect_s3_class(p, "ggplot")
})

test_that("plot_diversity_density and plot_mean_violin return ggplot", {
  skip_if_not_installed(c("ggplot2", "SummarizedExperiment"))
  mat <- matrix(rnorm(20), nrow = 5)
  rownames(mat) <- paste0("g", 1:5)
  colnames(mat) <- paste0("S", 1:4)
  se <- SummarizedExperiment::SummarizedExperiment(assays = list(diversity = mat))
  cd <- S4Vectors::DataFrame(sample_type = c("A", "A", "B", "B"))
  rownames(cd) <- colnames(mat)
  SummarizedExperiment::colData(se) <- cd
  d1 <- TSENAT::plot_diversity_density(se)
  d2 <- TSENAT::plot_mean_violin(se)
  expect_s3_class(d1, "ggplot")
  expect_s3_class(d2, "ggplot")
})

test_that("plot_ma_tsallis and plot_ma_expression_impl handle simple inputs", {
  skip_if_not_installed("ggplot2")
  x <- data.frame(genes = paste0("g", 1:6), mean = runif(6), log2_fold_change = rnorm(6))
  p1 <- TSENAT::plot_ma_tsallis(x)
  expect_s3_class(p1, "ggplot")

  # precomputed fold-change data.frame
  fc <- data.frame(genes = paste0("g", 1:6), log2_fold_change = rnorm(6), stringsAsFactors = FALSE)
  p2 <- TSENAT:::plot_ma_expression_impl(x = x, se = fc)
  expect_s3_class(p2, "ggplot")
})

test_that("plot_tsallis_q_curve and multq plots return ggplot", {
  skip_if_not_installed(c("ggplot2", "SummarizedExperiment"))
  mat <- matrix(runif(12), nrow = 3)
  rownames(mat) <- paste0("g", 1:3)
  colnames(mat) <- c("S1_q=0.1", "S1_q=1", "S2_q=0.1", "S2_q=1")
  se <- SummarizedExperiment::SummarizedExperiment(assays = list(diversity = mat))
  cd <- S4Vectors::DataFrame(sample_type = c("A", "A", "B", "B"))
  rownames(cd) <- colnames(mat)
  SummarizedExperiment::colData(se) <- cd
  pq <- TSENAT:::plot_tsallis_q_curve(se)
  expect_s3_class(pq, "ggplot")
  pv <- TSENAT::plot_tsallis_violin_multq(se)
  pd <- TSENAT::plot_tsallis_density_multq(se)
  expect_s3_class(pv, "ggplot")
  expect_s3_class(pd, "ggplot")
})

test_that("plot_volcano auto-detects x_col and returns ggplot", {
  skip_if_not_installed("ggplot2")
  df <- data.frame(gene = paste0("g", 1:10), mean_difference = rnorm(10), adjusted_p_values = runif(10))
  p <- TSENAT::plot_volcano(df)
  expect_s3_class(p, "ggplot")
})
