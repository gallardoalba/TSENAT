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

    # exercise gene = NULL + lm_res selection path
    lm_res <- data.frame(gene = c("G1", "G2"), adj_p_interaction = c(0.01, 0.2), stringsAsFactors = FALSE)
    p2 <- plot_tsallis_gene_profile(ts_se, gene = NULL, lm_res = lm_res, n_top = 1)
    expect_s3_class(p2, "ggplot")

    # show_samples = TRUE branch
    p3 <- plot_tsallis_gene_profile(ts_se, gene = "G1", show_samples = TRUE)
    expect_s3_class(p3, "ggplot")
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
    colnames(mat) <- c("a", "b")
    se <- SummarizedExperiment(assays = list(diversity = mat), colData = S4Vectors::DataFrame(sample_type = c("X", "Y")))
    samples <- infer_samples_from_se(se)
    expect_equal(samples, c("X", "Y"))

    # If no sample_type, but a binary column exists
    se2 <- SummarizedExperiment(assays = list(diversity = mat), colData = S4Vectors::DataFrame(cond = c("A", "B")))
    samples2 <- infer_samples_from_se(se2)
    expect_equal(samples2, c("A", "B"))
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
    md <- list(tx2gene = data.frame(Transcript = rownames(mat), Gen = c("g1", "g1", "g2"), stringsAsFactors = FALSE))
    S4Vectors::metadata(se) <- md
    out <- get_tx2gene_from_se(se, readcounts_mat = mat)
    expect_equal(out$type, "vector")
    expect_equal(length(out$mapping), nrow(mat))

    # rowData case
    se2 <- SummarizedExperiment(assays = list(diversity = mat), rowData = S4Vectors::DataFrame(genes = c("g1", "g1", "g2")))
    out2 <- get_tx2gene_from_se(se2, readcounts_mat = mat)
    expect_equal(out2$type, "vector")
    expect_equal(length(out2$mapping), nrow(mat))

    # fallback to rownames
    se3 <- SummarizedExperiment(assays = list(diversity = mat))
    out3 <- get_tx2gene_from_se(se3, readcounts_mat = mat)
    expect_equal(out3$type, "vector")
})

test_that("validate_control_in_samples picks 'Normal' when present or first level", {
    samples <- c("Tumor", "Normal", "Tumor")
    expect_equal(validate_control_in_samples(NULL, samples), "Normal")
    samples2 <- c("A", "B")
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
    colData_df <- S4Vectors::DataFrame(sample = colnames(counts), sample_type = c("Normal", "Tumor", "Normal"))
    se <- SummarizedExperiment(assays = list(readcounts = counts), rowData = S4Vectors::DataFrame(genes = rownames(counts)), colData = colData_df)
    x2 <- data.frame(genes = rownames(counts), mean = runif(5), adjusted_p_values = runif(5))
    p2 <- plot_ma_expression_impl(x2, se = se, samples = c("Normal", "Tumor", "Normal"))
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
    expect_error(plot_top_transcripts(mat, gene = "G1", tx2gene = data.frame(Transcript = c("a", "b"), Gen = c("G1", "G1"))), "counts.*rownames")
})

context("generate_plots more tests")

library(SummarizedExperiment)

skip_on_bioc()

test_that(".ptt_prepare_inputs errors when tx2gene missing", {
    skip_if_not_installed("ggplot2")
    counts <- matrix(1:6, nrow = 3)
    rownames(counts) <- paste0("tx", seq_len(nrow(counts)))
    colnames(counts) <- c("S1", "S2")
    expect_error(TSENAT:::.ptt_prepare_inputs(counts = counts, readcounts = NULL, samples = c("S1", "S2"), coldata = NULL, sample_type_col = "sample_type", tx2gene = NULL, res = NULL, top_n = 2, pseudocount = 1e-6, output_file = NULL), "tx2gene")
})

test_that(".ptt_prepare_inputs returns list with mapping when provided", {
    skip_if_not_installed("ggplot2")
    counts <- matrix(1:6, nrow = 3)
    rownames(counts) <- paste0("tx", 1:3)
    colnames(counts) <- c("S1", "S2")
    tx2 <- data.frame(Transcript = rownames(counts), Gen = c("G1", "G1", "G2"), stringsAsFactors = FALSE)
    prep <- TSENAT:::.ptt_prepare_inputs(counts = counts, readcounts = NULL, samples = c("S1", "S2"), coldata = NULL, sample_type_col = "sample_type", tx2gene = tx2, res = NULL, top_n = 2, pseudocount = 1e-6, output_file = NULL)
    expect_type(prep, "list")
    expect_true(all(c("counts", "samples", "mapping", "agg_fun") %in% names(prep)))
})

test_that(".ptt_make_plot_for_gene returns ggplot object", {
    skip_if_not_installed("ggplot2")
    counts <- matrix(rpois(6, 10), nrow = 3)
    rownames(counts) <- paste0("tx", 1:3)
    colnames(counts) <- c("S1", "S2")
    mapping <- data.frame(Transcript = rownames(counts), Gen = c("G1", "G1", "G2"), stringsAsFactors = FALSE)
    agg_fun <- function(x) median(x, na.rm = TRUE)
    p <- TSENAT:::.ptt_make_plot_for_gene("G1", mapping = mapping, counts = counts, samples = c("Normal", "Tumor"), top_n = 2, agg_fun = agg_fun, pseudocount = 1e-6, agg_label_unique = "label")
    expect_s3_class(p, "ggplot")
})

test_that(".ptt_combine_plots returns a plot-like object", {
    skip_if_not_installed("ggplot2")
    p1 <- ggplot2::ggplot() +
        ggplot2::geom_point(mapping = ggplot2::aes(x = 1:3, y = 3:1))
    p2 <- ggplot2::ggplot() +
        ggplot2::geom_point(mapping = ggplot2::aes(x = 1:3, y = c(1, 2, 3)))
    out <- TSENAT:::.ptt_combine_plots(list(p1, p2), output_file = NULL, agg_label_unique = "agg")
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
    expect_error(TSENAT:::plot_top_transcripts(mat, gene = "G1", tx2gene = data.frame(Transcript = c("a", "b"), Gen = c("G1", "G1"))), "counts.*rownames")
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

skip_on_bioc()

context("Additional generate_plots tests")

library(TSENAT)

# plot_ma_expression errors when provided precomputed se missing 'log2_fold_change'

test_that("plot_ma_expression errors when precomputed se lacks log2_fold_change", {
    x <- data.frame(genes = paste0("g", seq_len(3)), mean = runif(3), stringsAsFactors = FALSE)
    bad_se <- data.frame(genes = paste0("g", seq_len(3)), other = rnorm(3), stringsAsFactors = FALSE)
    expect_error(plot_ma_expression(x, se = bad_se), "`se` data.frame must contain 'log2_fold_change'")
})

# plot_ma_expression accepts precomputed fold changes provided as data.frame/matrix with rownames

test_that("plot_ma_expression accepts precomputed fold changes with rownames", {
    x <- data.frame(genes = paste0("g", seq_len(4)), mean = runif(4), stringsAsFactors = FALSE)
    fc <- data.frame(log2_fold_change = rnorm(4))
    rownames(fc) <- paste0("g", seq_len(4))
    p <- plot_ma_expression(x, se = fc)
    expect_s3_class(p, "gg")
    expect_s3_class(p, "ggplot")
})

# plot_top_transcripts writes output_file for single and multiple genes

test_that("plot_top_transcripts writes output files for single and multiple genes", {
    skip_if_not_installed("ggplot2")
    set.seed(10)
    counts <- matrix(rpois(6 * 4, lambda = 10), nrow = 6)
    rownames(counts) <- paste0("tx", 1:6)
    colnames(counts) <- paste0("S", 1:4)
    samples <- c("N", "N", "T", "T")
    tx2 <- data.frame(Transcript = rownames(counts), Gen = rep(c("G1", "G2"), each = 3), stringsAsFactors = FALSE)

    tf1 <- tempfile(fileext = ".png")
    plot_top_transcripts(counts, gene = "G1", samples = samples, tx2gene = tx2, output_file = tf1)
    expect_true(file.exists(tf1) && file.info(tf1)$size > 0)

    tf2 <- tempfile(fileext = ".png")
    plot_top_transcripts(counts, gene = c("G1", "G2"), samples = samples, tx2gene = tx2, output_file = tf2)
    expect_true(file.exists(tf2) && file.info(tf2)$size > 0)
})

# plot_top_transcripts supports metric = 'iqr'

test_that("plot_top_transcripts supports metric 'iqr'", {
    skip_if_not_installed("ggplot2")
    counts <- matrix(rpois(3 * 4, lambda = 5), nrow = 3)
    rownames(counts) <- paste0("tx", 1:3)
    colnames(counts) <- paste0("S", 1:4)
    samples <- c("N", "N", "T", "T")
    tx2 <- data.frame(Transcript = rownames(counts), Gen = rep("G1", 3), stringsAsFactors = FALSE)

    p <- plot_top_transcripts(counts, gene = "G1", samples = samples, tx2gene = tx2, metric = "iqr")
    expect_s3_class(p, "ggplot")
})

# plot_volcano auto-detects x_col when not provided

test_that("plot_volcano auto-detects a numeric x column when x_col is NULL", {
    skip_if_not_installed("ggplot2")
    df <- data.frame(
        genes = paste0("g", seq_len(10)),
        stat = rnorm(10),
        adjusted_p_values = p.adjust(runif(10)),
        stringsAsFactors = FALSE
    )
    p <- plot_volcano(df, x_col = NULL, padj_col = "adjusted_p_values")
    expect_s3_class(p, "ggplot")
})

# .plot_ma_core uses fc_df values when provided; verify y values in plot data correspond to fc_df

test_that(".plot_ma_core uses fc_df values when provided", {
    skip_if_not_installed("ggplot2")
    x <- data.frame(genes = paste0("g", seq_len(5)), mean = runif(5), log2_fold_change = rnorm(5), stringsAsFactors = FALSE)
    fc <- data.frame(genes = x$genes, log2_fold_change = rnorm(5, mean = 5, sd = 0.1), stringsAsFactors = FALSE)

    p <- .plot_ma_core(x, fc_df = fc)
    expect_s3_class(p, "ggplot")
    pb <- ggplot2::ggplot_build(p)
    plotted_y <- pb$data[[1]]$y
    # Reconstruct expected merged df per implementation
    fdf <- as.data.frame(fc, stringsAsFactors = FALSE)
    df <- as.data.frame(x, stringsAsFactors = FALSE)
    df <- merge(df, fdf[, c("genes", "log2_fold_change")], by = "genes", all.x = TRUE, suffixes = c("", ".fc"))
    if ("log2_fold_change.fc" %in% colnames(df)) df$log2_fold_change <- ifelse(!is.na(df$log2_fold_change.fc), df$log2_fold_change.fc, df$log2_fold_change)
    expected_y <- as.numeric(df$log2_fold_change)
    expect_equal(plotted_y, expected_y)
})

context("plot_top_helpers extra cases")

library(testthat)

# .ptt_select_genes_from_res
test_that(".ptt_select_genes_from_res errors on NULL or missing genes", {
    expect_error(.ptt_select_genes_from_res(NULL, 3), "Either 'gene' or 'res' must be provided")
    expect_error(.ptt_select_genes_from_res(data.frame(a = 1:3), 2), "must contain a 'genes' column")
})

test_that(".ptt_select_genes_from_res sorts by adjusted or raw p-values and returns unique genes", {
    res <- data.frame(genes = c("g1", "g2", "g1", "g3"), adjusted_p_values = c(0.2, 0.01, 0.05, NA))
    sel <- .ptt_select_genes_from_res(res, top_n = 3)
    expect_true(all(c("g2", "g1") %in% sel))
    expect_length(unique(sel), length(sel))

    # raw p-values fallback
    res2 <- data.frame(genes = c("a", "b", "c"), raw_p_values = c(0.5, 0.1, 0.3))
    sel2 <- .ptt_select_genes_from_res(res2, top_n = 2)
    expect_equal(sel2, c("b", "c"))
})

# .ptt_infer_samples_from_coldata
test_that(".ptt_infer_samples_from_coldata handles row-named coldata and Sample id column", {
    counts <- matrix(1:6, nrow = 2)
    colnames(counts) <- c("S1", "S2", "S3")[1:ncol(counts)]

    cdf <- data.frame(sample_type = c("A", "B", "A"), stringsAsFactors = FALSE)
    rownames(cdf) <- c("S1", "S2", "S3")
    out <- .ptt_infer_samples_from_coldata(cdf, counts, "sample_type")
    expect_equal(out, c("A", "B", "A")[1:ncol(counts)])

    # use Sample id column
    cdf2 <- data.frame(Sample = c("S1", "S2", "S3"), sample_type = c("A", "B", "A"), stringsAsFactors = FALSE)
    out2 <- .ptt_infer_samples_from_coldata(cdf2, counts, "sample_type")
    expect_equal(out2, c("A", "B", "A")[1:ncol(counts)])

    # mismatched sample ids should error
    cdf_bad <- data.frame(Sample = c("X", "Y", "Z"), sample_type = c("A", "B", "A"), stringsAsFactors = FALSE)
    expect_error(.ptt_infer_samples_from_coldata(cdf_bad, counts, "sample_type"), "coldata sample id column does not match")
    expect_error(.ptt_infer_samples_from_coldata(123, counts, "sample_type"), "must be a data.frame or path")
})

# .ptt_read_tx2gene
test_that(".ptt_read_tx2gene validates inputs and reads mapping", {
    expect_error(.ptt_read_tx2gene(NULL), "`tx2gene` must be provided")

    bad <- data.frame(X = 1:2)
    expect_error(.ptt_read_tx2gene(bad), "tx2gene must have columns 'Transcript' and 'Gen'")

    good <- data.frame(Transcript = c("t1", "t2"), Gen = c("g1", "g1"), stringsAsFactors = FALSE)
    out <- .ptt_read_tx2gene(good)
    expect_equal(out, good)

    tf <- tempfile(fileext = ".tsv")
    write.table(good, file = tf, sep = "\t", row.names = FALSE, quote = FALSE)
    outf <- .ptt_read_tx2gene(tf)
    expect_true(is.data.frame(outf))
    unlink(tf)
    expect_error(.ptt_read_tx2gene("no_such_file.tsv"), "tx2gene file not found")
})

# .ptt_make_agg
test_that(".ptt_make_agg returns an aggregation function and label", {
    maj <- .ptt_make_agg("median")
    expect_equal(maj$metric_choice, "median")
    expect_true(is.function(maj$agg_fun))
    expect_true(grepl("median", maj$agg_label_unique, ignore.case = TRUE))

    m2 <- .ptt_make_agg("iqr")
    expect_equal(m2$metric_choice, "iqr")
    expect_true(grepl("IQR", m2$agg_label_unique))
})

# .ptt_build_tx_long and .ptt_aggregate_df_long
test_that(".ptt_build_tx_long and aggregation pipeline works and errors appropriately", {
    counts <- matrix(1:12, nrow = 4)
    rownames(counts) <- paste0("tx", 1:4)
    colnames(counts) <- paste0("S", 1:3)
    mapping <- data.frame(Transcript = paste0("tx", 1:4), Gen = c("G1", "G1", "G2", "G2"), stringsAsFactors = FALSE)
    samples <- c("A", "B", "A")

    expect_error(.ptt_build_tx_long("NOPE", mapping, counts, samples, top_n = NULL), "No transcripts found")

    res <- .ptt_build_tx_long("G1", mapping, counts, samples, top_n = 1)
    expect_true(is.list(res))
    expect_true(all(c("df_long", "txs") %in% names(res)))
    expect_equal(length(unique(res$df_long$tx)), length(res$txs))

    summ <- .ptt_aggregate_df_long(res$df_long, agg_fun = function(x) mean(x, na.rm = TRUE), pseudocount = 0.1)
    expect_true("log2expr" %in% colnames(summ))
    expect_true(is.factor(summ$tx))
})

# .ptt_build_plot_from_summary and combine functions
test_that(".ptt_build_plot_from_summary generates ggplot and combine functions operate", {
    skip_if_not_installed("ggplot2")
    p <- .ptt_build_plot_from_summary(data.frame(tx = factor(c("a", "b")), group = c("A", "B"), log2expr = c(1, 2)), "Label")
    expect_s3_class(p, "gg")

    # patchwork combine
    if (rlang::is_installed("patchwork")) {
        skip_if_not_installed("patchwork")
        p2 <- p + p
        combined <- .ptt_combine_patchwork(list(p, p), "Label")
        expect_true(inherits(combined, "patchwork"))
    }

    if (rlang::is_installed("cowplot")) {
        skip_if_not_installed("cowplot")
        outp <- .ptt_combine_cowplot(list(p, p), output_file = NULL, agg_label_unique = "Label")
        expect_true(inherits(outp, "gtable") || inherits(outp, "ggplot") || inherits(outp, "grob"))
    }

    if (rlang::is_installed("grid")) {
        skip_if_not_installed("grid")
        # .ptt_combine_grid returns invisibly NULL when not writing file and should not
        # create an Rplots.pdf in the working directory
        rpf <- "Rplots.pdf"
        if (file.exists(rpf)) unlink(rpf)
        res_grid <- .ptt_combine_grid(list(p, p), output_file = NULL, agg_label_unique = "Label")
        expect_null(res_grid)
        expect_false(file.exists(rpf))
    }
})

context("plot_helpers extras")

library(testthat)

test_that(".tsenat_format_label handles various inputs", {
    expect_null(.tsenat_format_label(NULL))
    expect_equal(.tsenat_format_label("__FOO_bar  "), "Foo bar")
    expect_equal(.tsenat_format_label(" a "), "A")
    expect_equal(.tsenat_format_label("   "), "")
    expect_equal(.tsenat_format_label("SINGLE"), "Single")
})

test_that(".tsenat_prepare_ma_plot_df handles mean_cols length >=2 and significance detection", {
    df <- data.frame(genes = c("g1", "g2", "g3"), meanA = c(1, 2, 3), meanB = c(1.5, 1.5, 1.5), log2fc = c(0, 1.2, -0.5), padj = c(0.2, 0.01, NA), stringsAsFactors = FALSE)
    res <- .tsenat_prepare_ma_plot_df(df, fold_col = "log2fc", mean_cols = c("meanA", "meanB"), x_label = NULL, y_label = "Log2FC")
    expect_is(res, "list")
    # when mean_cols length>=2 and x_label is NULL, default to 'meanA vs meanB'
    expect_equal(res$x_label, "meanA vs meanB")
    expect_true("plot_df" %in% names(res))
    expect_equal(nrow(res$plot_df), 3)
    # gene 2 should be significant (abs(y)>0 and padj<0.05)
    sig <- res$plot_df$significant
    expect_equal(sig, c("non-significant", "significant", "non-significant"))
})

test_that(".tsenat_prepare_ma_plot_df handles single mean col and fallback mean/index", {
    df1 <- data.frame(genes = c("g1", "g2"), m = c(5, 6), fc = c(0, 2), stringsAsFactors = FALSE)
    r1 <- .tsenat_prepare_ma_plot_df(df1, fold_col = "fc", mean_cols = c("m"), x_label = NULL, y_label = NULL)
    expect_equal(r1$x_label, "m")
    expect_equal(r1$plot_df$x, as.numeric(c(5, 6)))

    df2 <- data.frame(genes = c("g1", "g2"), mean = c(3, 4), fc = c(1, 0), stringsAsFactors = FALSE)
    r2 <- .tsenat_prepare_ma_plot_df(df2, fold_col = "fc", mean_cols = character(0), x_label = NULL, y_label = NULL)
    expect_equal(r2$x_label, "Mean")

    df3 <- data.frame(genes = c("g1", "g2"), fc = c(1, 2), stringsAsFactors = FALSE)
    r3 <- .tsenat_prepare_ma_plot_df(df3, fold_col = "fc", mean_cols = character(0), x_label = NULL, y_label = NULL)
    expect_equal(r3$x_label, "Index")
    expect_equal(r3$plot_df$x, c(1, 2))
})


test_that(".tsenat_prepare_volcano_df detects _difference column and formats labels", {
    df <- data.frame(gene = c("a", "b", "c"), median_difference = c(0.2, -0.5, 0.6), adjusted_p_values = c(0.2, 0.01, 0.001), stringsAsFactors = FALSE)
    res <- .tsenat_prepare_volcano_df(df)
    expect_equal(res$x_col, "median_difference")
    expect_equal(res$padj_col, "adjusted_p_values")
    expect_true("df" %in% names(res))
    expect_match(res$x_label_formatted, "Median")
    expect_match(res$padj_label_formatted, "Adjusted p values|Adjusted p values")
})

test_that(".tsenat_prepare_volcano_df errors for missing columns and empty data", {
    df <- data.frame(g = 1:3, something = letters[1:3], stringsAsFactors = FALSE)
    # Because 'g' is numeric it will be chosen as x_col but the default padj
    # column 'adjusted_p_values' is missing and an informative error is raised
    expect_error(.tsenat_prepare_volcano_df(df), "Column 'adjusted_p_values' not found")

    df2 <- data.frame(x = c(NA, Inf), adjusted_p_values = c(NA, NA), stringsAsFactors = FALSE)
    expect_error(.tsenat_prepare_volcano_df(df2, x_col = "x"), "No valid points to plot")

    df3 <- data.frame(x = c(1, 2), adj = c(0.01, 0.02), stringsAsFactors = FALSE)
    expect_error(.tsenat_prepare_volcano_df(df3, x_col = "x", padj_col = "nope"), "Column 'nope' not found")
})

test_that(".tsenat_prepare_volcano_df handles padj <=0 and signficance logic", {
    df <- data.frame(g = 1:4, value = c(0.2, 0.5, -0.2, 1), adjusted_p_values = c(0, 1e-10, 0.5, 0.001), stringsAsFactors = FALSE)
    res <- .tsenat_prepare_volcano_df(df, x_col = "value")
    expect_true(all(res$df$padj > 0))
    # label_thresh default 0.1: check significance assignment
    sig <- res$df$significant
    expect_equal(sig, ifelse(abs(res$df$xval) >= 0.1 & res$df$padj < 0.05, "significant", "non-significant"))
})

skip_on_bioc()

context("Unit tests for plot_top_transcripts helpers")

library(TSENAT)

# .ptt_select_genes_from_res

test_that(".ptt_select_genes_from_res selects by adjusted_p_values and raw_p_values", {
    res1 <- data.frame(genes = c("A", "B", "C"), adjusted_p_values = c(0.05, 0.01, 0.2), stringsAsFactors = FALSE)
    expect_equal(TSENAT:::.ptt_select_genes_from_res(res1, top_n = 2), c("B", "A"))

    res2 <- data.frame(genes = c("X", "Y", "Z"), raw_p_values = c(0.2, 0.01, 0.05), stringsAsFactors = FALSE)
    expect_equal(TSENAT:::.ptt_select_genes_from_res(res2, top_n = 2), c("Y", "Z"))

    expect_error(TSENAT:::.ptt_select_genes_from_res(NULL, top_n = 2))
    expect_error(TSENAT:::.ptt_select_genes_from_res(data.frame(a = 1), top_n = 2))
})

# .ptt_infer_samples_from_coldata

test_that(".ptt_infer_samples_from_coldata infers samples from data.frame and file path and errors on mismatch", {
    counts <- matrix(1:8, ncol = 4)
    colnames(counts) <- paste0("S", 1:4)

    cdf <- data.frame(sample_type = c("N", "T", "N", "T"), stringsAsFactors = FALSE)
    rownames(cdf) <- colnames(counts)

    samp <- TSENAT:::.ptt_infer_samples_from_coldata(cdf, counts, sample_type_col = "sample_type")
    expect_equal(as.character(samp), as.character(cdf[colnames(counts), "sample_type"]))

    # write as file with sample id column
    tf <- tempfile(fileext = ".tsv")
    dff <- data.frame(sample = colnames(counts), sample_type = c("N", "T", "N", "T"), stringsAsFactors = FALSE)
    utils::write.table(dff, file = tf, sep = "\t", quote = FALSE, row.names = FALSE)

    samp2 <- TSENAT:::.ptt_infer_samples_from_coldata(tf, counts, sample_type_col = "sample_type")
    expect_equal(as.character(samp2), as.character(dff$sample_type))

    # mismatch
    badcdf <- data.frame(other = c("a", "b"))
    expect_error(TSENAT:::.ptt_infer_samples_from_coldata(badcdf, counts, sample_type_col = "sample_type"))
})

# .ptt_read_tx2gene

test_that(".ptt_read_tx2gene reads mapping from data.frame and file and errors on missing columns", {
    mapping <- data.frame(Transcript = c("t1", "t2"), Gen = c("G1", "G1"), stringsAsFactors = FALSE)
    out <- TSENAT:::.ptt_read_tx2gene(mapping)
    expect_equal(out, mapping)

    tf <- tempfile(fileext = ".tsv")
    utils::write.table(mapping, file = tf, sep = "\t", quote = FALSE, row.names = FALSE)
    out2 <- TSENAT:::.ptt_read_tx2gene(tf)
    expect_equal(out2$Transcript, mapping$Transcript)

    expect_error(TSENAT:::.ptt_read_tx2gene(data.frame(a = 1)))
})

# .ptt_make_agg

test_that(".ptt_make_agg returns correct aggregator and label", {
    med <- TSENAT:::.ptt_make_agg("median")
    expect_equal(med$metric_choice, "median")
    expect_equal(med$agg_fun(c(1, 2, NA)), stats::median(c(1, 2, NA), na.rm = TRUE))

    mn <- TSENAT:::.ptt_make_agg("mean")
    expect_equal(mn$agg_fun(c(1, 2, NA)), mean(c(1, 2, NA), na.rm = TRUE))

    varr <- TSENAT:::.ptt_make_agg("variance")
    expect_equal(varr$agg_fun(c(1, 2, 3, NA)), stats::var(c(1, 2, 3, NA), na.rm = TRUE))

    iq <- TSENAT:::.ptt_make_agg("iqr")
    expect_equal(iq$agg_fun(c(1, 2, 3, 4, NA)), stats::IQR(c(1, 2, 3, 4, NA), na.rm = TRUE))

    # check counter side-effect increments
    opt_before <- as.integer(getOption("TSENAT.plot_top_counter", 0))
    TSENAT:::.ptt_make_agg("median")
    expect_true(as.integer(getOption("TSENAT.plot_top_counter", 0)) >= opt_before + 1)
})

# .ptt_build_tx_long & .ptt_aggregate_df_long & .ptt_build_plot_from_summary

test_that("tx long building, aggregation and plot building behave correctly", {
    counts <- matrix(rpois(6 * 2, lambda = 10), nrow = 6)
    rownames(counts) <- paste0("tx", 1:6)
    colnames(counts) <- paste0("S", 1:2)
    mapping <- data.frame(Transcript = rownames(counts), Gen = rep("G1", 6), stringsAsFactors = FALSE)
    samples <- c("N", "T")

    built <- TSENAT:::.ptt_build_tx_long("G1", mapping, counts, samples, top_n = 3)
    expect_true(is.list(built))
    expect_true(all(c("df_long", "txs") %in% names(built)))
    expect_true(length(built$txs) <= 3)
    expect_true(all(c("tx", "sample", "expr", "group") %in% colnames(built$df_long)))

    df_summary <- TSENAT:::.ptt_aggregate_df_long(built$df_long, agg_fun = function(x) mean(x, na.rm = TRUE), pseudocount = 1e-6)
    expect_true(all(c("tx", "group", "expr", "log2expr") %in% colnames(df_summary)))
    expect_true(is.factor(df_summary$tx))

    skip_if_not_installed("ggplot2")
    p <- TSENAT:::.ptt_build_plot_from_summary(df_summary, agg_label_unique = "label")
    expect_s3_class(p, "gg")
})

# .ptt_combine_grid writes to file when output_file provided

test_that(".ptt_combine_grid writes a PNG file when output_file is given", {
    skip_if_not_installed("ggplot2")
    library(ggplot2)

    df <- data.frame(x = 1:3, y = rnorm(3))
    p1 <- ggplot(df, ggplot2::aes(x = x, y = y)) +
        ggplot2::geom_point()
    p2 <- ggplot(df, ggplot2::aes(x = x, y = -y)) +
        ggplot2::geom_point()

    tf <- tempfile(fileext = ".png")
    # call grid combiner directly
    TSENAT:::.ptt_combine_grid(list(p1, p2), output_file = tf, agg_label_unique = "agg")
    expect_true(file.exists(tf))
    expect_true(file.info(tf)$size > 0)
})

context("plot_tsallis_gene_profile extra cases")

library(testthat)

# vector-of-genes input returns a named list of ggplots
test_that("plot_tsallis_gene_profile accepts vector of genes and returns list of ggplots", {
    skip_if_not_installed("SummarizedExperiment")
    skip_if_not_installed("ggplot2")

    set.seed(101)
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

    plots <- plot_tsallis_gene_profile(ts_se, gene = c("G1", "G2"))
    expect_type(plots, "list")
    expect_equal(length(plots), 2)
    expect_true(all(c("G1", "G2") %in% names(plots)))
    lapply(plots, function(p) expect_s3_class(p, "ggplot"))
})

# lm_res supplied with gene = NULL picks top n_top genes
test_that("plot_tsallis_gene_profile uses lm_res when gene is NULL and returns up to n_top plots", {
    skip_if_not_installed("SummarizedExperiment")
    skip_if_not_installed("ggplot2")

    set.seed(202)
    # Build a SummarizedExperiment with a controlled interaction for first few genes
    base_samples <- paste0("S", 1:6)
    qvals <- c(0.1, 0.5)
    # column names like S1_q=0.1, S1_q=0.5, S2_q=0.1, S2_q=0.5, ...
    sample_cols <- as.character(t(outer(base_samples, qvals, function(s, q) paste0(s, "_q=", q))))

    n_genes <- 50
    mat <- matrix(rnorm(n_genes * length(sample_cols), sd = 0.1), nrow = n_genes, ncol = length(sample_cols))
    rownames(mat) <- paste0("G", seq_len(n_genes))
    colnames(mat) <- sample_cols

    # create a strong q:group interaction for first 5 genes: tumor samples have slope with q
    conds <- rep(c("Normal", "Tumor"), length.out = length(base_samples))
    for (g in 1:5) {
        for (j in seq_along(sample_cols)) {
            base_idx <- ((j - 1) %/% length(qvals)) + 1
            qv <- qvals[((j - 1) %% length(qvals)) + 1]
            cond <- conds[base_idx]
            mat[g, j] <- 0.2 + (cond == "Tumor") * (1.0 * qv) + rnorm(1, 0, 0.05)
        }
    }

    se <- SummarizedExperiment::SummarizedExperiment(assays = list(diversity = mat))
    coldata_df <- data.frame(
        Sample = base_samples,
        Condition = conds,
        stringsAsFactors = FALSE
    )
    ts_se <- map_metadata(se, coldata_df)

    # Run linear lm interaction (fast) to get lm_res
    lm_res <- calculate_lm_interaction(ts_se, sample_type_col = "sample_type", method = "linear", pvalue = "lrt", min_obs = 2)
    if (nrow(lm_res) == 0) {
        # fallback synthetic lm_res: ensure top genes exist for plotting
        pvals <- runif(50, min = 0.01, max = 1)
        pvals[1:5] <- runif(5, min = 0, max = 1e-6)
        lm_res <- data.frame(gene = paste0("G", seq_len(50)), p_interaction = pvals, adj_p_interaction = stats::p.adjust(pvals, method = "BH"), stringsAsFactors = FALSE)
    }
    # request top 5 instead of default 10 for speed/stability
    plots <- plot_tsallis_gene_profile(ts_se, gene = NULL, lm_res = lm_res, n_top = 5)
    expect_type(plots, "list")
    expect_lte(length(plots), 5)
    # each element must be ggplot
    lapply(plots, function(p) expect_s3_class(p, "ggplot"))
})

context("generate_plots.R coverage")

test_that("require_pkgs errors if packages are not installed", {
    # This test will fail if the package is actually installed. Use a highly
    # improbable package name to avoid needing to mock `requireNamespace`.
    expect_error(TSENAT:::require_pkgs("definitely_not_installed_pkg_12345"), "definitely_not_installed_pkg_12345 required")
})

test_that("infer_samples_from_se fallback logic works", {
    se <- SummarizedExperiment(
        assays = list(counts = matrix(1:8, ncol = 4)),
        colData = DataFrame(foo = c("a", "b", "c", "d"), bar = c("x", "y", "z", "w"))
    )
    # It should pick 'foo' as it has fewer unique values > 1
    expect_equal(TSENAT:::infer_samples_from_se(se), c("a", "b", "c", "d"))

    se2 <- SummarizedExperiment(
        assays = list(counts = matrix(1:4, nrow = 2)),
        colData = DataFrame(baz = c(1, 1), qux = c("a", "b"))
    )
    expect_equal(TSENAT:::infer_samples_from_se(se2), c("a", "b"))

    se3 <- SummarizedExperiment(assays = list(counts = matrix(1:4, nrow = 2)))
    colData(se3) <- NULL
    expect_null(TSENAT:::infer_samples_from_se(se3))
})

test_that("get_readcounts_from_se works with file path and fallback", {
    # test with a file path (with gene column)
    rc_df <- data.frame(gene = c("g1", "g2"), c1 = c(1, 2), c2 = c(3, 4))
    rc_file <- tempfile()
    write.table(rc_df, rc_file, sep = "\t", row.names = FALSE)
    se <- SummarizedExperiment()
    rc <- TSENAT:::get_readcounts_from_se(se, readcounts_arg = rc_file)
    expect_equal(nrow(rc), 2)
    expect_equal(ncol(rc), 2)

    # test with a single-column file (no gene column) -> returns matrix of values
    rc_single <- data.frame(V1 = c(1, 2))
    rc_file_single <- tempfile()
    # include a header so read.delim(..., header = TRUE) reads two rows
    write.table(rc_single, rc_file_single, sep = "\t", row.names = FALSE, col.names = TRUE)
    rc2 <- TSENAT:::get_readcounts_from_se(se, readcounts_arg = rc_file_single)
    expect_equal(as.vector(rc2), c(1, 2))

    # test with a data.frame
    rc_df_no_gene <- data.frame(c1 = c(1, 2), c2 = c(3, 4))
    rc <- TSENAT:::get_readcounts_from_se(se, readcounts_arg = rc_df_no_gene)
    expect_equal(nrow(rc), 2)
    expect_equal(ncol(rc), 2)

    # test with a matrix with no rownames
    rc_mat_no_rownames <- matrix(1:4, 2)
    rc <- TSENAT:::get_readcounts_from_se(se, readcounts_arg = rc_mat_no_rownames)
    expect_equal(nrow(rc), 2)

    # test error for invalid readcounts_arg
    expect_error(TSENAT:::get_readcounts_from_se(se, readcounts_arg = 123), "`readcounts` must be a matrix/data.frame or path to a file")

    # test fallback to first assay
    se_assay <- SummarizedExperiment(assays = list(my_counts = matrix(1:4, 2)))
    expect_warning(rc_assay <- TSENAT:::get_readcounts_from_se(se_assay), "Using first assay from SummarizedExperiment")
    expect_equal(nrow(rc_assay), 2)

    # metadata$readcounts should be preferred when present
    se_meta <- SummarizedExperiment(assays = list(my_counts = matrix(1:4, nrow = 2)))
    S4Vectors::metadata(se_meta) <- list(readcounts = matrix(5:8, nrow = 2, dimnames = list(c("tx1", "tx2"), c("s1", "s2"))))
    rc_meta <- TSENAT:::get_readcounts_from_se(se_meta)
    expect_equal(as.vector(rc_meta), c(5, 6, 7, 8))

    # preferred assay name 'readcounts' should be selected when present
    se_pref <- SummarizedExperiment(assays = list(readcounts = matrix(11:14, nrow = 2, dimnames = list(c("tx1", "tx2"), c("s1", "s2"))), counts = matrix(1:4, nrow = 2)))
    rc_pref <- TSENAT:::get_readcounts_from_se(se_pref)
    expect_equal(as.vector(rc_pref), c(11, 12, 13, 14))
})

test_that("get_tx2gene_from_se fallback works", {
    se <- SummarizedExperiment(assays = list(counts = matrix(1:4, nrow = 2, dimnames = list(c("tx1", "tx2"), c("s1", "s2")))))
    res <- TSENAT:::get_tx2gene_from_se(se, readcounts_mat = assay(se))
    expect_equal(res$mapping, c("tx1", "tx2"))

    # test with different column names
    md <- list(tx2gene = data.frame(tx = c("tx1", "tx2"), g = c("g1", "g1")))
    S4Vectors::metadata(se) <- md
    res2 <- TSENAT:::get_tx2gene_from_se(se, readcounts_mat = assay(se))
    expect_equal(res2$mapping, c("g1", "g1"))
})

test_that("plot_tsallis_gene_profile handles errors", {
    mat <- matrix(1:4, nrow = 2, dimnames = list(NULL, c("s1", "s2")))
    se <- SummarizedExperiment(assays = list(diversity = mat))
    colData(se) <- DataFrame(sample_type = c("a", "b"), row.names = c("s1", "s2"))
    expect_error(TSENAT::plot_tsallis_gene_profile(se, gene = "g1"), "Gene not found in assay: g1")

    # Create a SE where prepare_tsallis_long will return a Gene 'g1' so
    # plotting a different gene errors as expected.
    mat2 <- matrix(1, nrow = 1, ncol = 1, dimnames = list(NULL, c("s1_q=0.5")))
    se2 <- SummarizedExperiment(assays = list(diversity = mat2))
    SummarizedExperiment::rowData(se2)$genes <- "g1"
    SummarizedExperiment::colData(se2) <- DataFrame(sample_type = "A", row.names = "s1")
    expect_error(TSENAT::plot_tsallis_gene_profile(se2, gene = "g2"), "Gene not found in assay: g2")

    # gene is NULL, lm_res is not a data.frame
    expect_error(TSENAT::plot_tsallis_gene_profile(se, lm_res = "not a dataframe"), "'lm_res' must be a data.frame with a 'gene' column")

    # lm_res is a data.frame but without a 'gene' column
    lm_res_no_gene <- data.frame(p = c(0.1, 0.2))
    expect_error(TSENAT::plot_tsallis_gene_profile(se, lm_res = lm_res_no_gene), "'lm_res' must be a data.frame with a 'gene' column")

    # lm_res is a data.frame with 'gene' column but no p-value column
    lm_res_no_p <- data.frame(gene = c("g1", "g2"))
    expect_error(TSENAT::plot_tsallis_gene_profile(se, lm_res = lm_res_no_p), "'lm_res' must contain 'adj_p_interaction' or 'p_interaction' columns")
})

# Additional tests to exercise plotting branches
test_that("plot_tsallis_gene_profile plotting branches", {
    mat_full <- matrix(c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8), nrow = 2)
    colnames(mat_full) <- c("s1_q=0.1", "s2_q=0.1", "s1_q=1", "s2_q=1")
    rownames(mat_full) <- c("g1", "g2")
    se_full <- SummarizedExperiment(assays = list(diversity = mat_full))
    SummarizedExperiment::rowData(se_full)$genes <- rownames(mat_full)
    # colData must have one row per column in the assay matrix
    SummarizedExperiment::colData(se_full) <- DataFrame(sample_type = c("A", "B", "A", "B"))

    # Single gene plot with sample lines shown (skip sample_type mapping)
    p_single <- TSENAT::plot_tsallis_gene_profile(se_full, gene = "g1", show_samples = TRUE, sample_type_col = NULL)
    expect_s3_class(p_single, "ggplot")
    has_line <- any(vapply(p_single$layers, function(l) {
        if (!is.null(l$geom) && is.character(l$geom$geom_name)) {
            l$geom$geom_name == "line"
        } else {
            inherits(l$geom, "GeomLine")
        }
    }, logical(1)))
    expect_true(has_line)

    # Multi-gene via lm_res using adj_p_interaction
    lm_res <- data.frame(gene = c("g2", "g1", "g3"), adj_p_interaction = c(0.01, 0.1, 0.2))
    plots_list <- TSENAT::plot_tsallis_gene_profile(se_full, gene = NULL, lm_res = lm_res, n_top = 2, sample_type_col = NULL)
    expect_true(is.list(plots_list))
    expect_equal(names(plots_list), c("g2", "g1"))

    # lm_res with only p_interaction (no adj) should also work
    lm_res2 <- data.frame(gene = c("g1", "g2"), p_interaction = c(0.05, 0.01))
    plots_list2 <- TSENAT::plot_tsallis_gene_profile(se_full, gene = NULL, lm_res = lm_res2, n_top = 2, sample_type_col = NULL)
    expect_true(is.list(plots_list2))
    expect_equal(names(plots_list2), c("g2", "g1"))
})

test_that(".plot_ma_core handles more edge cases", {
    # No genes column, but rownames are present
    df <- data.frame(mean = runif(5), log2_fold_change = rnorm(5))
    rownames(df) <- paste0("g", 1:5)
    p <- TSENAT:::.plot_ma_core(df)
    expect_s3_class(p, "ggplot")

    # fc_df without 'log2_fold_change' column
    fc_df_bad <- data.frame(genes = paste0("g", 1:5))
    expect_error(TSENAT:::.plot_ma_core(df, fc_df = fc_df_bad), "Provided `fc_df` must contain 'log2_fold_change' column")

    # y_label_formatted branch
    p2 <- TSENAT:::.plot_ma_core(df, y_label = "log2")
    expect_s3_class(p2, "ggplot")
})

test_that("plot_ma_expression_impl handles errors", {
    x <- data.frame(genes = paste0("g", 1:5), mean = runif(5))
    # create a SummarizedExperiment with an assay but no colData so infer_samples_from_se returns NULL
    se <- SummarizedExperiment(assays = list(counts = matrix(1:4, nrow = 2, dimnames = list(c("tx1", "tx2"), c("s1", "s2")))))
    expect_error(TSENAT:::plot_ma_expression_impl(x, se = se), "Could not infer 'samples' from SummarizedExperiment; provide `samples`")

    # se as a matrix without log2_fold_change
    fc_mat_bad <- matrix(1:5)
    expect_error(TSENAT:::plot_ma_expression_impl(x, se = fc_mat_bad), "`se` data.frame must contain 'log2_fold_change' column")

    # Unsupported 'se' argument
    expect_error(TSENAT:::plot_ma_expression_impl(x, se = 123), "Unsupported 'se' argument")
})

test_that("plot_tsallis_q_curve handles single group and empty long df", {
    se <- SummarizedExperiment(assays = list(diversity = matrix(rnorm(4), 2, dimnames = list(NULL, c("s1_q=0.1", "s2_q=0.1")))))
    colData(se) <- DataFrame(sample_type = c("A", "A"), row.names = c("s1", "s2"))
    p <- plot_tsallis_q_curve(se)
    expect_s3_class(p, "ggplot")
    # check that legend is removed for single group
    expect_true(p$theme$legend.position == "none")

    # empty long df: create a SE with only NA tsallis values so the
    # helper returns no valid rows
    se_empty <- SummarizedExperiment(assays = list(diversity = matrix(NA_real_, nrow = 1, ncol = 1, dimnames = list(NULL, c("s1_q=0.1")))))
    SummarizedExperiment::rowData(se_empty)$genes <- "g1"
    SummarizedExperiment::colData(se_empty) <- DataFrame(sample_type = "A", row.names = "s1")
    expect_error(plot_tsallis_q_curve(se_empty), "No tsallis values found in SummarizedExperiment")

    # not a summarized experiment
    expect_error(plot_tsallis_q_curve(123), "requires a SummarizedExperiment")
})

test_that(".compute_transcript_fill_limits handles no transcripts found", {
    counts <- matrix(1:4, 2)
    rownames(counts) <- c("tx1", "tx2")
    mapping <- data.frame(Transcript = c("tx3"), Gen = c("g1"))
    samples <- c("a", "b")
    expect_error(TSENAT:::.compute_transcript_fill_limits(genes = "g1", mapping = mapping, counts = counts, samples = samples, top_n = 1, agg_fun = mean, pseudocount = 1), "No transcripts found for provided genes")

    # case where one gene has no txs, but other does
    mapping2 <- data.frame(Transcript = c("tx1", "tx4"), Gen = c("g2", "g3"))
    limits <- TSENAT:::.compute_transcript_fill_limits(genes = c("g1", "g2"), mapping = mapping2, counts = counts, samples = samples, top_n = 1, agg_fun = mean, pseudocount = 1)
    expect_is(limits, "numeric")
})

test_that(".draw_transcript_grid creates a temporary pdf in non-interactive sessions", {
    # This is hard to test directly, but we can check the logic.
    # We can't easily force a non-interactive session in a test.
    # We can check that it doesn't error when no device is open.
    grob <- grid::rectGrob()
    expect_silent(TSENAT:::.draw_transcript_grid(list(grob), "title", NULL, 1, grid::unit(1, "null")))

    # test with file (open a device so the function can close it)
    tf <- tempfile(fileext = ".png")
    png(tf, width = 400, height = 300)
    expect_silent(TSENAT:::.draw_transcript_grid(list(grob), "title", NULL, 1, grid::unit(1, "null"), to_file = tf))
    expect_true(file.exists(tf))
    if (file.exists(tf)) unlink(tf)
})

test_that("plot_volcano handles errors", {
    df <- data.frame(gene = c("a", "b"), p = c(0.1, 0.01))
    expect_error(plot_volcano(df), "Column 'adjusted_p_values' not found in diff_df")
})

test_that(".ptt_combine_plots fallbacks work", {
    p1 <- ggplot2::ggplot()
    # Ensure the function runs and falls back to any available backend;
    # don't rely on mocking namespace checks here.
    expect_silent(.ptt_combine_plots(list(p1), "label"))
})

test_that(".ptt_combine_plots treats single string second arg as label", {
    p1 <- ggplot2::ggplot()
    # explicit label
    out1 <- NULL
    out2 <- NULL
    expect_error(out1 <- .ptt_combine_plots(list(p1), output_file = NULL, agg_label_unique = "mylabel"), NA)
    expect_error(out2 <- .ptt_combine_plots(list(p1), "mylabel"), NA)
    expect_equal(class(out1), class(out2))
})

test_that(".ptt_prepare_inputs handles file paths and various errors", {
    counts <- matrix(1:4, 2)
    rownames(counts) <- paste0("tx", 1:2)
    colnames(counts) <- paste0("s", 1:2)
    samples <- c("a", "b")

    # coldata as file
    cd_file <- tempfile()
    write.table(data.frame(sample_id = c("s1", "s2"), sample_type = c("a", "b")), cd_file, sep = "\t", row.names = F)

    # tx2gene as file
    t2g_file <- tempfile()
    write.table(data.frame(Transcript = c("tx1", "tx2"), Gen = c("g1", "g1")), t2g_file, sep = "\t", row.names = F)

    prep <- .ptt_prepare_inputs(counts, samples = NULL, coldata = cd_file, sample_type_col = "sample_type", tx2gene = t2g_file, res = NULL, top_n = 1, pseudocount = 1)
    expect_equal(prep$samples, c("a", "b"))

    # bad coldata (no sample_id-like columns)
    bad_cd_file <- tempfile()
    write.table(data.frame(x = 1), bad_cd_file)
    expect_error(.ptt_prepare_inputs(counts, samples = NULL, coldata = bad_cd_file, tx2gene = t2g_file), "Could not match `coldata` rows to `counts` columns")

    # coldata sample_id column does not match counts column names
    cd_file_mismatch <- tempfile()
    write.table(data.frame(sample_id = c("s3", "s4"), sample_type = c("a", "b")), cd_file_mismatch, sep = "\t", row.names = FALSE)
    expect_error(.ptt_prepare_inputs(counts, samples = NULL, coldata = cd_file_mismatch, tx2gene = t2g_file), "coldata sample id column does not match column names of counts")

    # coldata file path not found
    expect_error(.ptt_prepare_inputs(counts, samples = NULL, coldata = "no_such_file.tsv", tx2gene = t2g_file), "coldata file not found")

    # tx2gene must be provided
    expect_error(.ptt_prepare_inputs(counts, samples = samples, tx2gene = NULL), "`tx2gene` must be provided")

    # tx2gene file not found
    expect_error(.ptt_prepare_inputs(counts, samples = samples, tx2gene = "no_such_tx2gene.tsv"), "tx2gene file not found")

    # tx2gene missing required columns
    bad_t2g <- tempfile()
    write.table(data.frame(A = 1, B = 2), bad_t2g, sep = "\t", row.names = FALSE)
    expect_error(.ptt_prepare_inputs(counts, samples = samples, tx2gene = bad_t2g), "tx2gene must have columns 'Transcript' and 'Gen'")

    # counts must have rownames
    counts_no_rownames <- matrix(1:4, 2)
    expect_error(.ptt_prepare_inputs(counts_no_rownames, samples = samples, tx2gene = t2g_file), "`counts` must have rownames corresponding to transcript identifiers")

    # counts must be matrix/data.frame
    expect_error(.ptt_prepare_inputs(123, samples = samples, tx2gene = t2g_file), "`counts` must be a matrix or data.frame")

    # samples length must match number of columns
    expect_error(.ptt_prepare_inputs(counts, samples = c("a"), tx2gene = t2g_file), "Length of `samples` must equal number of columns in `counts`")

    # SummarizedExperiment input with tx2gene in metadata
    se <- SummarizedExperiment(assays = list(counts = matrix(1:4, nrow = 2, dimnames = list(c("tx1", "tx2"), c("s1", "s2")))))
    S4Vectors::metadata(se) <- list(tx2gene = data.frame(Transcript = c("tx1", "tx2"), Gen = c("g1", "g1"), stringsAsFactors = FALSE))
    cd_file2 <- tempfile()
    write.table(data.frame(sample_id = c("s1", "s2"), sample_type = c("a", "b")), cd_file2, sep = "\t", row.names = FALSE)
    prep2 <- .ptt_prepare_inputs(se, readcounts = NULL, samples = NULL, coldata = cd_file2, sample_type_col = "sample_type", tx2gene = NULL, res = NULL, top_n = 1, pseudocount = 1, output_file = NULL)
    expect_equal(prep2$mapping$Gen, c("g1", "g1"))

    # no samples or coldata
    expect_error(.ptt_prepare_inputs(counts, tx2gene = t2g_file), "Either 'samples' or 'coldata' must be provided")
})

test_that("plot_diversity_density handles empty data gracefully", {
    # Bug fix: plot_diversity_density failed with cryptic ggplot2 error
    # when get_assay_long returned empty dataframe (all values NA)
    # Solution: Validate data before plotting and provide clear error message
    
    skip_if_not_installed("SummarizedExperiment")
    skip_if_not_installed("ggplot2")
    
    library(SummarizedExperiment)
    
    # Create SE with all NA values
    mat <- matrix(NA_real_, nrow = 5, ncol = 4)
    colnames(mat) <- c("S1_q=0.1", "S2_q=0.1", "S3_q=0.1", "S4_q=0.1")
    rownames(mat) <- paste0("G", 1:5)
    rowData_df <- S4Vectors::DataFrame(genes = rownames(mat))
    colData_df <- S4Vectors::DataFrame(samples = c("S1", "S2", "S3", "S4"))
    se <- SummarizedExperiment(
        assays = list(diversity = mat),
        rowData = rowData_df,
        colData = colData_df
    )
    
    # Should raise clear error about no non-NA values, not ggplot2 error
    expect_error(plot_diversity_density(se),
                 "No non-NA values found",
                 info = "Should have clear error message about NA values")
})

test_that("plot_diversity_density works with small dataset from calculate_diversity", {
    # Regression test: Ensure plot_diversity_density works when diversity values are calculated
    # Use unnormalized entropy to ensure all genes (including single-isoform) produce finite values
    skip_if_not_installed("SummarizedExperiment")
    skip_if_not_installed("ggplot2")
    
    library(SummarizedExperiment)
    library(ggplot2)
    
    data("tcga_brca_luma", package = "TSENAT")
    rc <- as.matrix(tcga_brca_luma[1:100, -1, drop = FALSE])
    gs <- tcga_brca_luma[1:100, 1]
    # Use norm=FALSE to ensure single-isoform genes produce finite values
    se <- calculate_diversity(rc, gs, q = 0.1, norm = FALSE)
    
    # Should succeed
    p <- plot_diversity_density(se)
    expect_s3_class(p, "ggplot")
})

