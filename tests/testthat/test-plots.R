# Comprehensive testing of all plotting functions
# Tests plot_ma, plot_top_transcripts, plot_volcano, plot_tsallis_q_curve,
# plot_tsallis_violin_multq

library(TSENAT)
skip_on_bioc()

context("plots: Visualization and Data Exploration")

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

    se <- SummarizedExperiment(
        assays = list(counts = counts),
        rowData = S4Vectors::DataFrame(genes = tx2$Gen),
        colData = S4Vectors::DataFrame(sample_type = samples)
    )
    p <- plot_top_transcripts(se,
        gene = "GENE1",
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

    se <- SummarizedExperiment(
        assays = list(counts = counts),
        rowData = S4Vectors::DataFrame(genes = tx2$Gen),
        colData = S4Vectors::DataFrame(sample_type = samples)
    )
    p <- plot_top_transcripts(se, res = res, top_n = 2)
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

    ts_se <- TSENAT:::.map_metadata(ts_se, coldata_df)

    p <- plot_tsallis_q_curve(ts_se)
    expect_true(inherits(p, "ggplot"))
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

test_that(".plot_ma_core errors when fold-change column missing or x axis missing", {
    skip_if_not_installed("ggplot2")
    df <- data.frame(genes = paste0("g", 1:4), val = runif(4))
    expect_error(.plot_ma_core(df), "Could not find a fold-change column")
    df2 <- data.frame(genes = paste0("g", 1:4), log2_fold_change = rnorm(4))
    expect_s3_class(.plot_ma_core(df2), "ggplot")
})

context("Visualization: Top Transcripts Plotting")

library(SummarizedExperiment)

test_that("plot_top_transcripts works on simple matrix input", {
    tx_counts <- matrix(sample(1:100, 24, replace = TRUE), nrow = 6)
    rownames(tx_counts) <- paste0("tx", seq_len(nrow(tx_counts)))
    colnames(tx_counts) <- paste0("S", seq_len(ncol(tx_counts)))

    tx2gene <- data.frame(Transcript = rownames(tx_counts), Gen = rep(paste0("G", seq_len(3)), each = 2), stringsAsFactors = FALSE)
    samples <- rep(c("Normal", "Tumor"), length.out = ncol(tx_counts))

    se <- SummarizedExperiment(
        assays = list(counts = tx_counts),
        rowData = S4Vectors::DataFrame(genes = tx2gene$Gen),
        colData = S4Vectors::DataFrame(sample_type = samples)
    )
    p <- plot_top_transcripts(se, gene = c("G1", "G2"), top_n = 2)
    expect_true(!is.null(p))
    # expect ggplot object or patchwork
    expect_true(inherits(p, "ggplot") || inherits(p, "patchwork") || inherits(p, "gtable") || inherits(p, "ggarrange"))
})

test_that("plot_top_transcripts errors when se is not SummarizedExperiment", {
    mat <- matrix(1:6, nrow = 2)
    expect_error(plot_top_transcripts(mat, gene = "G1"), "se must be a SummarizedExperiment")
})

context("Visualization: Generate Plots Additional Tests")

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

context("Visualization: Generate Plots Extra Tests")

library(SummarizedExperiment)

skip_on_bioc()

test_that("plot_ma_tsallis handles simple inputs", {
    skip_if_not_installed("ggplot2")
    x <- data.frame(genes = paste0("g", 1:6), mean = runif(6), log2_fold_change = rnorm(6))
    p1 <- TSENAT::plot_ma_tsallis(x)
    expect_s3_class(p1, "ggplot")
})


test_that("plot_tsallis_q_curve correctly handles multiple groups with different entropy values", {
    skip_if_not_installed(c("ggplot2", "SummarizedExperiment", "dplyr"))
    library(dplyr)
    
    # Create SE with two groups having different entropy profiles
    set.seed(42)
    n_genes <- 20
    n_q_vals <- 5
    n_samples_per_group <- 4
    
    # Create data where "normal" group has higher entropy than "tumor" group across all q-values
    normal_data <- matrix(rnorm(n_genes * n_q_vals * n_samples_per_group, mean = 0.7, sd = 0.1), 
                          nrow = n_genes)
    tumor_data <- matrix(rnorm(n_genes * n_q_vals * n_samples_per_group, mean = 0.4, sd = 0.1), 
                         nrow = n_genes)
    
    mat <- cbind(normal_data, tumor_data)
    
    # Create column names with multiple q-values
    q_vals <- seq(0.1, 0.5, by = 0.1)
    col_names <- c(
        paste0("S", 1:n_samples_per_group, "_q=", rep(q_vals, each = n_samples_per_group)),
        paste0("S", (n_samples_per_group+1):(2*n_samples_per_group), "_q=", rep(q_vals, each = n_samples_per_group))
    )
    
    colnames(mat) <- col_names
    rownames(mat) <- paste0("g", 1:n_genes)
    
    # Ensure matrix values are in [0, 1]
    mat <- pmax(pmin(mat, 1), 0)
    
    se <- SummarizedExperiment::SummarizedExperiment(assays = list(diversity = mat))
    rowData(se)$genes <- rownames(mat)
    
    # Set sample type
    sample_types <- c(rep("normal", n_samples_per_group * n_q_vals), 
                      rep("tumor", n_samples_per_group * n_q_vals))
    cd <- S4Vectors::DataFrame(sample_type = sample_types, row.names = colnames(mat))
    SummarizedExperiment::colData(se) <- cd
    
    # Generate plot
    p <- TSENAT:::plot_tsallis_q_curve(se, sample_type_col = "sample_type")
    
    # Verify plot is ggplot
    expect_s3_class(p, "ggplot")
    
    # Verify plot data contains both groups
    plot_data <- p$data
    expect_true("group" %in% colnames(plot_data))
    expect_true("normal" %in% plot_data$group)
    expect_true("tumor" %in% plot_data$group)
    
    # Verify q values are numeric and correct
    expect_true("q" %in% colnames(plot_data))
    expect_true(is.numeric(plot_data$q))
    expect_false(is.factor(plot_data$q))
    
    # Verify there are multiple q-values
    unique_q_vals <- unique(plot_data$q)
    expect_equal(length(unique_q_vals), length(q_vals))
    
    # Verify normal group has higher median entropy than tumor group (based on our data construction)
    normal_medians <- filter(plot_data, group == "normal") %>% pull(median)
    tumor_medians <- filter(plot_data, group == "tumor") %>% pull(median)
    expect_true(mean(normal_medians) > mean(tumor_medians))
})

test_that("plot_tsallis_q_curve preserves decimal q-values correctly", {
    skip_if_not_installed(c("ggplot2", "SummarizedExperiment", "dplyr"))
    library(dplyr)
    
    # Create SE with decimal q-values
    q_decimal_vals <- c(0.15, 0.35)
    n_samples <- 3
    n_genes <- 2
    n_cols <- n_samples * length(q_decimal_vals)
    
    mat <- matrix(rnorm(n_genes * n_cols, mean = 0.5, sd = 0.1), nrow = n_genes, ncol = n_cols)
    col_names <- character(n_cols)
    idx <- 1
    for (q in q_decimal_vals) {
        for (s in 1:n_samples) {
            col_names[idx] <- paste0("S", s, "_q=", q)
            idx <- idx + 1
        }
    }
    colnames(mat) <- col_names
    rownames(mat) <- c("g1", "g2")
    
    se <- SummarizedExperiment::SummarizedExperiment(assays = list(diversity = mat))
    rowData(se)$genes <- rownames(mat)
    
    # Set sample type
    cd <- S4Vectors::DataFrame(sample_type = rep(c("N", "T", "N"), length(q_decimal_vals)), row.names = colnames(mat))
    SummarizedExperiment::colData(se) <- cd
    
    # Generate plot
    p <- TSENAT:::plot_tsallis_q_curve(se, sample_type_col = "sample_type")
    
    # Verify plot data q values are numeric
    plot_data <- p$data
    expect_true(is.numeric(plot_data$q))
    
    # Verify q-values are preserved (should be approximately equal, accounting for floating point)
    plotted_q <- sort(unique(plot_data$q))
    expected_q <- sort(unique(q_decimal_vals))
    expect_equal(length(plotted_q), length(expected_q))
    
    # Check each q-value with tolerance for floating point
    for (i in seq_len(length(expected_q))) {
        expect_true(abs(plotted_q[i] - expected_q[i]) < 1e-10)
    }
})

test_that("plot_volcano auto-detects x_col and returns ggplot", {
    skip_if_not_installed("ggplot2")
    df <- data.frame(gene = paste0("g", 1:10), mean_difference = rnorm(10), padj = runif(10))
    p <- TSENAT::plot_volcano(df)
    expect_s3_class(p, "ggplot")
})

skip_on_bioc()

context("Visualization: Generate Plots Extended Tests")

library(TSENAT)

# plot_top_transcripts writes output_file for single and multiple genes

test_that("plot_top_transcripts writes output files for single and multiple genes", {
    skip_if_not_installed("ggplot2")
    set.seed(10)
    counts <- matrix(rpois(6 * 4, lambda = 10), nrow = 6)
    rownames(counts) <- paste0("tx", 1:6)
    colnames(counts) <- paste0("S", 1:4)
    samples <- c("N", "N", "T", "T")
    tx2 <- data.frame(Transcript = rownames(counts), Gen = rep(c("G1", "G2"), each = 3), stringsAsFactors = FALSE)

    se <- SummarizedExperiment(
        assays = list(counts = counts),
        rowData = S4Vectors::DataFrame(genes = tx2$Gen),
        colData = S4Vectors::DataFrame(sample_type = samples)
    )

    tf1 <- tempfile(fileext = ".png")
    plot_top_transcripts(se, gene = "G1", output_file = tf1)
    expect_true(file.exists(tf1) && file.info(tf1)$size > 0)

    tf2 <- tempfile(fileext = ".png")
    plot_top_transcripts(se, gene = c("G1", "G2"), output_file = tf2)
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

    se <- SummarizedExperiment(
        assays = list(counts = counts),
        rowData = S4Vectors::DataFrame(genes = tx2$Gen),
        colData = S4Vectors::DataFrame(sample_type = samples)
    )
    p <- plot_top_transcripts(se, gene = "G1", metric = "iqr")
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

    p <- TSENAT:::.plot_ma_core(x, fc_df = fc)
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

context("Visualization: Top Transcripts Helper Functions")

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

context("Visualization: Plot Helper Functions")

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

test_that(".tsenat_prepare_volcano_df errors when x_col is not found in data", {
    # Test the error: stop(sprintf("Column '%s' not found in diff_df", x_col))
    df <- data.frame(
        gene = c("g1", "g2", "g3"),
        log2fc = c(0.5, -0.3, 0.8),
        adjusted_p_values = c(0.01, 0.5, 0.001),
        stringsAsFactors = FALSE
    )
    
    # Explicitly provide non-existent x_col
    expect_error(
        .tsenat_prepare_volcano_df(df, x_col = "missing_column"),
        "Column 'missing_column' not found in diff_df"
    )
})

test_that(".tsenat_prepare_volcano_df errors when padj_col is not found in data", {
    # Test the error: stop(sprintf("Column '%s' not found in diff_df", padj_col))
    df <- data.frame(
        gene = c("g1", "g2", "g3"),
        log2fc = c(0.5, -0.3, 0.8),
        pvalue = c(0.01, 0.5, 0.001),
        stringsAsFactors = FALSE
    )
    
    # Use default padj_col which doesn't exist
    expect_error(
        .tsenat_prepare_volcano_df(df, x_col = "log2fc"),
        "Column 'adjusted_p_values' not found in diff_df"
    )
    
    # Explicitly provide non-existent padj_col
    expect_error(
        .tsenat_prepare_volcano_df(df, x_col = "log2fc", padj_col = "wrong_padj"),
        "Column 'wrong_padj' not found in diff_df"
    )
})

test_that(".tsenat_prepare_volcano_df handles all valid column combinations", {
    # Test with various valid column names to ensure error catching is precise
    df <- data.frame(
        gene = c("g1", "g2", "g3"),
        mean_difference = c(0.5, -0.3, 0.8),
        p_adj = c(0.01, 0.5, 0.001),
        stringsAsFactors = FALSE
    )
    
    # Should work with valid columns
    result <- .tsenat_prepare_volcano_df(df, x_col = "mean_difference", padj_col = "p_adj")
    expect_is(result, "list")
    expect_true("df" %in% names(result))
    expect_equal(result$x_col, "mean_difference")
    expect_equal(result$padj_col, "p_adj")
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

context("Visualization: Unit Tests for Plotting Helpers")

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

context("Visualization: Gene Profile Plotting (Edge Cases)")

library(testthat)

context("Visualization: generate_plots.R Comprehensive Coverage")

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
    expect_error(plot_volcano(df), "Column 'padj' not found in diff_df")
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

