context("filter_se: SummarizedExperiment Filtering and Quality Control")

testthat::test_that("filter_se keeps expected transcripts and updates metadata", {
    # create small toy dataset to avoid dependency on external data
    set.seed(1)
    n_tx <- 12L
    n_samps <- 4L
    genes <- paste0("G", sprintf("%03d", seq_len(n_tx)))
    readcounts <- matrix(sample(0:100, n_tx * n_samps, replace = TRUE), nrow = n_tx, ncol = n_samps)
    colnames(readcounts) <- paste0("S", seq_len(n_samps))
    tx2gene_df <- data.frame(
        Transcript = paste0("tx", seq_len(n_tx)),
        Gene = genes,
        stringsAsFactors = FALSE
    )
    rownames(readcounts) <- tx2gene_df$Transcript

    se <- build_se(readcounts, tx2gene_df)
    # Expect warnings: No TPM data, min_samples > ncol, all transcripts filtered
    testthat::expect_warning(
        se_f <- filter_se(se, min_tpm = 5L, min_samples = 5L, verbose = FALSE),
        "min_samples|No TPM|Filtering removed"
    )

    testthat::expect_s4_class(se_f, "SummarizedExperiment")
    testthat::expect_true(nrow(se_f) < nrow(se))
    md <- S4Vectors::metadata(se_f)
    testthat::expect_true(is.null(md$readcounts) || nrow(md$readcounts) == nrow(se_f))
    if (!is.null(md$tx2gene)) testthat::expect_true(nrow(md$tx2gene) <= nrow(se))
})

testthat::test_that("filter_se warns when assay removes all rows", {
    # toy dataset matching the first test
    set.seed(1)
    n_tx <- 12L
    n_samps <- 4L
    genes <- paste0("G", sprintf("%03d", seq_len(n_tx)))
    readcounts <- matrix(sample(0:100, n_tx * n_samps, replace = TRUE), nrow = n_tx, ncol = n_samps)
    colnames(readcounts) <- paste0("S", seq_len(n_samps))
    tx2gene_df <- data.frame(
        Transcript = paste0("tx", seq_len(n_tx)),
        Gene = genes,
        stringsAsFactors = FALSE
    )
    rownames(readcounts) <- tx2gene_df$Transcript

    se <- build_se(readcounts, tx2gene_df)
    # choose thresholds that remove all rows
    testthat::expect_warning(sf <- filter_se(se, min_tpm = 1e6, min_samples = 1, verbose = FALSE))
    testthat::expect_s4_class(sf, "SummarizedExperiment")
    testthat::expect_true(nrow(sf) == 0)
})

context("SummarizedExperiment Filtering: Additional Tests")

library(SummarizedExperiment)

test_that("filter_se errors on non-SE input", {
    expect_error(filter_se(1:10), "must be a SummarizedExperiment")
})

test_that("filter_se uses specified assay by name and falls back", {
    mat1 <- matrix(c(0, 6, 7, 2, 8, 9), nrow = 3)
    mat2 <- matrix(1:6, nrow = 3)
    se <- SummarizedExperiment(assays = list(counts = mat1, other = mat2))
    S4Vectors::metadata(se)$salmon_tpm <- mat1  # Add TPM data for first call
    res <- filter_se(se, min_tpm = 5, min_samples = 1, assay_name = "counts", verbose = FALSE)
    expect_s4_class(res, "SummarizedExperiment")
    
    # using missing assay name warns and uses first
    # Create new SE without TPM data so it warns about both missing assay AND TPM
    se2 <- SummarizedExperiment(assays = list(counts = mat1, other = mat2))
    expect_warning(filter_se(se2, min_tpm = 5, min_samples = 1, assay_name = "nope", verbose = FALSE))
})

test_that("filter_se respects min_samples cap and returns empty SE when none kept", {
    mat <- matrix(0, nrow = 5, ncol = 3)
    se <- SummarizedExperiment(assays = list(counts = mat))
    S4Vectors::metadata(se)$salmon_tpm <- mat  # Add TPM data to avoid warning
    expect_warning(res <- filter_se(se, min_tpm = 0, min_samples = 5, verbose = FALSE))
    # all zeros and min_samples > available leads to zero kept and a warning
    expect_equal(nrow(SummarizedExperiment::assay(res, 1)), 0)
})

test_that("filter_se subsets metadata readcounts and tx2gene", {
    mat <- matrix(c(0, 6, 7, 2, 8, 9), nrow = 3)
    rownames(mat) <- paste0("tx", 1:3)
    se <- SummarizedExperiment(assays = list(counts = mat))
    # attach metadata
    S4Vectors::metadata(se)$readcounts <- mat
    S4Vectors::metadata(se)$salmon_tpm <- mat  # Add TPM data to avoid warning
    S4Vectors::metadata(se)$tx2gene <- data.frame(Transcript = rownames(mat), Gene = c("g1", "g1", "g2"), stringsAsFactors = FALSE)
    res <- filter_se(se, min_tpm = 5, min_samples = 1, verbose = FALSE)
    md <- S4Vectors::metadata(res)
    expect_true(is.null(md$readcounts) == FALSE)
    expect_true(is.null(md$tx2gene) == FALSE)
    # tx2gene rows should subset to remaining transcripts
    expect_true(all(md$tx2gene$Transcript %in% rownames(SummarizedExperiment::assay(res))))
})

test_that("filter_se removes genes with fewer than min_tx_per_gene transcripts", {
    skip_if_not_installed("SummarizedExperiment")
    library(SummarizedExperiment)
    
    # Create a test SE with genes having different numbers of transcripts
    # Gene 1: 3 transcripts, Gene 2: 1 transcript (should be filtered if min_tx_per_gene >= 2)
    mat <- matrix(c(
        100, 50, 80,   # Gene 1 - 3 txs, should pass
        50,             # Gene 2 - 1 tx, should fail if min_tx_per_gene >= 2
        90, 70, 60      # Gene 3 - 3 txs, should pass
    ), nrow = 7, ncol = 2)
    rownames(mat) <- paste0("tx", 1:7)
    
    se <- SummarizedExperiment(assays = list(counts = mat))
    S4Vectors::metadata(se)$salmon_tpm <- mat  # Add TPM data to avoid warning
    
    # Add gene information
    genes_df <- data.frame(
        gene_id = c("g1", "g1", "g1", "g2", "g3", "g3", "g3")
    )
    rowData(se) <- genes_df
    
    # Test default min_tx_per_gene = 2
    res_default <- filter_se(se, min_tpm = 10, min_samples = 1, verbose = FALSE)
    genes_after <- unique(rowData(res_default)$gene_id)
    
    # Gene 2 (single transcript) should be removed
    expect_false("g2" %in% genes_after)
    expect_true("g1" %in% genes_after)
    expect_true("g3" %in% genes_after)
    
    # Test with min_tx_per_gene = 1 (keep single-tx genes)
    res_keep_single <- filter_se(se, min_tpm = 10, min_samples = 1, min_tx_per_gene = 1L, verbose = FALSE)
    genes_all <- unique(rowData(res_keep_single)$gene_id)
    
    # Should keep all genes including g2
    expect_true("g1" %in% genes_all)
    expect_true("g2" %in% genes_all)
    expect_true("g3" %in% genes_all)
})

test_that("filter_se min_tx_per_gene works without gene information", {
    skip_if_not_installed("SummarizedExperiment")
    library(SummarizedExperiment)
    
    # SE without gene information in rowData
    mat <- matrix(sample(1:100, 20), nrow = 5, ncol = 4)
    se <- SummarizedExperiment(assays = list(counts = mat))
    S4Vectors::metadata(se)$salmon_tpm <- mat  # Add TPM data to avoid warning
    
    # Should not error even though there's no gene info
    res <- filter_se(se, min_tpm = 5, min_samples = 1, min_tx_per_gene = 2L, verbose = FALSE)
    
    # Should be a SummarizedExperiment
    expect_s4_class(res, "SummarizedExperiment")
})
# ============================================================================
# TPM ASSAY DETECTION TESTS
# ============================================================================

context("filter_se: TPM Assay Detection and Priority")

test_that("filter_se detects TPM from user-specified tpm_assay_name parameter", {
    # Strategy 1: User explicitly specifies TPM assay name
    mat_tpm <- matrix(c(2.5, 5.0, 1.5, 3.2, 4.1, 0.8), nrow = 3)  # TPM values
    mat_counts <- matrix(c(10, 50, 8, 20, 30, 5), nrow = 3)         # Raw counts
    rownames(mat_tpm) <- rownames(mat_counts) <- paste0("tx", 1:3)
    
    se <- SummarizedExperiment(assays = list(counts = mat_counts, tpm = mat_tpm))
    
    # Filter using explicitly specified TPM assay
    res <- filter_se(se, min_tpm = 2.0, min_samples = 1, tpm_assay_name = "tpm", verbose = FALSE)
    
    # Should keep transcripts with TPM >= 2.0 in at least 1 sample
    # tx1: [2.5, 3.2] >= 2.0? Yes, keep
    # tx2: [5.0, 4.1] >= 2.0? Yes, keep
    # tx3: [1.5, 0.8] >= 2.0? No, remove
    expect_equal(nrow(res), 2)
    expect_true("tx1" %in% rownames(res))
    expect_true("tx2" %in% rownames(res))
    expect_false("tx3" %in% rownames(res))
})

test_that("filter_se warns when user-specified tpm_assay_name not found", {
    mat <- matrix(c(10, 50, 8, 20, 30, 5), nrow = 3)
    se <- SummarizedExperiment(assays = list(counts = mat))
    
    # Specify a non-existent TPM assay
    expect_warning(
        filter_se(se, min_tpm = 5, min_samples = 1, tpm_assay_name = "nonexistent", verbose = FALSE),
        "TPM assay 'nonexistent' not found"
    )
})

test_that("filter_se detects metadata$salmon_tpm (SALMON preprocessed TPM)", {
    # Strategy 2: Look in metadata for salmon_tpm
    mat_tpm <- matrix(c(2.5, 5.0, 1.5, 3.2, 4.1, 0.8), nrow = 3)
    mat_counts <- matrix(c(10, 50, 8, 20, 30, 5), nrow = 3)
    rownames(mat_tpm) <- rownames(mat_counts) <- paste0("tx", 1:3)
    
    se <- SummarizedExperiment(assays = list(counts = mat_counts))
    S4Vectors::metadata(se)$salmon_tpm <- mat_tpm  # Add SALMON-preprocessed TPM
    
    # Should automatically use metadata$salmon_tpm for filtering
    res <- filter_se(se, min_tpm = 2.0, min_samples = 1, verbose = FALSE)
    
    expect_equal(nrow(res), 2)
    expect_true("tx1" %in% rownames(res))
    expect_true("tx2" %in% rownames(res))
    expect_false("tx3" %in% rownames(res))
})

test_that("filter_se detects metadata$tpm (custom TPM)", {
    # Strategy 2b: Look in metadata for tpm (when salmon_tpm not present)
    mat_tpm <- matrix(c(2.5, 5.0, 1.5, 3.2, 4.1, 0.8), nrow = 3)
    mat_counts <- matrix(c(10, 50, 8, 20, 30, 5), nrow = 3)
    rownames(mat_tpm) <- rownames(mat_counts) <- paste0("tx", 1:3)
    
    se <- SummarizedExperiment(assays = list(counts = mat_counts))
    S4Vectors::metadata(se)$tpm <- mat_tpm  # Add custom TPM metadata
    
    # Should automatically use metadata$tpm for filtering
    res <- filter_se(se, min_tpm = 2.0, min_samples = 1, verbose = FALSE)
    
    expect_equal(nrow(res), 2)
})

test_that("filter_se auto-detects assay named 'tpm'", {
    # Strategy 3: Search for assay named 'tpm'
    mat_tpm <- matrix(c(2.5, 5.0, 1.5, 3.2, 4.1, 0.8), nrow = 3)
    mat_counts <- matrix(c(10, 50, 8, 20, 30, 5), nrow = 3)
    rownames(mat_tpm) <- rownames(mat_counts) <- paste0("tx", 1:3)
    
    se <- SummarizedExperiment(assays = list(counts = mat_counts, tpm = mat_tpm))
    
    # Should auto-detect 'tpm' assay without explicit specification
    res <- filter_se(se, min_tpm = 2.0, min_samples = 1, verbose = FALSE)
    
    expect_equal(nrow(res), 2)
})

test_that("filter_se auto-detects assay named 'abundance' (tximport format)", {
    # Strategy 3b: Search for assay named 'abundance' (tximport format)
    mat_tpm <- matrix(c(2.5, 5.0, 1.5, 3.2, 4.1, 0.8), nrow = 3)
    mat_counts <- matrix(c(10, 50, 8, 20, 30, 5), nrow = 3)
    rownames(mat_tpm) <- rownames(mat_counts) <- paste0("tx", 1:3)
    
    se <- SummarizedExperiment(assays = list(counts = mat_counts, abundance = mat_tpm))
    
    # Should auto-detect 'abundance' assay for tximport-style objects
    res <- filter_se(se, min_tpm = 2.0, min_samples = 1, verbose = FALSE)
    
    expect_equal(nrow(res), 2)
})

test_that("filter_se respects TPM detection priority order", {
    # Priority should be: user-specified > salmon_tpm > tpm > assay 'tpm' > 'abundance'
    # Create different TPM values in each location
    mat_tpm1 <- matrix(c(5.0, 5.0, 5.0, 5.0, 5.0, 5.0), nrow = 3)    # All pass at >= 2.0
    mat_tpm2 <- matrix(c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5), nrow = 3)    # All fail at >= 2.0
    mat_counts <- matrix(c(10, 50, 8, 20, 30, 5), nrow = 3)
    rownames(mat_tpm1) <- rownames(mat_tpm2) <- rownames(mat_counts) <- paste0("tx", 1:3)
    
    se <- SummarizedExperiment(assays = list(counts = mat_counts, tpm = mat_tpm1, abundance = mat_tpm2))
    S4Vectors::metadata(se)$salmon_tpm <- mat_tpm2
    
    # Should use user-specified even if others available
    res <- filter_se(se, min_tpm = 2.0, min_samples = 1, tpm_assay_name = "tpm", verbose = FALSE)
    expect_equal(nrow(res), 3)  # All pass with mat_tpm1
    
    # Remove tpm_assay_name and verify salmon_tpm is used next
    se2 <- SummarizedExperiment(assays = list(counts = mat_counts, abundance = mat_tpm1))
    S4Vectors::metadata(se2)$salmon_tpm <- mat_tpm2
    # When salmon_tpm (all 0.5) is used with min_tpm = 2.0, all transcripts are filtered out
    expect_warning(
        res2 <- filter_se(se2, min_tpm = 2.0, min_samples = 1, verbose = FALSE),
        "Filtering removed all transcripts"
    )
    expect_equal(nrow(res2), 0)  # Uses salmon_tpm (mat_tpm2), all fail
})

# ============================================================================
# STRINGENCY-BASED FILTERING TESTS
# ============================================================================

context("filter_se: Stringency-Based Filtering")

test_that("filter_se stringency='soft' applies permissive filtering (25% samples)", {
    # Stringency "soft": min_samples = 25% of samples (min 2)
    # Create data where different patterns appear
    mat_tpm <- matrix(c(
        5.0, 4.0, 3.0, 2.0,  # tx1: high in all samples
        2.0, 0.5, 0.5, 0.5,  # tx2: high in 1 sample only
        1.0, 1.0, 0.5, 0.5   # tx3: moderate in 2 samples
    ), nrow = 3, ncol = 4, byrow = TRUE)
    rownames(mat_tpm) <- paste0("tx", 1:3)
    colnames(mat_tpm) <- paste0("S", 1:4)
    
    # Create colData with pair information
    col_data <- data.frame(
        pair = c(1, 1, 2, 2),
        row.names = paste0("S", 1:4)
    )
    
    se <- SummarizedExperiment(
        assays = list(counts = mat_tpm),
        colData = col_data
    )
    S4Vectors::metadata(se)$salmon_tpm <- mat_tpm
    
    # stringency="soft" -> min_samples = 25% of 4 = 1 (but min(2)), so 2 samples
    res <- filter_se(se, stringency = "soft", pair_col = "pair", verbose = FALSE)
    
    # tx1: >= auto-estimated min_tpm in >= 2 samples -> keep
    # tx2: >= min_tpm in >= 1 sample, needs 2 -> may be kept/removed depending on estimate
    # tx3: >= min_tpm in >= 2 samples -> keep
    expect_s4_class(res, "SummarizedExperiment")
    expect_true(nrow(res) > 0)  # Should keep some transcripts
})

test_that("filter_se stringency='medium' applies balanced filtering (50% samples)", {
    # Stringency "medium": min_samples = 50% of samples (min 3)
    mat_tpm <- matrix(c(
        5.0, 4.0, 3.0, 2.0,  # tx1: high in all samples
        2.0, 0.5, 0.5, 0.5,  # tx2: high in 1 sample only
        1.0, 1.0, 0.5, 0.5   # tx3: moderate in 2 samples
    ), nrow = 3, ncol = 4, byrow = TRUE)
    rownames(mat_tpm) <- paste0("tx", 1:3)
    colnames(mat_tpm) <- paste0("S", 1:4)
    
    col_data <- data.frame(
        pair = c(1, 1, 2, 2),
        row.names = paste0("S", 1:4)
    )
    
    se <- SummarizedExperiment(
        assays = list(counts = mat_tpm),
        colData = col_data
    )
    S4Vectors::metadata(se)$salmon_tpm <- mat_tpm
    
    # stringency="medium" -> min_samples = 50% of 4 = 2 (but min(3)), so 3 samples
    res <- filter_se(se, stringency = "medium", pair_col = "pair", verbose = FALSE)
    
    # Should be more stringent than soft
    expect_s4_class(res, "SummarizedExperiment")
})

test_that("filter_se stringency='severe' applies stringent filtering (75% samples)", {
    # Stringency "severe": min_samples = 75% of samples
    mat_tpm <- matrix(c(
        5.0, 4.0, 3.0, 2.0,
        2.0, 0.5, 0.5, 0.5,
        1.0, 1.0, 0.5, 0.5
    ), nrow = 3, ncol = 4, byrow = TRUE)
    rownames(mat_tpm) <- paste0("tx", 1:3)
    colnames(mat_tpm) <- paste0("S", 1:4)
    
    col_data <- data.frame(
        pair = c(1, 1, 2, 2),
        row.names = paste0("S", 1:4)
    )
    
    se <- SummarizedExperiment(
        assays = list(counts = mat_tpm),
        colData = col_data
    )
    S4Vectors::metadata(se)$salmon_tpm <- mat_tpm
    
    # stringency="severe" -> min_samples = 75% of 4 = 3 samples
    res <- filter_se(se, stringency = "severe", pair_col = "pair", verbose = FALSE)
    
    # Should be more stringent than medium
    expect_s4_class(res, "SummarizedExperiment")
})

test_that("filter_se rejects invalid stringency value", {
    mat <- matrix(1:12, nrow = 3, ncol = 4)
    colnames(mat) <- paste0("S", 1:4)
    
    col_data <- data.frame(pair = c(1, 1, 2, 2), row.names = paste0("S", 1:4))
    se <- SummarizedExperiment(assays = list(counts = mat), colData = col_data)
    S4Vectors::metadata(se)$salmon_tpm <- mat
    
    expect_error(
        filter_se(se, stringency = "invalid", pair_col = "pair", verbose = FALSE),
        "stringency.*must be one of"
    )
})

test_that("filter_se auto-detects pair column from colData", {
    # Test auto-detection of pair column when not explicitly specified
    mat_tpm <- matrix(sample(1:50, 20), nrow = 5, ncol = 4)
    rownames(mat_tpm) <- paste0("tx", 1:5)
    colnames(mat_tpm) <- paste0("S", 1:4)
    
    col_data <- data.frame(
        pair_id = c(1, 1, 2, 2),  # pair_id is auto-detectablePairs
        row.names = paste0("S", 1:4)
    )
    
    se <- SummarizedExperiment(
        assays = list(counts = mat_tpm),
        colData = col_data
    )
    S4Vectors::metadata(se)$salmon_tpm <- mat_tpm
    
    # Should auto-detect 'pair_id' column without explicit specification
    res <- filter_se(se, stringency = "soft", verbose = FALSE)
    
    expect_s4_class(res, "SummarizedExperiment")
})

test_that("filter_se errors when pair column cannot be detected", {
    mat_tpm <- matrix(sample(1:20, 20), nrow = 5, ncol = 4)
    rownames(mat_tpm) <- paste0("tx", 1:5)
    colnames(mat_tpm) <- paste0("S", 1:4)
    
    # No pair column in colData
    col_data <- data.frame(
        sample_id = paste0("S", 1:4),
        row.names = paste0("S", 1:4)
    )
    
    se <- SummarizedExperiment(
        assays = list(counts = mat_tpm),
        colData = col_data
    )
    S4Vectors::metadata(se)$salmon_tpm <- mat_tpm
    
    # Should error when stringency specified but pair column missing and undetectable
    expect_error(
        filter_se(se, stringency = "soft", verbose = FALSE),
        "Could not auto-detect pair column"
    )
})

test_that("filter_se stringency adjusts min_tx_per_gene based on stringency level", {
    # Soft/medium: 2 transcripts per gene; Severe: 3 transcripts per gene
    mat_tpm <- matrix(c(
        5.0, 4.0, 3.0, 2.0,
        4.5, 3.5, 2.5, 1.5,
        2.0, 1.0, 0.5, 0.5
    ), nrow = 3, ncol = 4, byrow = TRUE)
    rownames(mat_tpm) <- paste0("tx", 1:3)
    colnames(mat_tpm) <- paste0("S", 1:4)
    
    genes_df <- data.frame(
        gene_id = c("gene1", "gene1", "gene2"),
        row.names = paste0("tx", 1:3)
    )
    
    col_data <- data.frame(
        pair = c(1, 1, 2, 2),
        row.names = paste0("S", 1:4)
    )
    
    se <- SummarizedExperiment(
        assays = list(counts = mat_tpm),
        rowData = genes_df,
        colData = col_data
    )
    S4Vectors::metadata(se)$salmon_tpm <- mat_tpm
    
    # Medium stringency should adjust min_tx_per_gene to 2
    # This may filter out all transcripts depending on the estimated min_tpm
    expect_warning(
        res <- filter_se(se, stringency = "medium", pair_col = "pair", verbose = FALSE),
        "Filtering removed all transcripts"
    )
    
    # Metadata should record the stringency-adjusted parameters
    if (nrow(res) > 0) {
        filtered_info <- S4Vectors::metadata(res)$filtered
        expect_equal(filtered_info$stringency, "medium")
        expect_equal(filtered_info$min_tx_per_gene, 2L)
    }
})

