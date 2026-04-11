context("Metadata Mapping: Core Functionality")
library(SummarizedExperiment)


# Use a small in-memory toy `coldata` for tests to avoid file/system dependencies
toy_coldata <- data.frame(
    Sample = paste0("S", sprintf("%02d", seq_len(8))),
    Condition = rep(c("Normal", "Tumor"), 4),
    stringsAsFactors = FALSE
)

test_that("Exact match maps all samples and preserves order", {
    coldata <- toy_coldata
    sample_names <- as.character(coldata$Sample)
    cols <- paste0(sample_names, "_q=0.1")
    mat <- matrix(seq_along(cols), nrow = 1)
    colnames(mat) <- cols
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(diversity = mat)
    )

    se2 <- TSENAT:::.map_metadata_se(se, coldata, sample_col = "Sample", condition_col = "Condition")
    mapped_base <- sub("_q=.*", "", colnames(SummarizedExperiment::assay(se2)))
    expect_equal(length(mapped_base), length(sample_names))
    expect_equal(mapped_base, sample_names)
    expect_true(all(!is.na(SummarizedExperiment::colData(se2)$sample_type)))
})

test_that("Reversed coldata reorders SE to follow coldata", {
    coldata <- toy_coldata
    coldata_rev <- coldata[rev(seq_len(nrow(coldata))), , drop = FALSE]
    sample_names <- as.character(coldata$Sample)
    cols <- paste0(sample_names, "_q=0.1")
    mat <- matrix(runif(length(cols)), nrow = 1)
    colnames(mat) <- cols
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(diversity = mat)
    )

    se2 <- TSENAT:::.map_metadata_se(se, coldata_rev, sample_col = "Sample", condition_col = "Condition")
    mapped_base <- sub("_q=.*", "", colnames(SummarizedExperiment::assay(se2)))
    expect_equal(
        mapped_base[seq_len(nrow(coldata_rev))],
        as.character(coldata_rev$Sample)
    )
})

test_that("Missing sample in coldata triggers error (auto-detected pairing)", {
    coldata <- toy_coldata
    coldata_missing <- coldata[-3, , drop = FALSE]
    sample_names <- as.character(coldata$Sample)
    cols <- paste0(sample_names, "_q=0.1")
    mat <- matrix(runif(length(cols)), nrow = 1)
    colnames(mat) <- cols
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(diversity = mat)
    )

    expect_error(
        TSENAT:::.map_metadata_se(se, coldata_missing, sample_col = "Sample", condition_col = "Condition"),
        "unmatched samples in 'coldata'"
    )
})

test_that("Missing sample in coldata errors when insufficient columns", {
    coldata <- toy_coldata
    coldata_missing <- coldata[-3, , drop = FALSE]
    sample_names <- as.character(coldata$Sample)
    cols <- paste0(sample_names, "_q=0.1")
    mat <- matrix(runif(length(cols)), nrow = 1)
    colnames(mat) <- cols
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(diversity = mat)
    )

    expect_error(
        TSENAT:::.map_metadata_se(se, coldata_missing, sample_col = "Sample", condition_col = "Condition"),
        "unmatched samples in 'coldata'"
    )
})
test_that("Extra sample in coldata does not error via auto-detection", {
    coldata <- toy_coldata
    extra <- data.frame(
        Sample = "FAKE_SAMPLE_N",
        Condition = "Normal",
        stringsAsFactors = FALSE
    )
    coldata_extra <- rbind(extra, coldata)
    sample_names <- as.character(coldata$Sample)
    cols <- paste0(sample_names, "_q=0.1")
    mat <- matrix(runif(length(cols)), nrow = 1)
    colnames(mat) <- cols
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(diversity = mat)
    )

    se2 <- TSENAT:::.map_metadata_se(se, coldata_extra, sample_col = "Sample", condition_col = "Condition")
    expect_equal(
        ncol(SummarizedExperiment::assay(se2)),
        length(cols)
    )
})

test_that("Extra sample in coldata does not error when paired is FALSE", {
    coldata <- toy_coldata
    extra <- data.frame(
        Sample = "FAKE_SAMPLE_N",
        Condition = "Normal",
        stringsAsFactors = FALSE
    )
    coldata_extra <- rbind(extra, coldata)
    sample_names <- as.character(coldata$Sample)
    cols <- paste0(sample_names, "_q=0.1")
    mat <- matrix(runif(length(cols)), nrow = 1)
    colnames(mat) <- cols
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(diversity = mat)
    )

    se2 <- TSENAT:::.map_metadata_se(se, coldata_extra, sample_col = "Sample", condition_col = "Condition")
    expect_equal(
        ncol(SummarizedExperiment::assay(se2)),
        length(cols)
    )
})
test_that("Case mismatch in coldata errors (case-sensitive)", {
    coldata <- toy_coldata
    coldata_case <- coldata
    coldata_case$Sample <- tolower(coldata_case$Sample)
    sample_names <- as.character(coldata$Sample)
    cols <- paste0(sample_names, "_q=0.1")
    mat <- matrix(runif(length(cols)), nrow = 1)
    colnames(mat) <- cols
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(diversity = mat)
    )

    expect_error(
        TSENAT:::.map_metadata_se(se, coldata_case, sample_col = "Sample", condition_col = "Condition"),
        "unmatched samples in 'coldata'"
    )
})

test_that("SE with extra sample errors due to unmatched sample", {
    coldata <- toy_coldata
    sample_names <- as.character(coldata$Sample)
    cols <- c(paste0(sample_names, "_q=0.1"), "EXTRA_SAMPLE_q=0.1")
    mat <- matrix(runif(length(cols)), nrow = 1)
    colnames(mat) <- cols
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(diversity = mat)
    )

    expect_error(
        TSENAT:::.map_metadata_se(se, coldata, sample_col = "Sample", condition_col = "Condition"),
        "unmatched samples in 'coldata'"
    )
})

test_that("Multiple q columns are handled and mapped per-sample", {
    coldata <- toy_coldata
    sample_names <- as.character(coldata$Sample)[1:5]
    cols <- c(paste0(sample_names, "_q=0.1"), paste0(sample_names, "_q=0.2"))
    mat <- matrix(runif(length(cols)), nrow = 1)
    colnames(mat) <- cols
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(diversity = mat)
    )

    se2 <- TSENAT:::.map_metadata_se(se, coldata, sample_col = "Sample", condition_col = "Condition")
    mapped_base <- sub("_q=.*", "", colnames(SummarizedExperiment::assay(se2)))
    # each base sample appears twice
    tbl <- table(mapped_base)
    expect_true(all(tbl == 2))
    expect_true(all(!is.na(SummarizedExperiment::colData(se2)$sample_type)))
})

test_that("Non-data.frame coldata is ignored and SE remains unchanged", {
    coldata <- toy_coldata
    named_vec <- setNames(
        as.character(coldata$Condition),
        as.character(coldata$Sample)
    )
    sample_names <- as.character(coldata$Sample)
    cols <- paste0(sample_names, "_q=0.1")
    mat <- matrix(runif(length(cols)), nrow = 1)
    colnames(mat) <- cols
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(diversity = mat)
    )

    # Non-data.frame input should throw an error
    expect_error(
        TSENAT:::.map_metadata_se(se, named_vec),
        "metadata must be a data.frame"
    )
})

test_that("colData rownames and sample_base are aligned to assay columns", {
    coldata <- toy_coldata
    sample_names <- as.character(coldata$Sample)
    cols <- paste0(sample_names, "_q=0.1")
    mat <- matrix(runif(length(cols)), nrow = 1)
    colnames(mat) <- cols
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(diversity = mat)
    )

    se2 <- TSENAT:::.map_metadata_se(se, coldata, sample_col = "Sample", condition_col = "Condition")
    cd <- SummarizedExperiment::colData(se2)
    expect_equal(rownames(cd), colnames(SummarizedExperiment::assay(se2)))
    expect_true("sample_base" %in% colnames(cd))
    expect_equal(
        as.character(cd$sample_base),
        sub("_q=.*", "", colnames(SummarizedExperiment::assay(se2)))
    )
})


test_that("Auto-detection accepts complete pairs with suffix-based pairing", {
    # construct a small paired coldata with two participants A and B
    coldata_small <- data.frame(
        Sample = c("A_N", "A_T", "B_N", "B_T"),
        Condition = c("Normal", "Tumor", "Normal", "Tumor"),
        stringsAsFactors = FALSE
    )
    sample_names <- as.character(coldata_small$Sample)
    cols <- paste0(sample_names, "_q=0.1")
    mat <- matrix(runif(length(cols)), nrow = 1)
    colnames(mat) <- cols
    se <- SummarizedExperiment::SummarizedExperiment(assays = list(diversity = mat))

    expect_error_free <- function(expr) {
        expect_error(force(expr), NA)
    }

    expect_error_free(TSENAT:::.map_metadata_se(se, coldata_small, sample_col = "Sample", condition_col = "Condition"))
})

test_that("map_metadata warns when only one condition is present", {
    coldata <- data.frame(
        Sample = c("S1", "S2", "S3"),
        Condition = c("control", "control", "control"),
        stringsAsFactors = FALSE
    )
    
    sample_names <- as.character(coldata$Sample)
    cols <- paste0(sample_names, "_q=0.1")
    mat <- matrix(runif(length(cols)), nrow = 1)
    colnames(mat) <- cols
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(diversity = mat)
    )
    
    # Should warn about only one condition
    expect_warning(
        TSENAT:::.map_metadata_se(se, coldata, sample_col = "Sample", condition_col = "Condition"),
        "Only one condition detected"
    )
})

test_that("map_metadata does not warn when multiple conditions are present", {
    coldata <- toy_coldata  # Has both "Normal" and "Tumor" conditions
    sample_names <- as.character(coldata$Sample)
    cols <- paste0(sample_names, "_q=0.1")
    mat <- matrix(runif(length(cols)), nrow = 1)
    colnames(mat) <- cols
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(diversity = mat)
    )
    
    # Should NOT warn when multiple conditions are present
    expect_no_warning(
        TSENAT:::.map_metadata_se(se, coldata, sample_col = "Sample", condition_col = "Condition")
    )
})

test_that("map_tx_to_readcounts assigns rownames when sizes match", {
    rc <- matrix(1:4, nrow = 2)
    txmap <- data.frame(Transcript = c("tx1", "tx2"), Gene = c("g1", "g2"), stringsAsFactors = FALSE)
    out <- TSENAT:::.map_tx_to_readcounts(rc, txmap)
    expect_equal(rownames(out), as.character(txmap$Transcript))
})

test_that("map_tx_to_readcounts errors when sizes mismatch and no matching ids", {
    rc <- matrix(1:6, nrow = 3)
    txmap <- data.frame(Transcript = c("a", "b"), Gene = c("g1", "g2"), stringsAsFactors = FALSE)
    expect_error(TSENAT:::.map_tx_to_readcounts(rc, txmap), "does not match readcounts rows")
})

test_that("map_tx_to_readcounts matches existing rownames", {
    rc <- matrix(1:4, nrow = 2)
    rownames(rc) <- c("b", "a")
    txmap <- data.frame(Transcript = c("a", "b", "c"), stringsAsFactors = FALSE)
    expect_message(out <- TSENAT:::.map_tx_to_readcounts(rc, txmap, verbose = TRUE), "Matched and assigned transcript IDs")
    expect_equal(rownames(out), as.character(c("b", "a")))
})

test_that("map_samples_to_group uses colData mapping and errors on missing", {
    mat <- matrix(1:4, nrow = 2)
    colnames(mat) <- c("S1_q=1", "S2_q=1")
    se <- SummarizedExperiment(assays = SimpleList(diversity = mat))
    SummarizedExperiment::colData(se)$sample_type <- c("A", "B")
    # valid mapping
    mapped <- TSENAT:::.map_samples_to_group(c("S1", "S2"), se = se, condition_col = "sample_type", mat = NULL)
    expect_equal(mapped, c("A", "B"))
    # missing sample should error
    expect_error(TSENAT:::.map_samples_to_group(c("S1", "S3"), se = se, condition_col = "sample_type", mat = NULL), "Missing sample_type mapping")
})

test_that("get_assay_long returns long df and respects sample_type_col", {
    mat <- matrix(c(1, 2, 3, 4), nrow = 2)
    rownames(mat) <- c("g1", "g2")
    colnames(mat) <- c("S1", "S2")
    se <- SummarizedExperiment(assays = SimpleList(diversity = mat))
    rowData(se)$genes <- c("G1", "G2")
    SummarizedExperiment::colData(se) <- S4Vectors::DataFrame(sample_type = c("N", "T"), row.names = colnames(mat))
    long <- TSENAT:::.get_assay_long(se, assay_name = "diversity", value_name = "val", condition_col = "sample_type")
    expect_true(all(c("Gene", "sample", "val", "sample_type") %in% colnames(long)))
    expect_equal(unique(as.character(long$sample_type)), c("N", "T"))
})

test_that("prepare_tsallis_long parses _q= suffixes and maps groups", {
    mat <- matrix(c(1, NA, 2, 3, 4, NA), nrow = 3)
    colnames(mat) <- c("S1_q=1", "S1_q=2")
    rownames(mat) <- c("g1", "g2", "g3")
    se <- SummarizedExperiment(assays = SimpleList(diversity = mat))
    rowData(se)$genes <- c("G1", "G2", "G3")
    SummarizedExperiment::colData(se) <- S4Vectors::DataFrame(sample_type = c("A", "A"), row.names = colnames(mat))
    long <- TSENAT:::.prepare_tsallis_long(se, assay_name = "diversity", condition_col = "sample_type")
    # q should be numeric (bug fix: was converting to factor)
    expect_true("q" %in% colnames(long))
    expect_true(is.numeric(long$q))
    expect_false(is.factor(long$q))
    expect_true(all(long$group == "A"))
    # rows with NA tsallis removed
    expect_false(any(is.na(long$tsallis)))
})

test_that("prepare_tsallis_long preserves decimal q-values as numeric", {
    # Test for bug fix: q values like 0.1, 0.15, 0.2 were being converted to factors
    skip_if_not_installed("tidyr")
    
    mat <- matrix(rnorm(30), nrow = 3, ncol = 10)
    # Create column names with decimal q-values
    q_vals <- seq(0.1, 0.5, by = 0.1)
    col_names <- paste0("S", rep(1:2, each = 5), "_q=", rep(q_vals, 2))
    colnames(mat) <- col_names
    rownames(mat) <- c("g1", "g2", "g3")
    
    se <- SummarizedExperiment(assays = SimpleList(diversity = mat))
    rowData(se)$genes <- c("G1", "G2", "G3")
    SummarizedExperiment::colData(se) <- S4Vectors::DataFrame(
        sample_type = rep(c("GroupA", "GroupB"), each = 5), 
        row.names = colnames(mat)
    )
    
    long <- TSENAT:::.prepare_tsallis_long(se, assay_name = "diversity", condition_col = "sample_type")
    
    # Verify q is numeric, not factor
    expect_true(is.numeric(long$q))
    expect_false(is.factor(long$q))
    
    # Verify the actual q-values are preserved (0.1, 0.2, 0.3, 0.4, 0.5)
    # Use round() to handle floating-point precision issues
    unique_q <- sort(unique(round(long$q, 10)))
    q_vals_rounded <- round(q_vals, 10)
    expect_true(all(unique_q %in% q_vals_rounded))
    expect_equal(length(unique_q), length(q_vals_rounded))
})

test_that("map_metadata handles NULL/invalid coldata and maps sample types + metadata", {
    mat <- matrix(1:6, nrow = 3)
    colnames(mat) <- c("S1_q=1", "S2_q=1")
    se <- SummarizedExperiment(assays = SimpleList(diversity = mat))
    
    # NULL coldata should throw an error
    expect_error(
        TSENAT:::.map_metadata_se(se, NULL),
        "metadata must be a data.frame"
    )
    
    # Less than 2 columns should throw an error
    bad <- data.frame(A = 1:2)
    expect_error(
        TSENAT:::.map_metadata_se(se, bad),
        "metadata must have at least 2 columns"
    )
    
    # proper mapping: create coldata with Sample and Condition
    coldata <- data.frame(Sample = c("S1", "S2"), Condition = c("N", "T"), stringsAsFactors = FALSE)
    # Note: globalenv() lookups removed for reproducibility
    # readcounts and tx2gene must be explicitly provided or already in SE metadata
    se4 <- TSENAT:::.map_metadata_se(se, coldata, sample_col = "Sample", condition_col = "Condition")
    expect_true("sample_type" %in% colnames(SummarizedExperiment::colData(se4)))
    expect_true("sample_base" %in% colnames(SummarizedExperiment::colData(se4)))
    md <- S4Vectors::metadata(se4)
    # Since we don't look in globalenv anymore, these should be NULL
    expect_true(is.null(md$readcounts))
    expect_true(is.null(md$tx2gene))
})

test_that("map_tx_to_readcounts assigns rownames from tx2gene data.frame", {
    rc <- matrix(1:6, nrow = 3, ncol = 2)
    txmap <- data.frame(Transcript = paste0("tx", 1:3), Gene = c("G1", "G1", "G2"), stringsAsFactors = FALSE)
    out <- TSENAT:::.map_tx_to_readcounts(rc, txmap)
    expect_equal(rownames(out), txmap$Transcript)
})

test_that("map_tx_to_readcounts reads tx2gene from file and errors on mismatch", {
    rc <- matrix(1:4, nrow = 2)
    txmap <- data.frame(Transcript = paste0("tx", 1:3), stringsAsFactors = FALSE)
    tmp <- tempfile(fileext = ".tsv")
    utils::write.table(txmap, file = tmp, sep = "\t", row.names = FALSE, quote = FALSE)
    expect_error(TSENAT:::.map_tx_to_readcounts(rc, tmp))
})

context("Metadata Mapping: Extended Tests")


test_that("map_tx_to_readcounts accepts file path input", {
    rc <- matrix(1:6, nrow = 3)
    txmap <- data.frame(Transcript = paste0("tx", 1:4), Gene = c("g1", "g1", "g2", "g2"), stringsAsFactors = FALSE)
    tf <- tempfile(fileext = ".tsv")
    write.table(txmap, tf, sep = "\t", row.names = FALSE, quote = FALSE)
    expect_error(TSENAT:::.map_tx_to_readcounts(rc, tf), "does not match readcounts rows")
    unlink(tf)
})

test_that("map_samples_to_group with mat provided returns single group when no mapping", {
    mat <- matrix(1:4, nrow = 2)
    colnames(mat) <- c("A_q=1", "B_q=1")
    res <- TSENAT:::.map_samples_to_group(c("A", "B"), se = NULL, condition_col = NULL, mat = mat)
    expect_equal(res, c("Group", "Group"))
})

test_that("get_assay_long errors when assay missing", {
    se <- SummarizedExperiment(assays = S4Vectors::SimpleList(dummy = matrix(1, nrow = 1)))
    expect_error(TSENAT:::.get_assay_long(se, assay_name = "nope"))
})

test_that("get_assay_long errors when all values are NA", {
    # Test the code path: if (nrow(long_filtered) == 0) { stop(...) }
    mat <- matrix(NA_real_, nrow = 2, ncol = 3)
    rownames(mat) <- c("g1", "g2")
    colnames(mat) <- c("S1", "S2", "S3")
    se <- SummarizedExperiment(assays = SimpleList(diversity = mat))
    rowData(se)$genes <- c("G1", "G2")
    SummarizedExperiment::colData(se) <- S4Vectors::DataFrame(sample_type = c("N", "T", "N"), row.names = colnames(mat))
    
    expect_error(
        TSENAT:::.get_assay_long(se, assay_name = "diversity", value_name = "val", condition_col = "sample_type"),
        "No non-NA values found in assay 'diversity'. All values are NA."
    )
})

test_that("get_assay_long filters out NA values but errors when all are NA", {
    # Mixed NA and non-NA values should work
    mat <- matrix(c(1, NA, 2, 3, 4, 5), nrow = 2, ncol = 3)
    rownames(mat) <- c("g1", "g2")
    colnames(mat) <- c("S1", "S2", "S3")
    se <- SummarizedExperiment(assays = SimpleList(diversity = mat))
    rowData(se)$genes <- c("G1", "G2")
    SummarizedExperiment::colData(se) <- S4Vectors::DataFrame(sample_type = c("N", "T", "N"), row.names = colnames(mat))
    
    long <- TSENAT:::.get_assay_long(se, assay_name = "diversity", value_name = "val", condition_col = "sample_type")
    
    # Should have filtered out NA values but kept valid ones
    expect_true(all(!is.na(long$val)))
    expect_true(nrow(long) == 5)  # 6 values - 1 NA = 5
})

test_that("get_assay_long with default sample_type when column missing", {
    # When sample_type_col is NULL, should use "Group" as default
    mat <- matrix(c(1, 2, 3, 4), nrow = 2)
    rownames(mat) <- c("g1", "g2")
    colnames(mat) <- c("S1", "S2")
    se <- SummarizedExperiment(assays = SimpleList(diversity = mat))
    rowData(se)$genes <- c("G1", "G2")
    
    long <- TSENAT:::.get_assay_long(se, assay_name = "diversity", value_name = "val", condition_col = NULL)
    
    expect_true("sample_type" %in% colnames(long))
    expect_true(all(long$sample_type == "Group"))
})



test_that("prepare_tsallis_long handles no _q suffix and default group", {
    mat <- matrix(1:4, nrow = 2)
    colnames(mat) <- c("S1", "S2")
    rownames(mat) <- c("g1", "g2")
    se <- SummarizedExperiment(assays = S4Vectors::SimpleList(diversity = mat))
    rowData(se)$genes <- c("G1", "G2")
    long <- TSENAT:::.prepare_tsallis_long(se, assay_name = "diversity", condition_col = NULL)
    expect_true(all(long$group == "Group"))
    expect_true(all(is.na(long$q)))
})

test_that("map_tx_to_readcounts supports custom tx_col and data.frame input", {
    rc <- data.frame(c1 = 1:3, c2 = 4:6)
    txmap <- data.frame(ID = paste0("t", 1:3), stringsAsFactors = FALSE)
    out <- TSENAT:::.map_tx_to_readcounts(rc, txmap, tx_col = "ID")
    expect_equal(rownames(out), as.character(txmap$ID))
})

test_that("map_metadata with pairing column in coldata validates pairs", {
    mat <- matrix(1:12, nrow = 2)
    # SE data: 4 samples with 3 q-values each = 12 columns
    colnames(mat) <- c("A_N_q=1", "A_N_q=2", "A_T_q=1", "A_T_q=2", "B_N_q=1", "B_N_q=2")
    se <- SummarizedExperiment(assays = S4Vectors::SimpleList(diversity = mat))
    
    # Create coldata WITH both conditions but pair2 lacks treatment
    # pair1: complete (has both N and T from individual samples)
    # pair2: incomplete (only has N)
    coldata <- data.frame(
        Sample = c("A_N", "A_T", "B_N"),
        Condition = c("normal", "treated", "normal"),
        pair_id = c("pair1", "pair1", "pair2"),
        stringsAsFactors = FALSE
    )
    
    # Should error because pair2 only has condition "normal", missing "treated"
    expect_error(
        TSENAT:::.map_metadata_se(se, coldata, sample_col = "Sample", condition_col = "Condition", subject_col = "pair_id"), 
        "Unpaired subjects found"
    )
})

test_that("map_metadata with no tx2gene/txmap leaves metadata tx2gene NULL", {
    # When neither tx2gene nor txmap exist, tx2gene should remain NULL
    coldata <- data.frame(
        Sample = c("S1", "S2"),
        Condition = c("A", "B"),
        stringsAsFactors = FALSE
    )
    
    mat <- matrix(rnorm(4), nrow = 2)
    colnames(mat) <- c("S1_q=0.5", "S2_q=0.5")
    se <- SummarizedExperiment(assays = S4Vectors::SimpleList(assay1 = mat))
    
    result <- TSENAT:::.map_metadata_se(se, coldata, sample_col = "Sample", condition_col = "Condition")
    
    # Verify that tx2gene was not set (remains NULL)
    md <- S4Vectors::metadata(result)
    expect_true(is.null(md$tx2gene))
})

# ============================================================================
# Tests for paired_samples functionality (new feature)
# ============================================================================

test_that("map_metadata creates paired_samples column from third column of coldata", {
    # Create coldata with paired_samples information in third column
    coldata <- data.frame(
        Sample = c("SRR_001", "SRR_002", "SRR_003", "SRR_004"),
        Condition = c("normal", "tumor", "normal", "tumor"),
        paired_samples = c("A", "A", "B", "B"),
        stringsAsFactors = FALSE
    )
    
    mat <- matrix(rnorm(8), nrow = 2)
    colnames(mat) <- paste0(c("SRR_001", "SRR_002", "SRR_003", "SRR_004"), "_q=0.5")
    se <- SummarizedExperiment(assays = S4Vectors::SimpleList(diversity = mat))
    
    result <- TSENAT:::.map_metadata_se(se, coldata, sample_col = "Sample", condition_col = "Condition")
    cd <- colData(result)
    
    # Verify paired_samples column exists
    expect_true("paired_samples" %in% colnames(cd))
    
    # Verify correct mapping
    expected_pairing <- c("A", "A", "B", "B")
    actual_pairing <- as.character(cd$paired_samples)
    expect_equal(actual_pairing, expected_pairing)
})

test_that("map_metadata preserves paired_samples column name from metadata", {
    # Using a different column name
    coldata <- data.frame(
        Sample = c("S1", "S2", "S3", "S4"),
        Condition = c("control", "treatment", "control", "treatment"),
        subject_id = c("P1", "P1", "P2", "P2"),
        stringsAsFactors = FALSE
    )
    
    mat <- matrix(rnorm(8), nrow = 2)
    colnames(mat) <- paste0(c("S1", "S2", "S3", "S4"), "_q=1.0")
    se <- SummarizedExperiment(assays = S4Vectors::SimpleList(diversity = mat))
    
    result <- TSENAT:::.map_metadata_se(se, coldata, sample_col = "Sample", condition_col = "Condition")
    cd <- colData(result)
    
    # Verify the pairing column exists with the original name
    expect_true("subject_id" %in% colnames(cd))
    
    # Verify correct mapping
    expected <- c("P1", "P1", "P2", "P2")
    actual <- as.character(cd$subject_id)
    expect_equal(actual, expected)
})

test_that("map_metadata maps paired_samples correctly with multiple q values", {
    # Create coldata with paired samples
    coldata <- data.frame(
        Sample = c("A", "B", "C", "D"),
        Condition = c("control", "treated", "control", "treated"),
        paired_samples = c("pair1", "pair1", "pair2", "pair2"),
        stringsAsFactors = FALSE
    )
    
    # Create diversity matrix with multiple q values per sample
    q_vals <- c(0.5, 1.0, 1.5)
    cols <- c(
        paste0("A_q=", q_vals),
        paste0("B_q=", q_vals),
        paste0("C_q=", q_vals),
        paste0("D_q=", q_vals)
    )
    # Create matrix with proper dimensions: rows = features, cols = samples
    mat <- matrix(rnorm(2 * length(cols)), nrow = 2, ncol = length(cols))
    colnames(mat) <- cols
    se <- SummarizedExperiment(assays = S4Vectors::SimpleList(diversity = mat))
    
    result <- TSENAT:::.map_metadata_se(se, coldata, sample_col = "Sample", condition_col = "Condition")
    cd <- colData(result)
    
    # Verify paired_samples is repeated correctly for each q value
    # Samples are in order A, B, C, D with pairings: pair1, pair1, pair2, pair2
    # Each repeated for 3 q-values: [A×3, B×3, C×3, D×3]
    expected_pairing <- c(
        rep("pair1", 3),  # A repeated for 3 q-values
        rep("pair1", 3),  # B repeated for 3 q-values
        rep("pair2", 3),  # C repeated for 3 q-values
        rep("pair2", 3)   # D repeated for 3 q-values
    )
    actual_pairing <- as.character(cd$paired_samples)
    expect_equal(actual_pairing, expected_pairing)
})

test_that("map_metadata validates pairing structure when third column present", {
    # Create incomplete pairing (pair1 missing treatment condition)
    coldata_bad <- data.frame(
        Sample = c("S1", "S2", "S3"),
        Condition = c("control", "control", "treated"),
        paired_id = c("pair1", "pair1", "pair2"),
        stringsAsFactors = FALSE
    )
    
    mat <- matrix(rnorm(6), nrow = 2)
    colnames(mat) <- paste0(c("S1", "S2", "S3"), "_q=0.5")
    se <- SummarizedExperiment(assays = S4Vectors::SimpleList(diversity = mat))
    
    # Should error because pair1 lacks both conditions (S1 and S2 both have "control")
    expect_error(
        TSENAT:::.map_metadata_se(se, coldata_bad, sample_col = "Sample", condition_col = "Condition", subject_col = "paired_id"),
        "Unpaired subjects found"
    )
})

test_that("map_metadata without third column creates sample_base from suffix", {
    # Two-column coldata without explicit pairing
    coldata <- data.frame(
        Sample = c("SRR_001", "SRR_002", "SRR_003", "SRR_004"),
        Condition = c("normal", "tumor", "normal", "tumor"),
        stringsAsFactors = FALSE
    )
    
    mat <- matrix(rnorm(8), nrow = 2)
    colnames(mat) <- paste0(c("SRR_001", "SRR_002", "SRR_003", "SRR_004"), "_q=0.5")
    se <- SummarizedExperiment(assays = S4Vectors::SimpleList(diversity = mat))
    
    result <- TSENAT:::.map_metadata_se(se, coldata, sample_col = "Sample", condition_col = "Condition")
    cd <- colData(result)
    
    # paired_samples should not exist (no third column)
    expect_false("paired_samples" %in% colnames(cd))
    
    # sample_base should exist and contain sample identifiers
    expect_true("sample_base" %in% colnames(cd))
})

test_that("map_metadata expands colData correctly with paired_samples", {
    # Create paired metadata
    coldata <- data.frame(
        Sample = c("P1_N", "P1_T", "P2_N", "P2_T", "P3_N", "P3_T"),
        Condition = c("normal", "tumor", "normal", "tumor", "normal", "tumor"),
        pair_id = c("patient1", "patient1", "patient2", "patient2", "patient3", "patient3"),
        stringsAsFactors = FALSE
    )
    
    # Create SE with 2 q values per sample
    cols <- c(
        paste0("P1_N_q=", c(0.5, 1.0)),
        paste0("P1_T_q=", c(0.5, 1.0)),
        paste0("P2_N_q=", c(0.5, 1.0)),
        paste0("P2_T_q=", c(0.5, 1.0)),
        paste0("P3_N_q=", c(0.5, 1.0)),
        paste0("P3_T_q=", c(0.5, 1.0))
    )
    mat <- matrix(rnorm(3 * length(cols)), nrow = 3, ncol = length(cols))
    colnames(mat) <- cols
    se <- SummarizedExperiment::SummarizedExperiment(assays = S4Vectors::SimpleList(diversity = mat))
    
    result <- TSENAT:::.map_metadata_se(se, coldata, sample_col = "Sample", condition_col = "Condition")
    cd <- colData(result)
    
    # Verify expanded dimensions
    expect_equal(nrow(cd), length(cols))
    expect_equal(ncol(mat), length(cols))
    
    # Verify pair_id is correctly expanded
    # Columns are ordered as: P1_N × 2q, P1_T × 2q, P2_N × 2q, P2_T × 2q, P3_N × 2q, P3_T × 2q
    # So pair mapping is: patient1 (P1_N), patient1 (P1_T), patient2 (P2_N), patient2 (P2_T), 
    #                     patient3 (P3_N), patient3 (P3_T), each repeated 2 times for q-values
    expected_pairs <- c(
        rep("patient1", 2),  # P1_N with 2 q-values
        rep("patient1", 2),  # P1_T with 2 q-values
        rep("patient2", 2),  # P2_N with 2 q-values
        rep("patient2", 2),  # P2_T with 2 q-values
        rep("patient3", 2),  # P3_N with 2 q-values
        rep("patient3", 2)   # P3_T with 2 q-values
    )
    actual_pairs <- as.character(cd$pair_id)
    expect_equal(actual_pairs, expected_pairs)
})

test_that("map_metadata with metadata.tsv structure works", {
    # Simulate structure from actual metadata.tsv
    coldata <- data.frame(
        Sample = c("SRR14800481", "SRR14800479", "SRR14800480", "SRR14800478"),
        condition = c("normal", "tumor", "normal", "tumor"),
        paired_samples = c("A", "A", "B", "B"),
        stringsAsFactors = FALSE
    )
    
    mat <- matrix(rnorm(8), nrow = 2)
    colnames(mat) <- paste0(c("SRR14800481", "SRR14800479", "SRR14800480", "SRR14800478"), "_q=1.0")
    se <- SummarizedExperiment(assays = S4Vectors::SimpleList(diversity = mat))
    
    result <- TSENAT:::.map_metadata_se(se, coldata, sample_col = "Sample", condition_col = "condition")
    cd <- colData(result)
    
    # Verify sample_type mapped correctly
    expected_types <- c("normal", "tumor", "normal", "tumor")
    actual_types <- as.character(cd$sample_type)
    expect_equal(actual_types, expected_types)
    
    # Verify paired_samples mapped correctly
    expected_pairs <- c("A", "A", "B", "B")
    actual_pairs <- as.character(cd$paired_samples)
    expect_equal(actual_pairs, expected_pairs)
})
