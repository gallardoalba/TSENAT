context("map_metadata behavior")

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

    se2 <- map_metadata(se, coldata)
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

    se2 <- map_metadata(se, coldata_rev)
    mapped_base <- sub("_q=.*", "", colnames(SummarizedExperiment::assay(se2)))
    expect_equal(
        mapped_base[seq_len(nrow(coldata_rev))],
        as.character(coldata_rev$Sample)
    )
})

test_that("Missing sample in coldata triggers error (unpaired)", {
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
        map_metadata(se, coldata_missing, paired = TRUE),
        "Unpaired samples found in coldata"
    )
})

test_that("Missing sample in coldata errors when paired is FALSE", {
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
        map_metadata(se, coldata_missing),
        "unmatched samples in 'coldata'"
    )
})
test_that("Extra sample in coldata triggers error (unpaired)", {
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

    expect_error(
        map_metadata(se, coldata_extra, paired = TRUE),
        "Unpaired samples found in coldata"
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

    se2 <- map_metadata(se, coldata_extra)
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
        map_metadata(se, coldata_case),
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
        map_metadata(se, coldata),
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

    se2 <- map_metadata(se, coldata)
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

    se2 <- map_metadata(se, named_vec)
    # mapping should be a no-op because map_coldata_to_se expects a data.frame
    expect_equal(
        colnames(SummarizedExperiment::assay(se2)),
        colnames(SummarizedExperiment::assay(se))
    )
    expect_null(SummarizedExperiment::colData(se2)$sample_type)
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

    se2 <- map_metadata(se, coldata)
    cd <- SummarizedExperiment::colData(se2)
    expect_equal(rownames(cd), colnames(SummarizedExperiment::assay(se2)))
    expect_true("sample_base" %in% colnames(cd))
    expect_equal(
        as.character(cd$sample_base),
        sub("_q=.*", "", colnames(SummarizedExperiment::assay(se2)))
    )
})


test_that("paired = TRUE accepts complete pairs", {
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

    expect_error_free(map_metadata(se, coldata_small, paired = TRUE))
})

context("map_functions additional tests")

library(SummarizedExperiment)

test_that("map_tx_to_readcounts assigns rownames when sizes match", {
    rc <- matrix(1:4, nrow = 2)
    txmap <- data.frame(Transcript = c("tx1", "tx2"), Gene = c("g1", "g2"), stringsAsFactors = FALSE)
    out <- map_tx_to_readcounts(rc, txmap)
    expect_equal(rownames(out), as.character(txmap$Transcript))
})

test_that("map_tx_to_readcounts errors when sizes mismatch and no matching ids", {
    rc <- matrix(1:6, nrow = 3)
    txmap <- data.frame(Transcript = c("a", "b"), Gene = c("g1", "g2"), stringsAsFactors = FALSE)
    expect_error(map_tx_to_readcounts(rc, txmap), "does not match readcounts rows")
})

test_that("map_tx_to_readcounts matches existing rownames", {
    rc <- matrix(1:4, nrow = 2)
    rownames(rc) <- c("b", "a")
    txmap <- data.frame(Transcript = c("a", "b", "c"), stringsAsFactors = FALSE)
    expect_message(out <- map_tx_to_readcounts(rc, txmap, verbose = TRUE), "Matched and assigned transcript IDs")
    expect_equal(rownames(out), as.character(c("b", "a")))
})

test_that("map_samples_to_group uses colData mapping and errors on missing", {
    mat <- matrix(1:4, nrow = 2)
    colnames(mat) <- c("S1_q=1", "S2_q=1")
    se <- SummarizedExperiment(assays = SimpleList(diversity = mat))
    SummarizedExperiment::colData(se)$sample_type <- c("A", "B")
    # valid mapping
    mapped <- map_samples_to_group(c("S1", "S2"), se = se, sample_type_col = "sample_type", mat = NULL)
    expect_equal(mapped, c("A", "B"))
    # missing sample should error
    expect_error(map_samples_to_group(c("S1", "S3"), se = se, sample_type_col = "sample_type", mat = NULL), "Missing sample_type mapping")
})

test_that("get_assay_long returns long df and respects sample_type_col", {
    mat <- matrix(c(1, 2, 3, 4), nrow = 2)
    rownames(mat) <- c("g1", "g2")
    colnames(mat) <- c("S1", "S2")
    se <- SummarizedExperiment(assays = SimpleList(diversity = mat))
    rowData(se)$genes <- c("G1", "G2")
    SummarizedExperiment::colData(se) <- S4Vectors::DataFrame(sample_type = c("N", "T"), row.names = colnames(mat))
    long <- get_assay_long(se, assay_name = "diversity", value_name = "val", sample_type_col = "sample_type")
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
    long <- prepare_tsallis_long(se, assay_name = "diversity", sample_type_col = "sample_type")
    # q should be a factor and group mapping present
    expect_true("q" %in% colnames(long))
    expect_true(is.factor(long$q))
    expect_true(all(long$group == "A"))
    # rows with NA tsallis removed
    expect_false(any(is.na(long$tsallis)))
})

test_that("map_metadata handles NULL/invalid coldata and maps sample types + metadata", {
    mat <- matrix(1:6, nrow = 3)
    colnames(mat) <- c("S1_q=1", "S2_q=1")
    se <- SummarizedExperiment(assays = SimpleList(diversity = mat))
    # NULL coldata returns same
    se2 <- map_metadata(se, NULL)
    expect_identical(se2, se)
    # invalid coldata returns same
    bad <- data.frame(A = 1:2, B = c("x", "y"))
    se3 <- map_metadata(se, bad)
    expect_identical(se3, se)
    # proper mapping: create coldata with Sample and Condition
    coldata <- data.frame(Sample = c("S1", "S2"), Condition = c("N", "T"), stringsAsFactors = FALSE)
    # attach readcounts and tx2gene into globalenv to test metadata attaching
    readcounts <- matrix(1:6, nrow = 3)
    tx2gene <- data.frame(Transcript = paste0("tx", 1:3), Gen = c("g1", "g1", "g2"), stringsAsFactors = FALSE)
    assign("readcounts", readcounts, envir = globalenv())
    assign("tx2gene", tx2gene, envir = globalenv())
    on.exit(
        {
            rm(readcounts, envir = globalenv())
            rm(tx2gene, envir = globalenv())
        },
        add = TRUE
    )
    se4 <- map_metadata(se, coldata)
    expect_true("sample_type" %in% colnames(SummarizedExperiment::colData(se4)))
    expect_true("sample_base" %in% colnames(SummarizedExperiment::colData(se4)))
    md <- S4Vectors::metadata(se4)
    expect_true(!is.null(md$readcounts))
    expect_true(!is.null(md$tx2gene))
})

test_that("map_tx_to_readcounts assigns rownames from tx2gene data.frame", {
    rc <- matrix(1:6, nrow = 3, ncol = 2)
    txmap <- data.frame(Transcript = paste0("tx", 1:3), Gene = c("G1", "G1", "G2"), stringsAsFactors = FALSE)
    out <- map_tx_to_readcounts(rc, txmap)
    expect_equal(rownames(out), txmap$Transcript)
})

test_that("map_tx_to_readcounts reads tx2gene from file and errors on mismatch", {
    rc <- matrix(1:4, nrow = 2)
    txmap <- data.frame(Transcript = paste0("tx", 1:3), stringsAsFactors = FALSE)
    tmp <- tempfile(fileext = ".tsv")
    utils::write.table(txmap, file = tmp, sep = "\t", row.names = FALSE, quote = FALSE)
    expect_error(map_tx_to_readcounts(rc, tmp))
})

context("map_functions more tests")

library(SummarizedExperiment)

test_that("map_tx_to_readcounts accepts file path input", {
    rc <- matrix(1:6, nrow = 3)
    txmap <- data.frame(Transcript = paste0("tx", 1:4), Gene = c("g1", "g1", "g2", "g2"), stringsAsFactors = FALSE)
    tf <- tempfile(fileext = ".tsv")
    write.table(txmap, tf, sep = "\t", row.names = FALSE, quote = FALSE)
    expect_error(map_tx_to_readcounts(rc, tf), "does not match readcounts rows")
    unlink(tf)
})

test_that("map_samples_to_group with mat provided returns single group when no mapping", {
    mat <- matrix(1:4, nrow = 2)
    colnames(mat) <- c("A_q=1", "B_q=1")
    res <- map_samples_to_group(c("A", "B"), se = NULL, sample_type_col = NULL, mat = mat)
    expect_equal(res, c("Group", "Group"))
})

test_that("get_assay_long errors when assay missing", {
    se <- SummarizedExperiment(assays = S4Vectors::SimpleList(dummy = matrix(1, nrow = 1)))
    expect_error(get_assay_long(se, assay_name = "nope"))
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
        get_assay_long(se, assay_name = "diversity", value_name = "val", sample_type_col = "sample_type"),
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
    
    long <- get_assay_long(se, assay_name = "diversity", value_name = "val", sample_type_col = "sample_type")
    
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
    
    long <- get_assay_long(se, assay_name = "diversity", value_name = "val", sample_type_col = NULL)
    
    expect_true("sample_type" %in% colnames(long))
    expect_true(all(long$sample_type == "Group"))
})



test_that("prepare_tsallis_long handles no _q suffix and default group", {
    mat <- matrix(1:4, nrow = 2)
    colnames(mat) <- c("S1", "S2")
    rownames(mat) <- c("g1", "g2")
    se <- SummarizedExperiment(assays = S4Vectors::SimpleList(diversity = mat))
    rowData(se)$genes <- c("G1", "G2")
    long <- prepare_tsallis_long(se, assay_name = "diversity", sample_type_col = NULL)
    expect_true(all(long$group == "Group"))
    expect_true(all(is.na(long$q)))
})

test_that("map_tx_to_readcounts supports custom tx_col and data.frame input", {
    rc <- data.frame(c1 = 1:3, c2 = 4:6)
    txmap <- data.frame(ID = paste0("t", 1:3), stringsAsFactors = FALSE)
    out <- map_tx_to_readcounts(rc, txmap, tx_col = "ID")
    expect_equal(rownames(out), as.character(txmap$ID))
})

test_that("map_metadata paired validation errors on unpaired bases", {
    mat <- matrix(1:8, nrow = 2)
    colnames(mat) <- c("B1_q=1", "B2_q=1", "B1_q=2", "B2_q=2")
    se <- SummarizedExperiment(assays = S4Vectors::SimpleList(diversity = mat))
    # create coldata missing one condition for B2
    coldata <- data.frame(Sample = c("B1", "B1", "B2"), Condition = c("N", "T", "N"), stringsAsFactors = FALSE)
    expect_error(map_metadata(se, coldata, paired = TRUE), "Unpaired samples found in coldata")
})
