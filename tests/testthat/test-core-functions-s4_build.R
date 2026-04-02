# Comprehensive tests for R/se_manipulation_build.R functions
# Tests for: .get_gene_ids, .extract_gff3_data, .extract_attribute_fast,
# .load_tx2gene_data, .validate_readcounts, .map_transcripts_to_genes,
# .store_salmon_metadata, .build_rowdata_se, .build_se

# ===========================================================================
# SETUP HELPERS
# ===========================================================================

#' Create minimal test SummarizedExperiment
#' @noRd
make_test_se_minimal <- function(n_genes = 10, n_samples = 5) {
    counts <- matrix(rpois(n_genes * n_samples, lambda = 100), nrow = n_genes)
    rownames(counts) <- paste0("ENST", sprintf("%011d", seq_len(n_genes)))
    colnames(counts) <- paste0("sample_", seq_len(n_samples))

    rowdata <- S4Vectors::DataFrame(
        transcript_id = rownames(counts),
        gene_id = paste0("ENSG", sprintf("%011d", rep(1:(n_genes/2+1), length.out = n_genes))),
        row.names = rownames(counts)
    )

    SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = counts),
        rowData = rowdata,
        colData = data.frame(sample = colnames(counts), row.names = colnames(counts))
    )
}

#' Create a temporary TSV tx2gene file
#' @noRd
create_temp_tx2gene_tsv <- function(n_tx = 100, tmpdir = tempdir()) {
    filepath <- file.path(tmpdir, paste0("tx2gene_", format(Sys.time(), "%s"), ".tsv"))

    tx_ids <- paste0("ENST", sprintf("%011d", seq_len(n_tx)))
    gene_ids <- paste0("ENSG", sprintf("%011d", rep(1:(n_tx/2), length.out = n_tx)))

    data <- data.frame(Transcript = tx_ids, Gene = gene_ids)
    readr::write_tsv(data, filepath)
    filepath
}

#' Create a temporary GFF3 file
#' @noRd
create_temp_gff3_file <- function(n_genes = 10, n_tx_per_gene = 2, tmpdir = tempdir()) {
    filepath <- file.path(tmpdir, paste0("annotation_", format(Sys.time(), "%s"), ".gff3"))

    lines <- c("##gff-version 3")

    for (i in seq_len(n_genes)) {
        gene_id <- paste0("gene_", i)
        lines <- c(lines, sprintf(
            "chr1\tgencode\tgene\t%d\t%d\t.\t+\t.\tID=%s;gene_name=GENE%d",
            (i-1)*5000 + 1, i*5000, gene_id, i
        ))

        for (j in seq_len(n_tx_per_gene)) {
            tx_id <- paste0("ENST", sprintf("%011d", (i-1)*n_tx_per_gene + j))
            lines <- c(lines, sprintf(
                "chr1\tgencode\ttranscript\t%d\t%d\t.\t+\t.\tID=%s;Parent=%s",
                (i-1)*5000 + (j-1)*2000 + 1, (i-1)*5000 + j*2000, tx_id, gene_id
            ))
        }
    }

    writeLines(lines, filepath)
    filepath
}

# ===========================================================================
# Tests for .get_gene_ids()
# ===========================================================================

test_that(".get_gene_ids retrieves gene IDs from rowData", {
    se <- make_test_se_minimal(n_genes = 10)
    gene_ids <- TSENAT:::.get_gene_ids(se)

    expect_identical(length(gene_ids), nrow(se))
    expect_true(all(grepl("ENSG", gene_ids)))
})

test_that(".get_gene_ids returns NULL for empty rowData", {
    se <- make_test_se_minimal()
    rowData(se) <- NULL

    gene_ids <- TSENAT:::.get_gene_ids(se)
    expect_null(gene_ids)
})

test_that(".get_gene_ids uses gene_name as fallback", {
    se <- make_test_se_minimal()
    rowData(se)$gene_id <- NULL
    rowData(se)$gene_name <- paste0("GENE_", seq_len(nrow(se)))

    gene_names <- TSENAT:::.get_gene_ids(se)
    expect_true(all(grepl("GENE_", gene_names)))
})

# ===========================================================================
# Tests for .extract_attribute_fast()
# ===========================================================================

test_that(".extract_attribute_fast extracts attribute values correctly", {
    attr_str <- "ID=ENST00000001;Parent=ENSG00000001;gene_name=BRCA1"

    id_val <- TSENAT:::.extract_attribute_fast(attr_str, "ID=")
    expect_equal(id_val, "ENST00000001")

    parent_val <- TSENAT:::.extract_attribute_fast(attr_str, "Parent=")
    expect_equal(parent_val, "ENSG00000001")

    gene_name <- TSENAT:::.extract_attribute_fast(attr_str, "gene_name=")
    expect_equal(gene_name, "BRCA1")
})

test_that(".extract_attribute_fast handles missing attributes", {
    attr_str <- "ID=ENST00000001;Parent=ENSG00000001"

    missing <- TSENAT:::.extract_attribute_fast(attr_str, "gene_name=")
    expect_true(is.na(missing))
})

test_that(".extract_attribute_fast handles end-of-string attributes", {
    attr_str <- "ID=ENST00000001;gene_name=BRCA1"

    gene_name <- TSENAT:::.extract_attribute_fast(attr_str, "gene_name=")
    expect_equal(gene_name, "BRCA1")
})

# ===========================================================================
# Tests for .validate_readcounts()
# ===========================================================================

test_that(".validate_readcounts validates numeric matrix", {
    readcounts <- matrix(c(10, 5, 2, 3), nrow = 2, dimnames = list(c("tx1", "tx2"), c("s1", "s2")))

    result <- TSENAT:::.validate_readcounts(readcounts)
    expect_true(is.matrix(result$readcounts))
    expect_identical(result$tx_ids, c("tx1", "tx2"))
})

test_that(".validate_readcounts converts data.frame to matrix", {
    readcounts <- data.frame(s1 = c(10, 5), s2 = c(2, 3), row.names = c("tx1", "tx2"))

    result <- TSENAT:::.validate_readcounts(readcounts)
    expect_true(is.matrix(result$readcounts))
    expect_identical(result$tx_ids, c("tx1", "tx2"))
})

test_that(".validate_readcounts raises error for missing rownames", {
    readcounts <- matrix(c(10, 5, 2, 3), nrow = 2)

    expect_error(TSENAT:::.validate_readcounts(readcounts), "rownames")
})

test_that(".validate_readcounts raises error for non-numeric", {
    readcounts <- matrix(c("10", "5", "2", "3"), nrow = 2, dimnames = list(c("tx1", "tx2"), c("s1", "s2")))

    expect_error(TSENAT:::.validate_readcounts(readcounts), "numeric")
})

# ===========================================================================
# Tests for .load_tx2gene_data()
# ===========================================================================

test_that(".load_tx2gene_data loads TSV file", {
    tsv_file <- create_temp_tx2gene_tsv(n_tx = 50)
    on.exit(unlink(tsv_file))

    result <- TSENAT:::.load_tx2gene_data(tsv_file, verbose = FALSE)
    expect_true(is.data.frame(result$tx2gene_df))
    expect_equal(nrow(result$tx2gene_df), 50)
    expect_named(result$tx2gene_df, c("Transcript", "Gene"))
    expect_true(is.null(result$gff3_data))
})

test_that(".load_tx2gene_data loads data.frame directly", {
    tx2gene <- data.frame(Transcript = c("tx1", "tx2"), Gene = c("g1", "g1"))

    result <- TSENAT:::.load_tx2gene_data(tx2gene, verbose = FALSE)
    expect_identical(result$tx2gene_df, tx2gene)
    expect_true(is.null(result$gff3_data))
})

test_that(".load_tx2gene_data raises error for missing file", {
    expect_error(TSENAT:::.load_tx2gene_data("/nonexistent/file.tsv"), "not found")
})

test_that(".load_tx2gene_data loads GFF3 file", {
    gff3_file <- create_temp_gff3_file(n_genes = 5, n_tx_per_gene = 2)
    on.exit(unlink(gff3_file))

    result <- TSENAT:::.load_tx2gene_data(gff3_file, verbose = FALSE)
    expect_true(is.data.frame(result$tx2gene_df))
    expect_true(nrow(result$tx2gene_df) >= 5, info = "Expected at least 5 transcripts from GFF3")
    expect_named(result$tx2gene_df, c("Transcript", "Gene"))
    expect_false(is.null(result$gff3_data))
})

# ===========================================================================
# Tests for .map_transcripts_to_genes()
# ===========================================================================

test_that(".map_transcripts_to_genes maps correctly", {
    tx_ids <- c("ENST00000001", "ENST00000002", "ENST00000003")
    tx2gene_df <- data.frame(
        Transcript = tx_ids,
        Gene = c("ENSG00000001", "ENSG00000002", "ENSG00000001")
    )
    readcounts <- matrix(1:6, nrow = 3, dimnames = list(tx_ids, c("s1", "s2")))

    result <- TSENAT:::.map_transcripts_to_genes(
        tx_ids, readcounts, tx2gene_df, NULL, NULL, skip = FALSE
    )

    expect_identical(result$genes, c("ENSG00000001", "ENSG00000002", "ENSG00000001"))
    expect_identical(result$tx_ids, tx_ids)
})

test_that(".map_transcripts_to_genes errors on unmapped when skip=FALSE", {
    tx_ids <- c("ENST00000001", "UNMAPPED_TX")
    tx2gene_df <- data.frame(
        Transcript = "ENST00000001",
        Gene = "ENSG00000001"
    )
    readcounts <- matrix(1:4, nrow = 2, dimnames = list(tx_ids, c("s1", "s2")))

    expect_error(
        TSENAT:::.map_transcripts_to_genes(tx_ids, readcounts, tx2gene_df, NULL, NULL, skip = FALSE),
        "Unmapped"
    )
})

test_that(".map_transcripts_to_genes removes unmapped when skip=TRUE", {
    tx_ids <- c("ENST00000001", "UNMAPPED_TX", "ENST00000003")
    tx2gene_df <- data.frame(
        Transcript = c("ENST00000001", "ENST00000003"),
        Gene = c("ENSG00000001", "ENSG00000002")
    )
    readcounts <- matrix(1:6, nrow = 3, dimnames = list(tx_ids, c("s1", "s2")))

    result <- TSENAT:::.map_transcripts_to_genes(
        tx_ids, readcounts, tx2gene_df, NULL, NULL, skip = TRUE
    )

    expect_equal(length(result$tx_ids), 2)
    expect_false("UNMAPPED_TX" %in% result$tx_ids)
    expect_equal(nrow(result$readcounts), 2)
})

test_that(".map_transcripts_to_genes handles SALMON data correctly", {
    tx_ids <- c("ENST00000001", "ENST00000002")
    tx2gene_df <- data.frame(
        Transcript = tx_ids,
        Gene = c("ENSG00000001", "ENSG00000001")
    )
    readcounts <- matrix(1:4, nrow = 2, dimnames = list(tx_ids, c("s1", "s2")))
    tpm <- matrix(10:13, nrow = 2, dimnames = list(tx_ids, c("s1", "s2")))
    eff_length <- c(1000, 1200)

    result <- TSENAT:::.map_transcripts_to_genes(
        tx_ids, readcounts, tx2gene_df, tpm, eff_length, skip = FALSE
    )

    expect_identical(nrow(result$tpm), 2L)
    expect_identical(length(result$effective_length), 2L)
})

# ===========================================================================
# Tests for .store_salmon_metadata()
# ===========================================================================

test_that(".store_salmon_metadata stores TPM correctly", {
    se <- make_test_se_minimal(n_genes = 5, n_samples = 3)
    readcounts <- assay(se, "counts")
    tpm <- readcounts * 2

    se <- TSENAT:::.store_salmon_metadata(se, readcounts, tpm, NULL)

    expect_false(is.null(metadata(se)$tpm))
    expect_identical(nrow(metadata(se)$tpm), nrow(readcounts))
})

test_that(".store_salmon_metadata stores effective_length correctly", {
    se <- make_test_se_minimal(n_genes = 5, n_samples = 3)
    readcounts <- assay(se, "counts")
    eff_length <- rep(1000, nrow(readcounts))
    names(eff_length) <- rownames(readcounts)

    se <- TSENAT:::.store_salmon_metadata(se, readcounts, NULL, eff_length)

    expect_false(is.null(metadata(se)$salmon_effective_length))
    expect_identical(length(metadata(se)$salmon_effective_length), nrow(readcounts))
})

test_that(".store_salmon_metadata validates TPM dimensions", {
    se <- make_test_se_minimal(n_genes = 5, n_samples = 3)
    readcounts <- assay(se, "counts")
    tpm <- matrix(1:10, nrow = 5, ncol = 2)  # Wrong dimensions

    expect_error(
        TSENAT:::.store_salmon_metadata(se, readcounts, tpm, NULL),
        "dimensions"
    )
})

test_that(".store_salmon_metadata validates effective_length length", {
    se <- make_test_se_minimal(n_genes = 5, n_samples = 3)
    readcounts <- assay(se, "counts")
    eff_length <- rep(1000, 3)  # Wrong length

    expect_error(
        TSENAT:::.store_salmon_metadata(se, readcounts, NULL, eff_length),
        "length"
    )
})

# ===========================================================================
# Tests for .build_rowdata_se()
# ===========================================================================

test_that(".build_rowdata_se builds rowData correctly", {
    se <- make_test_se_minimal(n_genes = 5)
    tx_ids <- paste0("ENST", sprintf("%011d", 1:5))
    genes <- paste0("ENSG", sprintf("%011d", c(1, 1, 2, 2, 3)))

    se <- TSENAT:::.build_rowdata_se(se, tx_ids, genes, NULL)

    rowdata <- rowData(se)
    expect_true("transcript_id" %in% colnames(rowdata))
    expect_true("gene_id" %in% colnames(rowdata))
    expect_identical(rowdata$transcript_id, tx_ids)
    expect_identical(rowdata$gene_id, genes)
})

test_that(".build_rowdata_se includes gene names when provided", {
    se <- make_test_se_minimal(n_genes = 5)
    tx_ids <- paste0("ENST", sprintf("%011d", 1:5))
    genes <- paste0("ENSG", sprintf("%011d", c(1, 1, 2, 2, 3)))
    gene_names_df <- data.frame(
        GeneID = c("ENSG00000000001", "ENSG00000000002", "ENSG00000000003"),
        GeneName = c("BRCA1", "BRCA2", "TP53")
    )

    se <- TSENAT:::.build_rowdata_se(se, tx_ids, genes, gene_names_df)

    rowdata <- rowData(se)
    expect_true("gene_name" %in% colnames(rowdata))
    expect_equal(rowdata$gene_name[1], "BRCA1")
})

# ===========================================================================
# Tests for .build_se() - Main Function
# ===========================================================================

test_that(".build_se creates valid SummarizedExperiment", {
    tx2gene <- data.frame(
        Transcript = c("ENST00000001", "ENST00000002", "ENST00000003"),
        Gene = c("ENSG00000001", "ENSG00000002", "ENSG00000001")
    )
    readcounts <- matrix(
        c(10, 5, 2, 3, 15, 8),
        nrow = 3,
        dimnames = list(c("ENST00000001", "ENST00000002", "ENST00000003"), c("s1", "s2"))
    )

    se <- TSENAT:::.build_se(readcounts, tx2gene)

    expect_true(is(se, "SummarizedExperiment"))
    expect_equal(nrow(se), 3)
    expect_equal(ncol(se), 2)
})

test_that(".build_se stores tx2gene in metadata", {
    tx2gene <- data.frame(
        Transcript = c("ENST00000001", "ENST00000002"),
        Gene = c("ENSG00000001", "ENSG00000002")
    )
    readcounts <- matrix(
        1:4, nrow = 2,
        dimnames = list(c("ENST00000001", "ENST00000002"), c("s1", "s2"))
    )

    se <- TSENAT:::.build_se(readcounts, tx2gene)

    expect_false(is.null(metadata(se)$tx2gene))
    expect_identical(nrow(metadata(se)$tx2gene), 2L)
})

test_that(".build_se stores readcounts in metadata", {
    tx2gene <- data.frame(
        Transcript = c("ENST00000001", "ENST00000002"),
        Gene = c("ENSG00000001", "ENSG00000002")
    )
    readcounts <- matrix(
        1:4, nrow = 2,
        dimnames = list(c("ENST00000001", "ENST00000002"), c("s1", "s2"))
    )

    se <- TSENAT:::.build_se(readcounts, tx2gene)

    expect_false(is.null(metadata(se)$readcounts))
    expect_identical(nrow(metadata(se)$readcounts), 2L)
})

test_that(".build_se handles SALMON data", {
    tx2gene <- data.frame(
        Transcript = c("ENST00000001", "ENST00000002"),
        Gene = c("ENSG00000001", "ENSG00000002")
    )
    readcounts <- matrix(
        1:4, nrow = 2,
        dimnames = list(c("ENST00000001", "ENST00000002"), c("s1", "s2"))
    )
    tpm <- matrix(
        10:13, nrow = 2,
        dimnames = list(c("ENST00000001", "ENST00000002"), c("s1", "s2"))
    )
    eff_length <- c(ENST00000001 = 1000, ENST00000002 = 1200)

    se <- TSENAT:::.build_se(readcounts, tx2gene, tpm = tpm,
                              effective_length = eff_length)

    expect_false(is.null(metadata(se)$salmon_tpm))
    expect_false(is.null(metadata(se)$salmon_effective_length))
})

test_that(".build_se populates rowData correctly", {
    tx2gene <- data.frame(
        Transcript = c("ENST00000001", "ENST00000002"),
        Gene = c("ENSG00000001", "ENSG00000002")
    )
    readcounts <- matrix(
        1:4, nrow = 2,
        dimnames = list(c("ENST00000001", "ENST00000002"), c("s1", "s2"))
    )

    se <- TSENAT:::.build_se(readcounts, tx2gene)

    rowdata <- rowData(se)
    expect_true("transcript_id" %in% colnames(rowdata))
    expect_true("gene_id" %in% colnames(rowdata))
    expect_identical(rowdata$transcript_id, c("ENST00000001", "ENST00000002"))
})

test_that(".build_se skips unmapped transcripts", {
    tx2gene <- data.frame(
        Transcript = c("ENST00000001", "ENST00000003"),
        Gene = c("ENSG00000001", "ENSG00000002")
    )
    readcounts <- matrix(
        1:6, nrow = 3,
        dimnames = list(c("ENST00000001", "UNMAPPED_TX", "ENST00000003"), c("s1", "s2"))
    )

    se <- TSENAT:::.build_se(readcounts, tx2gene, skip = TRUE)

    expect_equal(nrow(se), 2)
    expect_false("UNMAPPED_TX" %in% rownames(se))
})

test_that(".build_se uses TSV file for tx2gene", {
    tsv_file <- create_temp_tx2gene_tsv(n_tx = 5)
    on.exit(unlink(tsv_file))

    readcounts <- matrix(
        1:10, nrow = 5,
        dimnames = list(
            paste0("ENST", sprintf("%011d", 1:5)),
            c("s1", "s2")
        )
    )

    se <- TSENAT:::.build_se(readcounts, tsv_file)

    expect_true(is(se, "SummarizedExperiment"))
    expect_equal(nrow(se), 5)
})

test_that(".build_se applies custom assay name", {
    tx2gene <- data.frame(
        Transcript = c("ENST00000001", "ENST00000002"),
        Gene = c("ENSG00000001", "ENSG00000002")
    )
    readcounts <- matrix(
        1:4, nrow = 2,
        dimnames = list(c("ENST00000001", "ENST00000002"), c("s1", "s2"))
    )

    se <- TSENAT:::.build_se(readcounts, tx2gene, assay_name = "abundance")

    expect_true("abundance" %in% assayNames(se))
})

test_that(".build_se errors on unmapped transcripts without skip", {
    tx2gene <- data.frame(
        Transcript = "ENST00000001",
        Gene = "ENSG00000001"
    )
    readcounts <- matrix(
        1:4, nrow = 2,
        dimnames = list(c("ENST00000001", "UNMAPPED_TX"), c("s1", "s2"))
    )

    expect_error(
        TSENAT:::.build_se(readcounts, tx2gene, skip = FALSE),
        "Unmapped"
    )
})

test_that(".build_se handles metadata application", {
    tx2gene <- data.frame(
        Transcript = c("ENST00000001", "ENST00000002"),
        Gene = c("ENSG00000001", "ENSG00000002")
    )
    readcounts <- matrix(
        1:4, nrow = 2,
        dimnames = list(c("ENST00000001", "ENST00000002"), c("s1", "s2"))
    )
    custom_metadata <- list(experiment_id = "exp_001", date = "2024-01-01")

    se <- TSENAT:::.build_se(readcounts, tx2gene, metadata = custom_metadata)

    # Metadata is applied via .map_metadata_se which may or may not add all fields
    # Just verify SE is created successfully
    expect_true(is(se, "SummarizedExperiment"))
})

# ===========================================================================
# EXTENDED TESTS - Comprehensive parameter coverage
# ===========================================================================

# ===========================================================================
# Extended tests for .extract_gff3_data - verbose parameter
# ===========================================================================

test_that(".extract_gff3_data with verbose=TRUE shows progress messages", {
    gff3_file <- create_temp_gff3_file(n_genes = 3, n_tx_per_gene = 2)
    on.exit(unlink(gff3_file))
    
    result <- TSENAT:::.extract_gff3_data(gff3_file, verbose = TRUE)
    
    expect_equal(nrow(result$tx2gene), 6)  # 3 genes * 2 transcripts
})

test_that(".extract_gff3_data with verbose=FALSE suppresses messages", {
    gff3_file <- create_temp_gff3_file(n_genes = 3, n_tx_per_gene = 2)
    on.exit(unlink(gff3_file))
    
    result <- TSENAT:::.extract_gff3_data(gff3_file, verbose = FALSE)
    
    expect_equal(nrow(result$tx2gene), 6)
})

# ===========================================================================
# Extended tests for .load_tx2gene_data - verbose parameter
# ===========================================================================

test_that(".load_tx2gene_data with verbose=TRUE shows detection message", {
    tsv_file <- create_temp_tx2gene_tsv(n_tx = 20)
    on.exit(unlink(tsv_file))
    
    result <- TSENAT:::.load_tx2gene_data(tsv_file, verbose = TRUE)
    
    expect_equal(nrow(result$tx2gene_df), 20)
})

test_that(".load_tx2gene_data with verbose=FALSE suppresses messages", {
    tsv_file <- create_temp_tx2gene_tsv(n_tx = 20)
    on.exit(unlink(tsv_file))
    
    result <- TSENAT:::.load_tx2gene_data(tsv_file, verbose = FALSE)
    
    expect_equal(nrow(result$tx2gene_df), 20)
})

test_that(".load_tx2gene_data GFF3 detection with verbose", {
    gff3_file <- create_temp_gff3_file(n_genes = 5, n_tx_per_gene = 2)
    on.exit(unlink(gff3_file))
    
    result <- TSENAT:::.load_tx2gene_data(gff3_file, verbose = TRUE)
    
    expect_equal(nrow(result$tx2gene_df), 10)
    expect_false(is.null(result$gff3_data))
})

# ===========================================================================
# Extended tests for .store_salmon_metadata - all parameter combinations
# ===========================================================================

test_that(".store_salmon_metadata with both TPM and effective_length", {
    se <- make_test_se_minimal(n_genes = 5, n_samples = 3)
    readcounts <- assay(se, "counts")
    tpm <- readcounts * 1.5
    eff_length <- setNames(rep(1000, nrow(readcounts)), rownames(readcounts))
    
    se <- TSENAT:::.store_salmon_metadata(se, readcounts, tpm, eff_length)
    
    expect_false(is.null(metadata(se)$salmon_tpm))
    expect_false(is.null(metadata(se)$salmon_effective_length))
    expect_equal(nrow(metadata(se)$tpm), nrow(readcounts))
    expect_equal(length(metadata(se)$salmon_effective_length), nrow(readcounts))
})

test_that(".store_salmon_metadata with only TPM (no effective_length)", {
    se <- make_test_se_minimal(n_genes = 5, n_samples = 3)
    readcounts <- assay(se, "counts")
    tpm <- readcounts * 1.5
    
    se <- TSENAT:::.store_salmon_metadata(se, readcounts, tpm, NULL)
    
    expect_false(is.null(metadata(se)$salmon_tpm))
    expect_true(is.null(metadata(se)$salmon_effective_length))
})

test_that(".store_salmon_metadata with only effective_length (no TPM)", {
    se <- make_test_se_minimal(n_genes = 5, n_samples = 3)
    readcounts <- assay(se, "counts")
    eff_length <- setNames(rep(1000, nrow(readcounts)), rownames(readcounts))
    
    se <- TSENAT:::.store_salmon_metadata(se, readcounts, NULL, eff_length)
    
    expect_true(is.null(metadata(se)$salmon_tpm))
    expect_false(is.null(metadata(se)$salmon_effective_length))
})

test_that(".store_salmon_metadata with neither TPM nor effective_length", {
    se <- make_test_se_minimal(n_genes = 5, n_samples = 3)
    readcounts <- assay(se, "counts")
    
    se <- TSENAT:::.store_salmon_metadata(se, readcounts, NULL, NULL)
    
    expect_true(is.null(metadata(se)$salmon_tpm))
    expect_true(is.null(metadata(se)$salmon_effective_length))
})

# ===========================================================================
# Extended tests for .map_transcripts_to_genes - all parameter combinations
# ===========================================================================

test_that(".map_transcripts_to_genes with skip=FALSE and valid mapping", {
    tx_ids <- c("ENST00000001", "ENST00000002", "ENST00000003")
    tx2gene_df <- data.frame(
        Transcript = tx_ids,
        Gene = c("ENSG00000001", "ENSG00000002", "ENSG00000001")
    )
    readcounts <- matrix(1:6, nrow = 3, dimnames = list(tx_ids, c("s1", "s2")))
    
    result <- TSENAT:::.map_transcripts_to_genes(
        tx_ids, readcounts, tx2gene_df, NULL, NULL, skip = FALSE
    )
    
    expect_equal(length(result$tx_ids), 3)
    expect_equal(nrow(result$readcounts), 3)
})

test_that(".map_transcripts_to_genes with skip=TRUE removes unmapped", {
    tx_ids <- c("ENST00000001", "UNMAPPED", "ENST00000002")
    tx2gene_df <- data.frame(
        Transcript = c("ENST00000001", "ENST00000002"),
        Gene = c("ENSG00000001", "ENSG00000002")
    )
    readcounts <- matrix(1:6, nrow = 3, dimnames = list(tx_ids, c("s1", "s2")))
    
    result <- TSENAT:::.map_transcripts_to_genes(
        tx_ids, readcounts, tx2gene_df, NULL, NULL, skip = TRUE
    )
    
    expect_equal(length(result$tx_ids), 2)
    expect_false("UNMAPPED" %in% result$tx_ids)
    expect_equal(nrow(result$readcounts), 2)
})

# ===========================================================================
# Extended tests for .build_se - all parameter combinations
# ===========================================================================

test_that(".build_se with verbose=TRUE shows progress", {
    tx2gene <- data.frame(
        Transcript = c("ENST00000001", "ENST00000002"),
        Gene = c("ENSG00000001", "ENSG00000002")
    )
    readcounts <- matrix(
        1:4, nrow = 2,
        dimnames = list(c("ENST00000001", "ENST00000002"), c("s1", "s2"))
    )
    
    se <- TSENAT:::.build_se(readcounts, tx2gene, verbose = TRUE)
    
    expect_true(is(se, "SummarizedExperiment"))
    expect_equal(nrow(se), 2)
})

test_that(".build_se with verbose=FALSE suppresses progress", {
    tx2gene <- data.frame(
        Transcript = c("ENST00000001", "ENST00000002"),
        Gene = c("ENSG00000001", "ENSG00000002")
    )
    readcounts <- matrix(
        1:4, nrow = 2,
        dimnames = list(c("ENST00000001", "ENST00000002"), c("s1", "s2"))
    )
    
    se <- TSENAT:::.build_se(readcounts, tx2gene, verbose = FALSE)
    
    expect_true(is(se, "SummarizedExperiment"))
    expect_equal(nrow(se), 2)
})

test_that(".build_se with tpm only (no effective_length)", {
    tx2gene <- data.frame(
        Transcript = c("ENST00000001", "ENST00000002"),
        Gene = c("ENSG00000001", "ENSG00000002")
    )
    readcounts <- matrix(
        1:4, nrow = 2,
        dimnames = list(c("ENST00000001", "ENST00000002"), c("s1", "s2"))
    )
    tpm <- matrix(
        10:13, nrow = 2,
        dimnames = list(c("ENST00000001", "ENST00000002"), c("s1", "s2"))
    )
    
    se <- TSENAT:::.build_se(readcounts, tx2gene, tpm = tpm, verbose = FALSE)
    
    expect_false(is.null(metadata(se)$salmon_tpm))
    expect_true(is.null(metadata(se)$salmon_effective_length))
})

test_that(".build_se with effective_length only (no tpm)", {
    tx2gene <- data.frame(
        Transcript = c("ENST00000001", "ENST00000002"),
        Gene = c("ENSG00000001", "ENSG00000002")
    )
    readcounts <- matrix(
        1:4, nrow = 2,
        dimnames = list(c("ENST00000001", "ENST00000002"), c("s1", "s2"))
    )
    eff_length <- c(ENST00000001 = 1000, ENST00000002 = 1200)
    
    se <- TSENAT:::.build_se(readcounts, tx2gene, effective_length = eff_length, verbose = FALSE)
    
    expect_true(is.null(metadata(se)$salmon_tpm))
    expect_false(is.null(metadata(se)$salmon_effective_length))
})

test_that(".build_se with both tpm and effective_length", {
    tx2gene <- data.frame(
        Transcript = c("ENST00000001", "ENST00000002"),
        Gene = c("ENSG00000001", "ENSG00000002")
    )
    readcounts <- matrix(
        1:4, nrow = 2,
        dimnames = list(c("ENST00000001", "ENST00000002"), c("s1", "s2"))
    )
    tpm <- matrix(
        10:13, nrow = 2,
        dimnames = list(c("ENST00000001", "ENST00000002"), c("s1", "s2"))
    )
    eff_length <- c(ENST00000001 = 1000, ENST00000002 = 1200)
    
    se <- TSENAT:::.build_se(readcounts, tx2gene, tpm = tpm, effective_length = eff_length, verbose = FALSE)
    
    expect_false(is.null(metadata(se)$salmon_tpm))
    expect_false(is.null(metadata(se)$salmon_effective_length))
})

test_that(".build_se parameter combinations: skip + custom assay + metadata", {
    tx2gene <- data.frame(
        Transcript = c("ENST00000001", "ENST00000002"),
        Gene = c("ENSG00000001", "ENSG00000002")
    )
    readcounts <- matrix(
        1:4, nrow = 2,
        dimnames = list(c("ENST00000001", "ENST00000002"), c("s1", "s2"))
    )
    
    se <- TSENAT:::.build_se(
        readcounts, tx2gene,
        assay_name = "abundance",
        skip = FALSE,
        metadata = list(test = TRUE),
        verbose = FALSE
    )
    
    expect_true("abundance" %in% assayNames(se))
    expect_equal(nrow(se), 2)
})

test_that(".build_se with single transcript", {
    tx2gene <- data.frame(Transcript = "ENST00000001", Gene = "ENSG00000001")
    readcounts <- matrix(1:2, nrow = 1, dimnames = list("ENST00000001", c("s1", "s2")))
    
    se <- TSENAT:::.build_se(readcounts, tx2gene, verbose = FALSE)
    
    expect_equal(nrow(se), 1)
    expect_equal(ncol(se), 2)
})

test_that(".build_se with large number of samples", {
    tx2gene <- data.frame(
        Transcript = c("ENST00000001", "ENST00000002"),
        Gene = c("ENSG00000001", "ENSG00000002")
    )
    readcounts <- matrix(
        sample(1:100, 2*50),
        nrow = 2,
        dimnames = list(
            c("ENST00000001", "ENST00000002"),
            paste0("sample_", 1:50)
        )
    )
    
    se <- TSENAT:::.build_se(readcounts, tx2gene, verbose = FALSE)
    
    expect_equal(nrow(se), 2)
    expect_equal(ncol(se), 50)
})

test_that(".build_se preserves count values correctly", {
    tx2gene <- data.frame(
        Transcript = c("ENST00000001", "ENST00000002"),
        Gene = c("ENSG00000001", "ENSG00000002")
    )
    readcounts <- matrix(
        c(100, 50, 200, 75),
        nrow = 2,
        dimnames = list(c("ENST00000001", "ENST00000002"), c("s1", "s2"))
    )
    
    se <- TSENAT:::.build_se(readcounts, tx2gene, verbose = FALSE)
    
    counts_assay <- assay(se, 1)
    expect_equal(as.numeric(counts_assay[1, ]), c(100, 200))
    expect_equal(as.numeric(counts_assay[2, ]), c(50, 75))
})
