context("Diversity Calculation: Preserved Counts Assay")

test_that("calculate_diversity preserves original counts assay", {
    # Create a simple SummarizedExperiment with isoform-level counts
    # 8 rows = 2 isoforms per gene (4 genes total)
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = matrix(c(100, 50, 30, 20, 40, 60, 80, 10,
                                        90, 55, 35, 25, 45, 65, 75, 15), nrow = 8, ncol = 2)),
        rowData = data.frame(gene_name = rep(paste0("Gene", 1:4), each = 2)),
        colData = data.frame(sample = c("S1", "S2"), row.names = c("S1", "S2"))
    )
    # Isoform identifiers: Iso1a, Iso1b, Iso2a, Iso2b, Iso3a, Iso3b, Iso4a, Iso4b
    rownames(se) <- paste0("Iso", rep(1:4, each = 2), c("a", "b"))
    genes <- rep(paste0("Gene", 1:4), each = 2)
    
    # Apply calculate_diversity with genes parameter
    div_se <- calculate_diversity(se, genes = genes, q = 2, norm = TRUE)
    
    # Check that counts assay is preserved
    expect_true("counts" %in% SummarizedExperiment::assayNames(div_se))
    # Counts should be aggregated to gene level (4 genes x 2 samples)
    counts_out <- SummarizedExperiment::assay(div_se, "counts")
    expect_equal(nrow(counts_out), 4)  # 4 genes
    expect_equal(ncol(counts_out), 2)  # 2 samples
})

test_that("calculate_diversity preserves counts with different normalization", {
    # Create a SummarizedExperiment with isoform-level data
    # 20 rows = 2 isoforms per gene (10 genes total)
    set.seed(789)
    counts_mat <- matrix(rpois(80, 20), nrow = 20, ncol = 4)
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = counts_mat),
        colData = data.frame(sample = paste0("S", 1:4), row.names = paste0("S", 1:4))
    )
    rownames(se) <- paste0("Iso", rep(1:10, each = 2), c("a", "b"))
    genes <- rep(paste0("Gene", 1:10), each = 2)
    
    # Calculate diversity with norm = FALSE
    div_se_raw <- calculate_diversity(se, genes = genes, q = 1, norm = FALSE)
    
    # Check structure and counts preservation
    expect_true("counts" %in% SummarizedExperiment::assayNames(div_se_raw))
    counts_out <- SummarizedExperiment::assay(div_se_raw, "counts")
    expect_equal(nrow(counts_out), 10)  # 10 genes
    expect_equal(ncol(counts_out), 4)   # 4 samples
})

test_that("calculate_diversity preserves counts for bootstrap compatibility", {
    # Create realistic test data with isoform-level counts
    # 20 rows = 2 isoforms per gene (10 genes total)
    set.seed(123)
    counts_mat <- matrix(rpois(100, 20), nrow = 20, ncol = 5)
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = counts_mat),
        colData = data.frame(
            group = rep(c("A", "B"), c(2, 3)),
            row.names = paste0("S", 1:5)
        )
    )
    rownames(se) <- paste0("Iso", rep(1:10, each = 2), c("a", "b"))
    genes <- rep(paste0("Gene", 1:10), each = 2)
    
    # Apply calculate_diversity
    div_se <- calculate_diversity(se, genes = genes, q = 2, norm = TRUE)
    
    # Test that bootstrap works with the diversity-transformed SE
    # This verifies that counts assay is accessible and usable
    expect_no_error({
        result <- calculate_tsallis_entropy_bootstrap(
            se = div_se, 
            x = SummarizedExperiment::assay(div_se, "counts")[1, ],
            q = 2, 
            nboot = 100,
            seed = 42
        )
    })
    
    expect_is(result, "tsenat_bootstrap_ci")
    expect_true(!is.na(result$estimate))
})

test_that("calculate_diversity preserves counts for jackknife compatibility", {
    # Create realistic test data with isoform-level counts
    # 20 rows = 2 isoforms per gene (10 genes total)
    set.seed(456)
    counts_mat <- matrix(rpois(120, 25), nrow = 20, ncol = 6)
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = counts_mat),
        colData = data.frame(
            group = rep(c("Control", "Treatment"), each = 3),
            row.names = paste0("S", 1:6)
        )
    )
    rownames(se) <- paste0("Iso", rep(1:10, each = 2), c("a", "b"))
    genes <- rep(paste0("Gene", 1:10), each = 2)
    
    # Apply calculate_diversity
    div_se <- calculate_diversity(se, genes = genes, q = 1, norm = TRUE)
    
    # Test that jackknife works with the diversity-transformed SE
    # This verifies that counts assay is accessible and usable
    expect_no_error({
        result <- jackknife_entropy_outliers(
            x = SummarizedExperiment::assay(div_se, "counts")[1, ],
            q = 1,
            norm = TRUE
        )
    })
    
    expect_is(result, "tsenat_jackknife")
    expect_true(!is.na(result$estimate))
    expect_equal(length(result$jackknife_estimates), 6)  # Leave-one-out for 6 samples
})

test_that("diversity assay and counts assay coexist without conflict", {
    # Create test data with isoform-level counts
    # 4 rows = 2 isoforms per gene (2 genes total)
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = matrix(c(100, 50, 80, 20, 60, 40, 70, 30), nrow = 4, ncol = 2)),
        colData = data.frame(sample = c("S1", "S2"), row.names = c("S1", "S2"))
    )
    rownames(se) <- c("Iso1a", "Iso1b", "Iso2a", "Iso2b")
    genes <- c("Gene1", "Gene1", "Gene2", "Gene2")
    
    # Calculate diversity
    div_se <- calculate_diversity(se, genes = genes, q = 2, norm = TRUE)
    
    # Check that both diversity and counts are present
    assay_names <- SummarizedExperiment::assayNames(div_se)
    expect_true("diversity" %in% assay_names)
    expect_true("counts" %in% assay_names)
    
    # Both should have the same dimensions (gene-level for output)
    diversity_assay <- SummarizedExperiment::assay(div_se, "diversity")
    counts_assay <- SummarizedExperiment::assay(div_se, "counts")
    expect_equal(dim(diversity_assay), dim(counts_assay))
    expect_equal(nrow(diversity_assay), 2)  # 2 genes
})

test_that("Hill numbers preserves counts assay", {
    # Create test data with isoform-level counts
    # 10 rows = 2 isoforms per gene (5 genes total)
    set.seed(999)
    counts_mat <- matrix(rpois(80, 15), nrow = 10, ncol = 8)
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = counts_mat),
        colData = data.frame(sample = paste0("S", 1:8), row.names = paste0("S", 1:8))
    )
    rownames(se) <- paste0("Iso", rep(1:5, each = 2), c("a", "b"))
    genes <- rep(paste0("Gene", 1:5), each = 2)
    
    # Calculate Hill numbers (D instead of S)
    hill_se <- calculate_diversity(se, genes = genes, q = 1.5, what = "D", norm = TRUE)
    
    # Check that counts assay is preserved
    expect_true("counts" %in% SummarizedExperiment::assayNames(hill_se))
    expect_true("hill" %in% SummarizedExperiment::assayNames(hill_se))
    
    # Verify counts structure
    counts_out <- SummarizedExperiment::assay(hill_se, "counts")
    expect_equal(nrow(counts_out), 5)  # 5 genes
    expect_equal(ncol(counts_out), 8)  # 8 samples
})

test_that("metadata still includes readcounts reference for backward compatibility", {
    # Create test data with isoform-level counts
    # 10 rows = 2 isoforms per gene (5 genes total)
    set.seed(321)
    counts_mat <- matrix(rpois(60, 20), nrow = 10, ncol = 6)
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = counts_mat),
        colData = data.frame(sample = paste0("S", 1:6), row.names = paste0("S", 1:6))
    )
    rownames(se) <- paste0("Iso", rep(1:5, each = 2), c("a", "b"))
    genes <- rep(paste0("Gene", 1:5), each = 2)
    
    # Calculate diversity
    div_se <- calculate_diversity(se, genes = genes, q = 2, norm = TRUE)
    
    # Check metadata still has readcounts for backward compatibility
    metadata <- S4Vectors::metadata(div_se)
    expect_true("readcounts" %in% names(metadata))
    expect_true(!is.null(metadata$readcounts))
})
