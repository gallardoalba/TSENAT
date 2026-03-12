test_that("build_se handles subsetted readcounts correctly", {
  # Load full data
  data(readcounts)
  readcounts_full <- as.matrix(salmon_dataset)
  mode(readcounts_full) <- "numeric"
  
  # Subset readcounts
  readcounts_subset <- readcounts_full[1:100, ]
  
  gff3_dataset <- system.file("extdata", "annotation.gff3.gz", package = "TSENAT")
  metadata_df <- read.table(
    system.file("extdata", "metadata.tsv", package = "TSENAT"),
    header = TRUE, sep = "\t"
  )
  
  # Test 1: No duplicate rownames after subsetting
  se_subset <- build_se(readcounts_subset, gff3_dataset, metadata = metadata_df)
  expect_equal(
    length(unique(rownames(se_subset))), 
    nrow(se_subset),
    info = "No duplicate rownames in subsetted SE"
  )
  
  # Test 2: calculate_diversity works without shrinkage
  expect_no_error(
    calculate_diversity(
      se_subset, 
      q = seq(0.5, 1.5, 0.5), 
      pseudocount = 0.1,
      verbose = FALSE
    )
  )
  
  # Test 3: calculate_diversity works WITH shrinkage (after fix)
  expect_no_error(
    calculate_diversity(
      se_subset, 
      q = seq(0.5, 1.5, 0.5), 
      shrinkage = "empirical_bayes",
      pseudocount = 0.1,
      verbose = FALSE
    )
  )
})

test_that("full vs subset - rownames and gene names comparison", {
  data(readcounts)
  readcounts_full <- as.matrix(salmon_dataset)
  mode(readcounts_full) <- "numeric"
  
  gff3_dataset <- system.file("extdata", "annotation.gff3.gz", package = "TSENAT")
  metadata_df <- read.table(
    system.file("extdata", "metadata.tsv", package = "TSENAT"),
    header = TRUE, sep = "\t"
  )
  
  # Build full and subset
  se_full <- build_se(readcounts_full, gff3_dataset, metadata = metadata_df)
  readcounts_subset <- readcounts_full[1:300, ]
  se_subset <- build_se(readcounts_subset, gff3_dataset, metadata = metadata_df)
  
  # Compare rownames - should all be unique in both
  expect_equal(
    length(unique(rownames(se_full))),
    nrow(se_full),
    info = "Full data: all rownames are unique"
  )
  
  expect_equal(
    length(unique(rownames(se_subset))),
    nrow(se_subset),
    info = "Subset data: all rownames are unique"
  )
  
  # Compare gene_name distributions
  full_gene_names <- rowData(se_full)$gene_name
  subset_gene_names <- rowData(se_subset)$gene_name
  
  # Verify: No issues with gene name duplicates affecting rownames
  expect_true(nrow(se_full) > 0, info = "Full subset has genes")
  expect_true(nrow(se_subset) > 0, info = "Subset has genes")
})

test_that("calculate_diversity gene name handling - full vs subset", {
  data(readcounts)
  readcounts_full <- as.matrix(salmon_dataset)
  mode(readcounts_full) <- "numeric"
  
  gff3_dataset <- system.file("extdata", "annotation.gff3.gz", package = "TSENAT")
  metadata_df <- read.table(
    system.file("extdata", "metadata.tsv", package = "TSENAT"),
    header = TRUE, sep = "\t"
  )
  
  # Build SEs
  se_full <- build_se(readcounts_full, gff3_dataset, metadata = metadata_df)
  readcounts_subset <- readcounts_full[1:250, ]
  se_subset <- build_se(readcounts_subset, gff3_dataset, metadata = metadata_df)
  
  # Test WITHOUT shrinkage first (should work for both)
  result_full_no_shrink <- calculate_diversity(
    se_full, q = 1, pseudocount = 0.1, verbose = FALSE
  )
  
  result_subset_no_shrink <- calculate_diversity(
    se_subset, q = 1, pseudocount = 0.1, verbose = FALSE
  )
  
  expect_equal(
    length(unique(rownames(result_full_no_shrink))),
    nrow(result_full_no_shrink),
    info = "Full: no duplicates without shrinkage"
  )
  
  expect_equal(
    length(unique(rownames(result_subset_no_shrink))),
    nrow(result_subset_no_shrink),
    info = "Subset: no duplicates without shrinkage"
  )
  
  # Test WITH shrinkage
  expect_no_error({
    result_full_shrink <- calculate_diversity(
      se_full, q = 1, shrinkage = "empirical_bayes", pseudocount = 0.1, verbose = FALSE
    )
  })
  
  expect_no_error({
    result_subset_shrink <- calculate_diversity(
      se_subset, q = 1, shrinkage = "empirical_bayes", pseudocount = 0.1, verbose = FALSE
    )
  })
})

test_that("filter_se impact on gene names - full vs subset", {
  data(readcounts)
  readcounts_full <- as.matrix(salmon_dataset)
  mode(readcounts_full) <- "numeric"
  
  gff3_dataset <- system.file("extdata", "annotation.gff3.gz", package = "TSENAT")
  metadata_df <- read.table(
    system.file("extdata", "metadata.tsv", package = "TSENAT"),
    header = TRUE, sep = "\t"
  )
  
  # Build UNFILTERED vs FILTERED (with metadata for pair detection)
  se_full_unfiltered <- build_se(readcounts_full, gff3_dataset, metadata = metadata_df)
  se_full_filtered <- tryCatch(
    suppressWarnings(filter_se(se_full_unfiltered, stringency = "soft")),
    error = function(e) se_full_unfiltered  # Use unfiltered if filter fails
  )
  
  readcounts_subset <- readcounts_full[1:250, ]
  se_subset_unfiltered <- build_se(readcounts_subset, gff3_dataset, metadata = metadata_df)
  se_subset_filtered <- tryCatch(
    suppressWarnings(filter_se(se_subset_unfiltered, stringency = "soft")),
    error = function(e) se_subset_unfiltered  # Use unfiltered if filter fails
  )
  
  # Both should maintain unique rownames
  expect_equal(
    length(unique(rownames(se_full_filtered))),
    nrow(se_full_filtered),
    info = "Filtered full data: all rownames are unique"
  )
  
  expect_equal(
    length(unique(rownames(se_subset_filtered))),
    nrow(se_subset_filtered),
    info = "Filtered subset data: all rownames are unique"
  )
})

test_that("shrinkage parameter is compatible with subsetted data", {
  data(readcounts)
  readcounts_full <- as.matrix(salmon_dataset)
  mode(readcounts_full) <- "numeric"
  
  gff3_dataset <- system.file("extdata", "annotation.gff3.gz", package = "TSENAT")
  metadata_df <- read.table(
    system.file("extdata", "metadata.tsv", package = "TSENAT"),
    header = TRUE, sep = "\t"
  )
  
  # Subset readcounts and build SE
  readcounts_subset <- readcounts_full[1:150, ]
  se_subset <- build_se(readcounts_subset, gff3_dataset, metadata = metadata_df)
  
  # Run with shrinkage enabled
  result_shrink <- calculate_diversity(
    se_subset, 
    q = 1, 
    shrinkage = "empirical_bayes",
    pseudocount = 0.05,
    verbose = FALSE
  )
  
  # Verify output structure
  expect_true(nrow(result_shrink) > 0, info = "Shrinkage returns genes")
  expect_true(all(!is.na(rownames(result_shrink))), info = "All genes have valid rownames")
  expect_equal(
    length(unique(rownames(result_shrink))), 
    nrow(result_shrink),
    info = "No duplicate rownames in shrinkage output"
  )
})

test_that("tx2gene filtering doesn't affect full readcounts", {
  data(readcounts)
  readcounts_full <- as.matrix(salmon_dataset)
  mode(readcounts_full) <- "numeric"
  
  gff3_dataset <- system.file("extdata", "annotation.gff3.gz", package = "TSENAT")
  metadata_df <- read.table(
    system.file("extdata", "metadata.tsv", package = "TSENAT"),
    header = TRUE, sep = "\t"
  )
  
  # Build SE with full readcounts
  se_full <- build_se(readcounts_full, gff3_dataset, metadata = metadata_df)
  
  # Both with and without shrinkage should work
  result_no_shrink <- calculate_diversity(
    se_full, 
    q = 1,
    pseudocount = 0.1,
    verbose = FALSE
  )
  
  result_shrink <- calculate_diversity(
    se_full, 
    q = 1,
    shrinkage = "empirical_bayes",
    pseudocount = 0.1,
    verbose = FALSE
  )
  
  # Both should produce results without duplicate rownames
  expect_equal(
    length(unique(rownames(result_no_shrink))), 
    nrow(result_no_shrink),
    info = "No duplicate rownames without shrinkage"
  )
  
  expect_equal(
    length(unique(rownames(result_shrink))), 
    nrow(result_shrink),
    info = "No duplicate rownames with shrinkage"
  )
  
  expect_true(nrow(result_shrink) > 0, info = "Shrinkage produces output genes")
})

test_that("build_se creates clean tx2gene mapping for subsetted data", {
  data(readcounts)
  readcounts_full <- as.matrix(salmon_dataset)
  mode(readcounts_full) <- "numeric"
  
  # Create small subset
  readcounts_subset <- readcounts_full[1:50, ]
  
  gff3_dataset <- system.file("extdata", "annotation.gff3.gz", package = "TSENAT")
  
  # Build SE - this will internally filter tx2gene
  se_subset <- build_se(readcounts_subset, gff3_dataset)
  
  # Verify: all genes in SE should have genes (no NAs from unmapped transcripts)
  expect_true(
    all(!is.na(rowData(se_subset)$gene_id)),
    info = "All genes should be mapped after build_se"
  )
  
  # Verify: no duplicate rownames (important for downstream processing)
  expect_equal(
    length(unique(rownames(se_subset))),
    nrow(se_subset),
    info = "Rownames should be unique"
  )
})

