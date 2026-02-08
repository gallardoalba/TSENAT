context("map_metadata behavior")

# Resolve the path to the packaged example coldata robustly.
# When tests are run during `R CMD check` the file may not be found
# via `system.file()` in some environments; fall back to the source
# copy at `inst/extdata/coldata.tsv` when present.
coldata_path <- system.file("extdata/coldata.tsv", package = "TSENAT")
if (is.null(coldata_path) || identical(coldata_path, "") || !file.exists(coldata_path)) {
  fallback <- file.path(getwd(), "inst", "extdata", "coldata.tsv")
  if (file.exists(fallback)) {
    coldata_path <- fallback
  }
}
if (!file.exists(coldata_path)) {
  stop("Required test data not found: inst/extdata/coldata.tsv")
}

test_that("Exact match maps all samples and preserves order", {
  coldata <- utils::read.delim(coldata_path, stringsAsFactors = FALSE)
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
  coldata <- utils::read.delim(coldata_path, stringsAsFactors = FALSE)
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
  coldata <- utils::read.delim(coldata_path, stringsAsFactors = FALSE)
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
  coldata <- utils::read.delim(coldata_path, stringsAsFactors = FALSE)
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
  coldata <- utils::read.delim(coldata_path, stringsAsFactors = FALSE)
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
  coldata <- utils::read.delim(coldata_path, stringsAsFactors = FALSE)
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
  coldata <- utils::read.delim(coldata_path, stringsAsFactors = FALSE)
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
  coldata <- utils::read.delim(
    system.file("extdata/coldata.tsv", 
                 package = "TSENAT"),
    stringsAsFactors = FALSE
  )
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
  coldata <- utils::read.delim(
    system.file("extdata/coldata.tsv", 
                 package = "TSENAT"),
    stringsAsFactors = FALSE
  )
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
  coldata <- utils::read.delim(
    system.file("extdata/coldata.tsv", 
                 package = "TSENAT"),
    stringsAsFactors = FALSE
  )
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
  coldata <- utils::read.delim(coldata_path, stringsAsFactors = FALSE)
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
  expect_equal(as.character(cd$sample_base), 
    sub("_q=.*", "", colnames(SummarizedExperiment::assay(se2))))
})


test_that("paired = TRUE accepts complete pairs", {
  # construct a small paired coldata with two participants A and B
  coldata_small <- data.frame(Sample = c("A_N", "A_T", "B_N", "B_T"),
                              Condition = c("Normal", "Tumor", "Normal", "Tumor"),
                              stringsAsFactors = FALSE)
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
