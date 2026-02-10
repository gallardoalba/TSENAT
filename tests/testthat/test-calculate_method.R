context("Aggregation of methods")

test_that("Tsallis aggregation returns expected shape", {
    read_count_matrix <- rbind(
        matrix(
            rpois(
                36,
                6
            ),
            ncol = 6
        ),
        matrix(0,
            nrow = 2,
            ncol = 6
        )
    )
    colnames(read_count_matrix) <- paste0(
        "Sample",
        seq_len(ncol(read_count_matrix))
    )
    genes <- c("A", "B", "B", "C", "C", "C", "D", "D")
    qvec <- c(1, 2)

    res <- calculate_method(read_count_matrix, genes, norm = TRUE, q = qvec)
    expect_true(is.data.frame(res))
    expect_true("Gene" %in% colnames(res))
    expect_equal(ncol(res), 1 + (ncol(read_count_matrix) * length(qvec)))
})

test_that("calculate_method returns expected column names for multiple q", {
    # small synthetic transcript matrix: 4 transcripts, 2 genes (2 transcripts each), 2 samples
    mat <- matrix(c(
        10, 5,
        0, 0,
        2, 8,
        3, 7
    ), nrow = 4, byrow = TRUE)
    colnames(mat) <- c("S1", "S2")
    genes <- c("g1", "g1", "g2", "g2")
    res <- calculate_method(mat, genes, norm = TRUE, q = c(0.5, 1), what = "S")
    # expected column names: S1_q=0.5, S1_q=1, S2_q=0.5, S2_q=1
    expect_true(all(c("S1_q=0.5", "S1_q=1", "S2_q=0.5", "S2_q=1") %in% colnames(res)))
    expect_equal(res$Gene, unique(genes))
})

test_that("calculate_method handles missing sample names by producing columns nevertheless", {
    mat <- matrix(rep(1, 6), nrow = 3)
    colnames(mat) <- NULL
    genes <- letters[1:3]
    res <- calculate_method(mat, genes, q = 1, what = "S")
    # should return a data.frame with columns (Gene + one column)
    expect_true(is.data.frame(res))
    expect_true(ncol(res) >= 2)
})

context("calculate_method helpers and edge cases")

test_that("single-transcript genes produce zeros (and NAs for zero counts) via helper", {
    mat <- matrix(c(
        5, 0, # gene A, single transcript
        2, 3, # gene B, transcript 1
        4, 4 # gene B, transcript 2
    ), nrow = 3, byrow = TRUE)
    colnames(mat) <- c("S1", "S2")
    genes <- c("A", "B", "B")

    # Directly test the helper's output ordering: for each sample, q varies inner
    out_vec <- .tsenat_tsallis_row(x = mat, genes = genes, gene = "A", q = c(1, 2), norm = TRUE, what = "S")
    # ordering is: S1_q=1, S1_q=2, S2_q=1, S2_q=2
    # When normalization is TRUE and n == 1, normalization divides by zero
    # yielding NaN for the finite-count case; zero-counts produce NA
    expect_true(all(is.nan(unname(out_vec[1:2]))))
    expect_true(all(is.na(out_vec[3:4])))

    # calculate_method filters genes with non-finite entries; gene 'A' should be excluded
    res <- calculate_method(mat, genes, q = c(1, 2), what = "S")
    expect_false("A" %in% res$Gene)
})


test_that("genes with all NA counts in a sample produce NA in the helper output and are excluded from results", {
    # use NA_real_ to ensure numeric NA values (not logical)
    mat <- matrix(c(
        NA_real_, NA_real_, # transcript 1 (gene A)
        NA_real_, NA_real_ # transcript 2 (gene A)
    ), nrow = 2, byrow = TRUE)
    colnames(mat) <- c("S1", "S2")
    genes <- c("A", "A")

    out_vec <- .tsenat_tsallis_row(x = mat, genes = genes, gene = "A", q = 1, norm = TRUE, what = "S")
    expect_true(all(is.na(out_vec)))

    res <- calculate_method(mat, genes, q = 1, what = "S")
    expect_false("A" %in% res$Gene)
})


test_that("calculate_method returns correct 'D' (Hill numbers) values for multiple q", {
    # small synthetic transcript matrix: 4 transcripts, 2 genes (2 transcripts each), 2 samples
    mat <- matrix(c(
        10, 5,
        0, 0,
        2, 8,
        3, 7
    ), nrow = 4, byrow = TRUE)
    colnames(mat) <- c("S1", "S2")
    genes <- c("g1", "g1", "g2", "g2")

    qvec <- c(0.5, 1)
    res <- calculate_method(mat, genes, norm = TRUE, q = qvec, what = "D")
    # check that expected column names exist
    expect_true(all(c("S1_q=0.5", "S1_q=1", "S2_q=0.5", "S2_q=1") %in% colnames(res)))

    # for each gene and sample compare to direct calculate_tsallis_entropy(..., what = "D")
    for (g in unique(genes)) {
        idx <- which(genes == g)
        for (j in seq_len(ncol(mat))) {
            expected <- calculate_tsallis_entropy(mat[idx, j], q = qvec, norm = TRUE, what = "D")
            # expected is a numeric vector (named when length(q) > 1)
            for (k in seq_along(qvec)) {
                colname <- paste0(colnames(mat)[j], "_q=", qvec[k])
                got <- res[res$Gene == g, colname]
                # allow NA equality where appropriate, otherwise numeric equality
                if (all(is.na(expected[k]))) {
                    expect_true(is.na(got))
                } else {
                    expect_equal(as.numeric(got), unname(expected[k]))
                }
            }
        }
    }
})


test_that("helper falls back to NA when calculate_tsallis_entropy returns wrong-length vector", {
    mat <- matrix(c(1, 10), nrow = 1, byrow = TRUE)
    colnames(mat) <- c("S1", "S2")
    genes <- c("A")

    ns <- asNamespace("TSENAT")
    orig <- get("calculate_tsallis_entropy", envir = ns)
    stub <- function(x, q, norm, what) {
        # if the sample value is 1 return a scalar (wrong length) to trigger fallback
        if (length(x) == 1 && x == 1) {
            return(42)
        }
        orig(x, q = q, norm = norm, what = what)
    }
    assignInNamespace("calculate_tsallis_entropy", stub, ns = "TSENAT")
    on.exit(assignInNamespace("calculate_tsallis_entropy", orig, ns = "TSENAT"), add = TRUE)

    out_vec <- .tsenat_tsallis_row(x = mat, genes = genes, gene = "A", q = c(1, 2), norm = TRUE, what = "S")
    # S1 (first two entries) should be NA due to wrong-length return; S2 should match orig
    expect_true(all(is.na(out_vec[1:2])))
    expected2 <- orig(mat[1, 2], q = c(1, 2), norm = TRUE, what = "S")
    expect_equal(unname(out_vec[3:4]), unname(expected2))
})


test_that("helper falls back to NA when calculate_tsallis_entropy returns non-finite values", {
    mat <- matrix(c(2, 3), nrow = 1, byrow = TRUE)
    colnames(mat) <- c("S1", "S2")
    genes <- c("A")

    ns <- asNamespace("TSENAT")
    orig <- get("calculate_tsallis_entropy", envir = ns)
    stub <- function(x, q, norm, what) {
        # if the sample value is 2 return non-finite values to trigger fallback
        if (length(x) == 1 && x == 2) {
            return(c(Inf, Inf))
        }
        orig(x, q = q, norm = norm, what = what)
    }
    assignInNamespace("calculate_tsallis_entropy", stub, ns = "TSENAT")
    on.exit(assignInNamespace("calculate_tsallis_entropy", orig, ns = "TSENAT"), add = TRUE)

    out_vec <- .tsenat_tsallis_row(x = mat, genes = genes, gene = "A", q = c(1, 2), norm = TRUE, what = "S")
    expect_true(all(is.na(out_vec[1:2])))
    expected2 <- orig(mat[1, 2], q = c(1, 2), norm = TRUE, what = "S")
    expect_equal(unname(out_vec[3:4]), unname(expected2))
})


test_that("helper returns NA vector when gene has no matching transcripts", {
    mat <- matrix(c(5, 10), nrow = 1, byrow = TRUE)
    colnames(mat) <- c("S1", "S2")
    genes <- c("A")

    out_vec <- .tsenat_tsallis_row(x = mat, genes = genes, gene = "Z", q = c(1, 2), norm = TRUE, what = "S")
    expect_true(all(is.na(out_vec)))
})


test_that("helper returns empty numeric vector when there are no samples (ncol == 0)", {
    mat <- matrix(nrow = 1, ncol = 0)
    genes <- c("A")

    out_vec <- .tsenat_tsallis_row(x = mat, genes = genes, gene = "A", q = c(1, 2), norm = TRUE, what = "S")
    # accept any empty return (length 0) as the correct behavior
    expect_true(length(out_vec) == 0)
})
