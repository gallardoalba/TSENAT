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
