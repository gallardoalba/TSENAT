# Ensure package is loaded so helpers and testthat context are available
if (!requireNamespace("pkgload", quietly = TRUE)) install.packages("pkgload", repos = 'https://cloud.r-project.org')
if (!requireNamespace("testthat", quietly = TRUE)) install.packages("testthat", repos = 'https://cloud.r-project.org')
suppressPackageStartupMessages(pkgload::load_all())

test_files <- list.files('tests/testthat', pattern = '\\.[Rr]$', full.names = TRUE)
fc <- covr::file_coverage('R/calc_lm_helpers.R', test_files)
cat('calc_lm_helpers.R coverage:', covr::percent_coverage(fc), '\n')
cov_df <- as.data.frame(fc)
unc <- cov_df[cov_df$value == 0, ]
cat('uncovered lines:', nrow(unc), '\n')
if (nrow(unc) > 0) print(head(unc, 200)) else cat('None\n')
