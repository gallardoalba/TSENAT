#!/usr/bin/env Rscript
# Script to run coverage, generate HTML report, and optionally upload to Codecov
args <- commandArgs(trailingOnly = TRUE)
upload <- any(grepl("--upload", args))
inst <- function(pkgs){
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if(length(missing)) install.packages(missing, repos = 'https://cloud.r-project.org')
}
inst(c('covr','DT','htmltools'))
library(covr)
cov <- covr::package_coverage()
pc <- covr::percent_coverage(cov)
cat(sprintf("Coverage: %.4f%%\n", pc))
# write report to coverage/coverage.html
dir.create('coverage', showWarnings = FALSE)
report_file <- file.path('coverage', 'coverage.html')
if (requireNamespace("DT", quietly = TRUE) && requireNamespace("htmltools", quietly = TRUE)) {
  covr::report(cov, file = report_file)
  cat("Wrote report:", report_file, "\n")
} else {
  cat("CI: Skipping HTML coverage report: 'DT' and/or 'htmltools' not available.\n")
}
if(upload){
  if(identical(Sys.getenv('CODECOV_TOKEN', unset = ''), '')){
    message('No CODECOV_TOKEN set; attempting upload may fail for private repos')
  }
  covr::codecov(coverage = cov)
  cat('Uploaded coverage via covr::codecov()\n')
}
