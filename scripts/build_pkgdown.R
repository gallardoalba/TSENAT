#!/usr/bin/env Rscript
# Build pkgdown site with warnings suppressed to avoid spurious stack imbalance messages
options(warn = -1)
if (!requireNamespace("pkgdown", quietly = TRUE)) {
  install.packages("pkgdown", repos = "https://cloud.r-project.org")
}
pkgdown::build_site()
invisible()
