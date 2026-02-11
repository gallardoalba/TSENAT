#!/usr/bin/env Rscript
library(covr)
cov <- covr::package_coverage()
uncovered <- covr::tally_coverage(cov, by = "line")
print(uncovered[uncovered$value == 0, ])
