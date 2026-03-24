#!/usr/bin/env Rscript
# Regenerate Rd files from roxygen comments
roxygen2::roxygenise(roclets = "rd")
cat("Rd file regeneration completed successfully!\n")
