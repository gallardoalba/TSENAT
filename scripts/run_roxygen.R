if (!requireNamespace("roxygen2", quietly = TRUE)) install.packages("roxygen2", repos = "https://cloud.r-project.org")
roxygen2::roxygenise(".")
cat("roxygenise completed\n")
