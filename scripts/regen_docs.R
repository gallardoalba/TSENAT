if (!requireNamespace("roxygen2", quietly = TRUE)) {
  install.packages("roxygen2", repos = "https://cran.rstudio.com")
}
message("Running roxygen2::roxygenize('.') to regenerate man/ files...")
roxygen2::roxygenize('.')
message('Done.')
