files <- c("R/diversity_functions.R","R/calculate_method.R","R/calculate_diversity.R","R/calculate_difference.R","R/infer_sample_groups.R","R/generate_plots.R")
for (f in files) {
  cat("Sourcing:", f, "\n")
  tryCatch({
    source(f)
    cat("OK\n")
  }, error = function(e) {
    cat("ERROR while sourcing:", f, "\n")
    message(e)
    quit(status = 1)
  })
}
cat("All sourced\n")
