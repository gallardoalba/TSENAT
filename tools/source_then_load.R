files <- c("R/diversity_functions.R","R/calculate_method.R","R/calculate_diversity.R","R/calculate_difference.R","R/infer_sample_groups.R","R/generate_plots.R")
for (f in files) {
  cat("Sourcing:", f, "\n")
  try(source(f), silent = TRUE)
  cat("Now try loading SummarizedExperiment...\n")
  ok <- tryCatch({ library(SummarizedExperiment); TRUE }, error = function(e) { cat('Load failed after sourcing', f, '\n'); message(e); FALSE })
  if (!ok) quit(status = 1)
  cat('Load succeeded after', f, '\n')
  # detach to reset
  if ('package:SummarizedExperiment' %in% search()) detach('package:SummarizedExperiment', unload = TRUE, character.only = TRUE)
}
cat('All files sourced and SummarizedExperiment loaded successfully after each\n')
