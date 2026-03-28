# Package initialization and finalization

.onLoad <- function(libname, pkgname) {
  # Ensure S3 method registration for print methods
  # This makes sure R recognizes our custom print methods
  registerS3method("print", "tsenat_bootstrap_ci", print.tsenat_bootstrap_ci)
  registerS3method("print", "tsenat_bootstrap_ci_list", print.tsenat_bootstrap_ci_list)
  registerS3method("print", "tsenat_divergence_bootstrap_ci", print.tsenat_divergence_bootstrap_ci)
  registerS3method("print", "tsenat_jackknife", print.tsenat_jackknife)
  registerS3method("print", "tsenat_jackknife_list", print.tsenat_jackknife_list)
  registerS3method("print", "rank_assumptions", print.rank_assumptions)
  registerS3method("print", "rank_correlation_ci", print.rank_correlation_ci)
  registerS3method("print", "gtable", print.gtable)
  registerS3method("summary", "tsenat_bootstrap_ci", summary.tsenat_bootstrap_ci)
  registerS3method("summary", "tsenat_divergence_bootstrap_ci", summary.tsenat_divergence_bootstrap_ci)
  
  # Wrap ggplot2's print method to suppress spurious base R plotting warnings
  # These warnings come from pretty() with degenerate axis ranges during vignette rendering
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    tryCatch({
      .ggplot_print_original <- utils::getS3method("print", "ggplot")
      if (!is.null(.ggplot_print_original)) {
        registerS3method("print", "ggplot", function(x, ...) {
          # Use withCallingHandlers to intercept and filter specific warnings
          withCallingHandlers(
            .ggplot_print_original(x, ...),
            warning = function(w) {
              # Suppress warnings from base R's pretty() when it encounters
              # degenerate axis ranges (e.g., min() on empty vector)
              # These are false positives that don't indicate actual problems
              msg <- conditionMessage(w)
              if (grepl("no finite argument|ning\u00fan argumento finito", msg, ignore.case = TRUE)) {
                invokeRestart("muffleWarning")
              }
            }
          )
        })
      }
    }, error = function(e) {
      # Silently skip if we can't register the wrapper
      NULL
    })
  }
}
