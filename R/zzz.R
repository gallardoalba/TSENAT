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
  registerS3method("print", "gtable", print.gtable)
}
