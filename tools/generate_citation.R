#!/usr/bin/env Rscript
## Generate or print a suggested citation for the package and write CITATION file if missing
pkg_dir <- normalizePath(".")
pkg_name <- basename(pkg_dir)
message("Attempting to show citation for package located at: ", pkg_dir)
cite <- tryCatch(utils::citation(pkg_name), error = function(e) NULL)
if (!is.null(cite)) {
  print(cite)
  if (!file.exists("CITATION")) {
    message("Writing CITATION file to repository root")
    writeLines(capture.output(print(cite)), con = "CITATION")
  } else {
    message("CITATION already exists; not overwriting.")
  }
} else {
  message("Package not installed locally. To generate a full CITATION, install the package and re-run this script.")
  message("As an alternative, add a CITATION file or use 'usethis::use_citation()' interactively.")
}
