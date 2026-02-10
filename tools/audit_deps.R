#!/usr/bin/env Rscript
## Audit package dependencies and report outdated packages
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
message("Checking package dependencies (installed vs available)...")
deps <- tryCatch(remotes::package_deps(dependencies = TRUE), error = function(e) { message("package_deps failed: ", e$message); NULL })
if (is.null(deps)) quit(status = 0)
print(deps)
outdated <- deps[deps$installed < deps$remote, , drop = FALSE]
if (nrow(outdated) > 0) {
  message("Outdated packages detected:")
  print(outdated[, c("package", "installed", "remote")])
} else {
  message("All dependencies appear up-to-date.")
}
