#!/usr/bin/env Rscript
## Build pkgdown and deploy to gh-pages by invoking the shell deploy script.
## Usage: Rscript tools/deploy_ghpages.R

tryCatch({
  cat("Building pkgdown site (via tools/build_pkgdown.R)\n")
  res_build <- system2("Rscript", c("tools/build_pkgdown.R"))
  if (res_build != 0) stop("build_pkgdown.R failed with exit code: ", res_build)

  cat("Calling deploy shell script (tools/deploy_ghpages.sh)\n")
  if (!file.exists("tools/deploy_ghpages.sh")) stop("deploy_ghpages.sh not found")
  res_deploy <- system2("bash", c("tools/deploy_ghpages.sh"))
  if (res_deploy != 0) stop("deploy_ghpages.sh failed with exit code: ", res_deploy)

  cat("Done.\n")

}, error = function(e) {
  message("ERROR: ", conditionMessage(e))
  quit(status = 1)
})
