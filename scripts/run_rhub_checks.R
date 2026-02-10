#!/usr/bin/env Rscript
## Submit checks to R-hub for common platforms. Requires R-hub account/email.
if (!requireNamespace("rhub", quietly = TRUE)) install.packages("rhub")
email <- Sys.getenv("R_HUB_EMAIL")
if (identical(email, "")) message("Warning: R_HUB_EMAIL not set; rhub may email results to package maintainer if not provided.")
platforms <- c("debian-clang-devel", "debian-gcc-release", "windows-x86_64-devel")
for (p in platforms) {
  message("Submitting check to R-hub: ", p)
  res <- tryCatch(rhub::check(platform = p, email = if (nzchar(email)) email else NULL), error = function(e) e)
  print(res)
}
