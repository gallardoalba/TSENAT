#!/usr/bin/env Rscript
## Run formatting and linting for the package
if (!requireNamespace("styler", quietly = TRUE)) install.packages("styler")
if (!requireNamespace("lintr", quietly = TRUE)) install.packages("lintr")
message("Styling package source with styler::style_pkg()")
tryCatch(styler::style_pkg(), error = function(e) message("styler failed: ", e$message))
message("Running lintr::lint_package() — this may produce many messages")
tryCatch(lintr::lint_package(), error = function(e) message("lintr failed: ", e$message))
message("Linters and style complete.")
