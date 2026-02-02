# Check that code is styled. If styling would change files, exit with failure.
# Ensure `styler` is available
if (!requireNamespace("styler", quietly = TRUE)) install.packages("styler", repos = "https://cloud.r-project.org")
# Load styler so internal helpers are registered
library(styler)
# Run styler with 4-space indentation to format in-place
# Only style R source and tests to avoid needing knitr for vignettes
style_fun <- function(...) styler::tidyverse_style(...,
                                                  indent_by = 4,
                                                  include_roxygen_examples = FALSE)
# Style R source files
styler::style_dir(path = "R", style = style_fun)
# Style tests and top-level helper R files
if (dir.exists("tests")) styler::style_dir(path = "tests", style = style_fun)
if (dir.exists("scripts")) styler::style_dir(path = "scripts", style = style_fun)
# If git is available, fail if there are unstaged/uncommitted changes introduced
# by styler
res <- system("git status --porcelain", intern = TRUE)
if (length(res) > 0 && nzchar(paste(res, collapse = "\n"))) {
  cat("Styler reformatted files; please run `Rscript scripts/check_style.R` and commit the changes before merging.\n")
  system("git --no-pager diff")
  quit(status = 1)
}
cat("No style changes\n")
