# Check that code is styled. If styling would change files, exit with failure.
if (!requireNamespace("styler", quietly = TRUE)) install.packages("styler", repos = "https://cloud.r-project.org")
# Run styler to format in-place (allow roxygen examples; requires roxygen2 to be
# installed)
styler::style_pkg()
# If git is available, fail if there are unstaged/uncommitted changes introduced
# by styler
res <- system("git status --porcelain", intern = TRUE)
if (length(res) > 0 && nzchar(paste(res, collapse = "\n"))) {
  cat("Styler reformatted files; please run `Rscript scripts/check_style.R` and commit the changes before merging.\n")
  system("git --no-pager diff")
  quit(status = 1)
}
cat("No style changes\n")
