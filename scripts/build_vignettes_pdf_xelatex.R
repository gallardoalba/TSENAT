# Render vignettes to PDF using xelatex engine
if (!requireNamespace("rmarkdown", quietly = TRUE)) install.packages("rmarkdown", repos = "https://cloud.r-project.org")
if (!dir.exists("inst")) dir.create("inst")
if (!dir.exists("inst/doc")) dir.create("inst/doc", recursive = TRUE)
vfiles <- list.files("vignettes", pattern = "\\.Rmd$", full.names = TRUE)
if (length(vfiles) == 0) {
  cat("No Rmd vignettes found\n")
  quit(status = 0)
}
for (f in vfiles) {
  cat("Attempting PDF (xelatex) for", f, "-> inst/doc\n")
  if (nzchar(Sys.which("xelatex"))) {
    cat("xelatex found at:", Sys.which("xelatex"), "\n")
  } else {
    cat("Warning: xelatex not found on PATH; attempt may fail.\n")
  }
  pdf_name <- paste0(tools::file_path_sans_ext(basename(f)), ".pdf")
  tryCatch({
    rmarkdown::render(f,
      output_format = rmarkdown::pdf_document(latex_engine = "xelatex"),
      output_file = pdf_name,
      output_dir = "inst/doc",
      clean = TRUE,
      quiet = FALSE
    )
    cat("Rendered PDF:", file.path("inst/doc", pdf_name), "\n")
  }, error = function(e) {
    cat("PDF (xelatex) rendering failed for", f, ":", conditionMessage(e), "\n")
  })
}
cat("Done xelatex attempts\n")
