# Build vignettes into inst/doc/
if (!requireNamespace("rmarkdown", quietly = TRUE)) install.packages("rmarkdown", repos = "https://cloud.r-project.org")
# If available, load the package source so vignette examples use the current
# development code rather than an installed version. This avoids "unused
# argument" errors caused by API drift between the installed package and the
# working tree.
if (requireNamespace("pkgload", quietly = TRUE)) {
    try(pkgload::load_all("."), silent = TRUE)
}
if (!dir.exists("inst")) dir.create("inst")
if (!dir.exists("inst/doc")) dir.create("inst/doc", recursive = TRUE)
vfiles <- list.files("vignettes", pattern = "\\.Rmd$", full.names = TRUE)
if (length(vfiles) == 0) {
    cat("No Rmd vignettes found\n")
    quit(status = 0)
}
for (f in vfiles) {
    cat("Rendering", f, "-> inst/doc\n")
    tryCatch(
        {
            # render HTML (default) into inst/doc
            rmarkdown::render(f, output_dir = "inst/doc", clean = TRUE, quiet = FALSE)
            # attempt to also render PDF; LaTeX must be available on the system
            pdf_name <- paste0(tools::file_path_sans_ext(basename(f)), ".pdf")
            tryCatch(
                {
                    rmarkdown::render(f, output_format = "pdf_document", output_file = pdf_name, output_dir = "inst/doc", clean = TRUE, quiet = FALSE)
                    cat("Rendered PDF:", file.path("inst/doc", pdf_name), "\n")
                },
                error = function(e_pdf) {
                    warning("PDF rendering failed for ", f, ": ", conditionMessage(e_pdf), "\n")
                }
            )
        },
        error = function(e) {
            cat("ERROR rendering", f, ":", conditionMessage(e), "\n")
            quit(status = 1)
        }
    )
}
cat("Vignettes rendered into inst/doc/\n")
