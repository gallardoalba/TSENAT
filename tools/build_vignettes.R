#!/usr/bin/env Rscript
# Build vignettes (HTML and PDF) into inst/doc/
#
# Usage: Rscript tools/build_vignettes.R [--pdf]
# If `--pdf` is provided, a PDF version of each vignette will also be
# rendered (requires a working LaTeX installation). Without `--pdf` only
# HTML vignettes are built.

if (!requireNamespace("rmarkdown", quietly = TRUE)) {
    install.packages("rmarkdown", repos = "https://cloud.r-project.org")
}

# Load package source so vignette examples use current development code
if (requireNamespace("pkgload", quietly = TRUE)) {
    try(pkgload::load_all("."), silent = TRUE)
}

# Create output directory
if (!dir.exists("inst")) dir.create("inst")
if (!dir.exists("inst/doc")) dir.create("inst/doc", recursive = TRUE)

# Find all Rmd vignettes
args <- commandArgs(trailingOnly = TRUE)
pdf_flag <- any(args %in% c("--pdf", "-p"))
if (any(args %in% c("--help", "-h"))) {
    cat("Usage: Rscript tools/build_vignettes.R [--pdf]\n")
    cat("  --pdf / -p: also render PDF versions (requires LaTeX)\n")
    quit(status = 0)
}

vfiles <- list.files("vignettes", pattern = "\\.Rmd$", full.names = TRUE)

if (length(vfiles) == 0) {
    cat("✗ No Rmd vignettes found in vignettes/\n")
    quit(status = 0)
}

cat("Building", length(vfiles), "vignette(s)...\n")

for (f in vfiles) {
    fname <- basename(f)
    cat("\n  Rendering", fname, "\n")
    
    tryCatch(
        {
            # Render HTML (default) into inst/doc
            cat("    → HTML...\n")
            rmarkdown::render(
                f,
                output_dir = "inst/doc",
                clean = TRUE,
                quiet = FALSE
            )
            
            # Optionally render PDF only when requested
            if (pdf_flag) {
                pdf_name <- paste0(tools::file_path_sans_ext(basename(f)), ".pdf")
                tryCatch(
                    {
                        cat("    → PDF (XeLaTeX)...\n")
                        rmarkdown::render(
                            f,
                            output_format = "pdf_document",
                            output_file = pdf_name,
                            output_dir = "inst/doc",
                            clean = TRUE,
                            quiet = FALSE
                        )
                        cat("    ✓ Rendered:", file.path("inst/doc", pdf_name), "\n")
                    },
                    error = function(e_pdf) {
                        cat("    ⚠ PDF rendering skipped:", conditionMessage(e_pdf), "\n")
                    }
                )
            } else {
                cat("    → PDF skipped (use --pdf to enable)\n")
            }
        },
        error = function(e) {
            cat("  ✗ ERROR rendering", fname, ":", conditionMessage(e), "\n")
            quit(status = 1)
        }
    )
}

cat("\n✓ Vignettes rendered successfully into inst/doc/\n")
