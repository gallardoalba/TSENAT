#!/usr/bin/env Rscript
# Render a specific vignette to HTML and PDF formats
#
# Usage: Rscript scripts/render_vignette.R [vignette_file]
#        (defaults to vignettes/TSENAT.Rmd)

tryCatch(
    {
        if (!requireNamespace("pkgload", quietly = TRUE)) {
            install.packages("pkgload", repos = "https://cloud.r-project.org")
        }
        if (!requireNamespace("rmarkdown", quietly = TRUE)) {
            install.packages("rmarkdown", repos = "https://cloud.r-project.org")
        }
        
        cat("Loading package...\n")
        pkgload::load_all(".")
        
        vignette_file <- "vignettes/TSENAT.Rmd"
        cat("Rendering", vignette_file, "\n")
        
        cat("  → HTML...\n")
        rmarkdown::render(
            vignette_file,
            output_dir = "inst/doc",
            clean = TRUE,
            envir = new.env()
        )
        
        cat("  → PDF (XeLaTeX)...\n")
        tryCatch(
            {
                rmarkdown::render(
                    vignette_file,
                    output_format = rmarkdown::pdf_document(latex_engine = "xelatex"),
                    output_dir = "inst/doc",
                    clean = TRUE,
                    envir = new.env()
                )
            },
            error = function(e) {
                cat("  ⚠ PDF rendering skipped:", conditionMessage(e), "\n")
            }
        )
        
        cat("✓ Vignette rendering complete\n")
    },
    error = function(e) {
        cat("✗ ERROR:", conditionMessage(e), "\n")
        quit(status = 1)
    }
)
