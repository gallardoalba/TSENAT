#!/usr/bin/env Rscript
# Render a specific vignette to HTML and PDF formats
#
# Usage: Rscript tools/render_vignette.R [vignette_file]
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
        
        args <- commandArgs(trailingOnly = TRUE)

        # Default vignette file; can be overridden by first positional arg
        vignette_file <- if (length(args) >= 1 && !startsWith(args[1], "--")) args[1] else "vignettes/TSENAT.Rmd"
        cat("Rendering", vignette_file, "\n")

        # Enable PDF rendering only if explicit flag provided or env var set
        pdf_flag <- any(args %in% c("--pdf", "-p", "--with-pdf")) ||
            identical(tolower(Sys.getenv("TSENAT_RENDER_PDF", "false")), "true")

        cat("  → HTML...\n")
        rmarkdown::render(
            vignette_file,
            output_dir = "inst/doc",
            clean = TRUE,
            envir = new.env()
        )
        
        if (isTRUE(pdf_flag)) {
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
        } else {
            cat("  → PDF rendering disabled (pass --pdf or set TSENAT_RENDER_PDF=true to enable)\n")
        }
        
        cat("✓ Vignette rendering complete\n")
    },
    error = function(e) {
        cat("✗ ERROR:", conditionMessage(e), "\n")
        quit(status = 1)
    }
)
