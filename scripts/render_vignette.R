#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
tryCatch({
  pkgload::load_all('.')
  if (!requireNamespace('rmarkdown', quietly = TRUE)) stop('rmarkdown required')
  cat('Rendering HTML...\n')
  rmarkdown::render('vignettes/TSENAT.Rmd', output_dir = 'inst/doc', clean = TRUE, envir = new.env())
  cat('HTML render complete.\n')
  cat('Rendering PDF (XeLaTeX)...\n')
  rmarkdown::render('vignettes/TSENAT.Rmd', output_format = rmarkdown::pdf_document(latex_engine = 'xelatex'), output_dir = 'inst/doc', clean = TRUE, envir = new.env())
  cat('PDF render complete.\n')
}, error = function(e) {
  cat('ERROR:', conditionMessage(e), '\n')
  quit(status = 1)
})
