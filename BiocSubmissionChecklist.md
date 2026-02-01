Bioconductor submission checklist

1. DESCRIPTION
- Replace `Author:` and `Maintainer:` with a single `Authors@R:` field.
  - Include maintainer with role 'cre' (and optionally ORCID in comment).
  - Consider adding 'fnd' role if funded.
- Add appropriate `biocViews:` entries (suggested: Software, GeneExpression, DifferentialExpression).
- Ensure `License:` is a valid Bioconductor license (e.g., `GPL-3`).
- Ensure `VignetteBuilder: knitr` and `VignetteBuilder` package present in Suggests if vignettes included.

2. VIGNETTES
- Keep vignette source (`vignettes/*.Rmd`) in repo.
- For Bioc submission, include built vignettes in `inst/doc/` (HTML and PDF) or accept the NOTE; to include them:
  - Build vignettes locally: `R CMD build --no-build-vignettes .` then `R CMD check` with vignettes built, or use `devtools::build_manual()`.
  - Place generated HTML/PDF into `inst/doc/` before creating the release tarball.

3. NAMESPACE & MAN
- Ensure `roxygen2` docs regenerate cleanly (`roxygen2::roxygenise('.')`).
- Fix any man page long lines (>80 chars) and run `R CMD check` locally.

4. TESTS
- Ensure `tests/testthat` covers main functionality and that tests use `skip_on_bioc()` where needed.

5. CODING STYLE
- Run `styler::style_dir('R')` to conform to Bioconductor style.
- Replace `sapply()` with `vapply()` or `lapply()` where the return type is known.
- Replace `1:` patterns with `seq_len()` or `seq_along()`.

6. CI
- Use GitHub Actions provided by r-lib or Bioconductor template and ensure `remotes::install_deps(dependencies = c('Depends','Imports','LinkingTo'))` is used.
- Optionally install `BiocManager::install('BiocCheck')` and run `BiocCheck::BiocCheck('.')` in CI for pre-submission checks.

7. SUPPORT SITE
- Register a support site at: https://bioconductor.org/support/ (you must add the maintainer email/address there).
- After registration, ensure the email appears on the support site so `BiocCheck` can validate it.

8. SUBMISSION
- Tag a release (e.g., `git tag -a vX.Y.Z -m 'Release X.Y.Z'`) and push tags.
- Follow Bioconductor submission instructions: https://bioconductor.org/developers/package-submission/

9. FINAL CHECKS
- Run `devtools::check()` or `R CMD check --as-cran` locally and resolve ERRORS/WARNINGS/NOTEs.
- Run `BiocCheck::BiocCheck('.')` and resolve ERRORS (support-site and Authors@R are common blockers).

Helpful commands:

```sh
# style code
Rscript -e "if (!requireNamespace('styler', quietly=TRUE)) install.packages('styler'); styler::style_dir('R')"
# regenerate docs
Rscript -e "if (!requireNamespace('roxygen2', quietly=TRUE)) install.packages('roxygen2'); roxygen2::roxygenise('.')"
# run BiocCheck
Rscript -e "if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager'); BiocManager::install('BiocCheck'); BiocCheck::BiocCheck('.')"
```

If you want, I can prepare a suggested `Authors@R` snippet based on current `DESCRIPTION` contents and apply it for you.