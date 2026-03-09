# tools — Developer and Maintenance Helpers

This directory contains small, focused helper tools to simplify local
development, testing, documentation and release tasks for the TSENAT R
package. Each tool is executable with `Rscript` and is designed to do
one thing well.

## Quick Reference

| Tool | Purpose | Example |
|---|---:|---|
| `run_roxygen.R` | Regenerate `man/` from roxygen comments | `Rscript tools/run_roxygen.R` |
| `run_tests.R` | Run package tests (testthat) | `Rscript tools/run_tests.R` |
| `build_vignettes.R` | Build vignettes (HTML + PDF) → `inst/doc/` | `Rscript tools/build_vignettes.R` |
| `build_pkgdown.R` | Build pkgdown site → `docs/` | `Rscript tools/build_pkgdown.R` |
| `render_vignette.R` | Render a single vignette (HTML + PDF) | `Rscript tools/render_vignette.R` |
| `install_deps.R` | Install package dependencies (remotes::install_deps) | `Rscript tools/install_deps.R` |
| `run_bioccheck.R` | Run BiocCheck (Bioconductor) | `Rscript tools/run_bioccheck.R` |
| `check_style.R` | Check / (optionally) fix code style (styler) | `Rscript tools/check_style.R` |
| `fix_indent_multiple_of_4.R` | Normalize indentation in files | `Rscript tools/fix_indent_multiple_of_4.R <file>` |
| `apply_style.R` | Check or apply code style fixes (styler) | `Rscript tools/apply_style.R` or `Rscript tools/apply_style.R --apply` |

## New maintenance helpers

- `check_package.R` — Run a quick integrity suite (roxygen sync, style,
  tests, dependency checks, doc completeness):
```bash
Rscript tools/check_package.R
```

- `bump_version.R` — Bump package `Version` (DESCRIPTION), update
  `_pkgdown.yml` (if present) and prepend a NEWS entry. Usage:
```bash
Rscript tools/bump_version.R [major|minor|patch]
# default: patch
```

- `prepare_release.R` — Orchestrates pre-release tasks:
  checks → build vignettes → build docs → verify metadata. It prints a

```bash
Rscript tools/prepare_release.R
```

## Usage guidelines (recommended workflows)

### Recently added maintenance helpers

- `run_linters.R` — Format and lint package source using `styler` and `lintr`.
```bash
Rscript tools/run_linters.R
```

- `audit_deps.R` — Audit dependencies and print packages that are out-of-date
  (uses `remotes::package_deps()`).
```bash
Rscript tools/audit_deps.R
```

- `run_rhub_checks.R` — Submit the package to R-hub for platform checks.
  Set `R_HUB_EMAIL` in your environment to receive notifications.
```bash
R_HUB_EMAIL=you@example.org Rscript tools/run_rhub_checks.R
```

- `generate_citation.R` — Print the package citation or write a `CITATION`
  file if none exists (requires the package to be installed to generate a
  full citation automatically).
```bash
Rscript tools/generate_citation.R
```


### Quick development loop

```bash
# regenerate docs after changing roxygen comments
Rscript tools/run_roxygen.R

# run tests
Rscript tools/run_tests.R

# check code style and fix (if desired)
Rscript tools/check_style.R
```
Rscript tools/apply_style.R
# to apply fixes in-place:
Rscript tools/apply_style.R --apply
### Prepare documentation & website

```bash
# build vignettes and pkgdown site locally
Rscript tools/build_vignettes.R
Rscript tools/build_pkgdown.R
Rscript tools/build_vignettes.R
Rscript tools/build_vignettes.R --pdf
```

### Before tagging a release

```bash
Rscript tools/deploy_ghpages.R
# run full package checks
Rscript tools/check_package.R

# bump version (patch, minor or major), update NEWS.md
Rscript tools/bump_version.R patch

Rscript tools/prepare_release.R
```

## CI notes

- CI uses `.circleci/config.yml` to run the canonical build/test/deploy
  steps. These tools are primarily for local development and help
  reproduce CI steps locally.

## Conventions and tips

- Run tools from the package root.
- Tools will attempt to install missing R packages from CRAN/Bioconductor
  when necessary.
- Tools are intended as developer conveniences — they are not a
  replacement for a proper CI pipeline.