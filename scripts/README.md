# scripts — Developer and Maintenance Helpers

This directory contains small, focused helper scripts to simplify local
development, testing, documentation and release tasks for the TSENAT R
package. Each script is executable with `Rscript` and is designed to do
one thing well.

## Quick Reference

| Script | Purpose | Example |
|---|---:|---|
| `run_roxygen.R` | Regenerate `man/` from roxygen comments | `Rscript scripts/run_roxygen.R` |
| `run_tests.R` | Run package tests (testthat) | `Rscript scripts/run_tests.R` |
| `build_vignettes.R` | Build vignettes (HTML + PDF) → `inst/doc/` | `Rscript scripts/build_vignettes.R` |
| `build_pkgdown.R` | Build pkgdown site → `docs/` | `Rscript scripts/build_pkgdown.R` |
| `render_vignette.R` | Render a single vignette (HTML + PDF) | `Rscript scripts/render_vignette.R` |
| `install_deps.R` | Install package dependencies (remotes::install_deps) | `Rscript scripts/install_deps.R` |
| `run_bioccheck.R` | Run BiocCheck (Bioconductor) | `Rscript scripts/run_bioccheck.R` |
| `check_style.R` | Check / (optionally) fix code style (styler) | `Rscript scripts/check_style.R` |
| `fix_indent_multiple_of_4.R` | Normalize indentation in files | `Rscript scripts/fix_indent_multiple_of_4.R <file>` |
| `generate_tcga_rdata.R` | Generate synthetic example dataset (inst/extdata) | `Rscript scripts/generate_tcga_rdata.R` |

## New maintenance helpers

- `check_package.R` — Run a quick integrity suite (roxygen sync, style,
  tests, dependency checks, doc completeness):

```bash
Rscript scripts/check_package.R
```

- `bump_version.R` — Bump package `Version` (DESCRIPTION), update
  `_pkgdown.yml` (if present) and prepend a NEWS entry. Usage:

```bash
Rscript scripts/bump_version.R [major|minor|patch]
# default: patch
```

- `prepare_release.R` — Orchestrates pre-release tasks:
  checks → build vignettes → build docs → verify metadata. It prints a
  short release checklist on success.

```bash
Rscript scripts/prepare_release.R
```

## Usage guidelines (recommended workflows)

### Quick development loop

```bash
# regenerate docs after changing roxygen comments
Rscript scripts/run_roxygen.R

# run tests
Rscript scripts/run_tests.R

# check code style and fix (if desired)
Rscript scripts/check_style.R
```

### Prepare documentation & website

```bash
# build vignettes and pkgdown site locally
Rscript scripts/build_vignettes.R
Rscript scripts/build_pkgdown.R
# preview docs in ./docs/
```

### Before tagging a release

```bash
# run full package checks
Rscript scripts/check_package.R

# bump version (patch, minor or major), update NEWS.md
Rscript scripts/bump_version.R patch

# run the prepare-release workflow
Rscript scripts/prepare_release.R
```

## CI notes

- CI uses `.circleci/config.yml` to run the canonical build/test/deploy
  steps. These scripts are primarily for local development and help
  reproduce CI steps locally.

## Conventions and tips

- Run scripts from the package root.
- Scripts will attempt to install missing R packages from CRAN/Bioconductor
  when necessary.
- Scripts are intended as developer conveniences — they are not a
  replacement for a proper CI pipeline.