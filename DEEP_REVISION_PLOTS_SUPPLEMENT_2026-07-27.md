# Supplementary Re-Evaluation: TSENAT plots\_\*.R — Additional Findings

**Date**: 2026-07-27  
**Scope**: Second-pass deep analysis of all 10 `plots_*.R` files  
**Status**: Re-Evaluation Complete ✅

> This supplements `DEEP_REVISION_PLOTS_2026-07-27.md` with findings
> from a deeper re-read.

------------------------------------------------------------------------

## Table of Contents

1.  [Critical New Bugs](#id_1--critical-new-bugs)
2.  [Additional Bugs](#id_2--additional-bugs)
3.  [Additional Formatting
    Inconsistencies](#id_3--additional-formatting-inconsistencies)
4.  [Additional Architectural
    Issues](#id_4--additional-architectural-issues)
5.  [Updated Prioritized Action
    Items](#id_5--updated-prioritized-action-items)

------------------------------------------------------------------------

## 1. 🔴 Critical New Bugs

### B1.1 — **Exported function `plot_diversity_violin_density()` defined in TWO files**

**Files**: - `R/plots_transcript.R:1049` - `R/plots_violin_density.R:47`

``` r

# R/plots_transcript.R:1049
#' @export
plot_diversity_violin_density <- function(se, assay_name = "diversity", title = NULL,
    output_file = NULL) { ... }

# R/plots_violin_density.R:47
#' @export
plot_diversity_violin_density <- function(se, assay_name = "diversity", title = NULL,
    output_file = NULL) { ... }
```

**Severity**: **CRITICAL** — This is a silent function overwrite.
Whichever file R sources **last** wins. Both are `@export`, both have
the same signature, both have identical roxygen documentation. The user
has no way to know which version is active.

**Differences between the two versions**:

| Aspect | `plots_transcript.R` version | `plots_violin_density.R` version |
|----|----|----|
| Data prep | `.prepare_long_format()` | `.prepare_long_format()` |
| Save mechanism | Manual `ggsave()` with `.calculate_plot_dims()` | `.save_plot_standard()` |
| `create.dir` in ggsave | `TRUE` | Not set (defaults to `TRUE`) |

**Fix**: Delete one copy. `plots_violin_density.R` is the canonical home
for this function. Remove the duplicate from `plots_transcript.R`.

------------------------------------------------------------------------

### B1.2 — **Internal helpers `.plot_diversity_violin_singleq()` and `.plot_diversity_density_singleq()` duplicated**

**Files**: Same as B1.1

| Function | `plots_transcript.R` | `plots_violin_density.R` |
|----|----|----|
| `.plot_diversity_violin_singleq()` | Line 898 | Line 188 |
| `.plot_diversity_density_singleq()` | Line 961 | Line 128 |

**Severity**: **CRITICAL** — Same silent overwrite issue. The
`plots_violin_density.R` versions include empty
`suppressPackageStartupMessages({})` blocks that the
`plots_transcript.R` versions don’t have. Both versions are functionally
identical otherwise.

**Fix**: Delete the copies from `plots_transcript.R`. These are clearly
intended to live in `plots_violin_density.R` (the file is named for
them).

------------------------------------------------------------------------

## 2. 🟡 Additional Bugs

### B2.1 — `q_values` computed but never used in `.plot_sait()`

**File**: `R/plots_sait.R:111`

``` r

q_values <- .plot_gam_extract_q_values(model_data)
```

The variable `q_values` is assigned but never referenced anywhere else
in the function. If `model_data` is `NULL`,
`.plot_gam_extract_q_values()` presumably returns `NULL` silently, but
the result is discarded.

**Severity**: **LOW** — Dead computation. Likely leftover from a
refactoring where q-value validation was intended but never implemented.

**Fix**: Either use the extracted q-values for validation (compare
against colnames of `mat`) or remove the line.

------------------------------------------------------------------------

### B2.2 — `base_title` computed but never used in `plot_diversity_violin_density()`

**File**: `R/plots_violin_density.R:84` and `R/plots_transcript.R:1086`

``` r

base_title <- title %||% sprintf("Tsallis entropy at q = %g", q_val)
```

`base_title` is computed from the user-supplied `title` parameter (or
auto-generated), but is **never used**. The actual title is hardcoded
two lines later:

``` r

title_grob <- .create_title_grob("Tsallis Entropy Distribution by Group", 
    subtitle = "Violin and density plots across samples", ...)
```

The `title` parameter is documented but silently ignored.

**Severity**: **MEDIUM** — User-facing parameter is ignored. If a user
passes `title = "My Custom Title"`, it has no effect.

**Fix**: Pass `base_title` (or `title`) to `.create_title_grob()`
instead of the hardcoded string.

------------------------------------------------------------------------

### B2.3 — `@param q_value` documented but not a function parameter (3 functions)

**Files**: - `R/plots_violin_density.R:121` —
`.plot_diversity_density_singleq()` - `R/plots_violin_density.R:181` —
`.plot_diversity_violin_singleq()` - `R/plots_transcript.R:953` —
`.plot_diversity_density_singleq()`

All three functions document:

``` r

#' @param q_value The specific q value to plot (numeric, e.g., 1, 2, 0.5).
```

But none have `q_value` in their signatures:

``` r

.plot_diversity_density_singleq <- function(se, assay_name = "diversity", title = NULL) { ... }
.plot_diversity_violin_singleq <- function(se, assay_name = "diversity", title = NULL) { ... }
```

**Severity**: **LOW** — Documentation-only bug, but could confuse
developers. The functions auto-detect q from SE metadata instead of
accepting it as a parameter.

**Fix**: Either add `q_value` as an optional parameter or remove it from
the documentation.

------------------------------------------------------------------------

### B2.4 — Empty `suppressPackageStartupMessages({})` blocks (5 occurrences)

**Locations**:

| File                     | Line | Function                            |
|--------------------------|------|-------------------------------------|
| `plots_spectrum.R`       | 269  | `.spectrum_plot_global()`           |
| `plots_tsallis_q.R`      | 229  | `.plot_tsallis_gene_specific()`     |
| `plots_tsallis_q.R`      | 350  | `.plot_tsallis_bootstrap_ci()`      |
| `plots_violin_density.R` | 129  | `.plot_diversity_density_singleq()` |
| `plots_violin_density.R` | 189  | `.plot_diversity_violin_singleq()`  |

Each is an empty block:

``` r

suppressPackageStartupMessages({
})
```

This is dead code — it suppresses nothing because no library calls are
inside the braces.

**Severity**: **LOW** — No effect, but clutters the code and suggests an
incomplete cleanup.

**Fix**: Remove all 5 empty blocks.

------------------------------------------------------------------------

### B2.5 — `gene_display_name = NULL` parameter never passed non-NULL

**File**: `R/plots_gam_helpers.R:525`

``` r
.plot_gam_make_plot <- function(gene, gene_display_name = NULL, gene_name_map, mat,
    sample_to_group, condition_col) {
```

The only caller is `.plot_sait()` at line 142:

``` r

p <- .plot_gam_make_plot(g, gene_name_map = gene_name_map, mat = mat, 
    sample_to_group = sample_to_group, condition_col = condition_col)
```

`gene_display_name` is never passed, so it’s always `NULL`. The
parameter exists for theoretical flexibility that is never exercised.

**Severity**: **LOW** — Dead parameter. Not harmful but misleading.

**Fix**: Either remove the parameter or leave a comment explaining it’s
reserved for future use.

------------------------------------------------------------------------

## 3. 🟡 Additional Formatting Inconsistencies

### 3.1 — q-extraction block copy-pasted 7 times

The following 18-line pattern appears **7 times** across 3 files:

``` r

# Try to extract q from SE metadata first (best source for single-q SE)
q_val <- NA
if (!is.null(S4Vectors::metadata(se)$q) && length(S4Vectors::metadata(se)$q) > 0) {
    q_vals <- unique(as.numeric(S4Vectors::metadata(se)$q))
    if (length(q_vals) > 0 && !all(is.na(q_vals))) {
        q_val <- q_vals[1]
    }
}

# Fallback: use prepare_tsallis_long for data transformation
long <- .prepare_tsallis_long(se, assay_name = assay_name)

if (nrow(long) == 0)
    stop("No data found in the long format dataframe")

# If still no q, extract from data
if (is.na(q_val)) {
    q_values <- unique(long$q)
    q_values <- q_values[!is.na(q_values)]
    if (length(q_values) > 0) {
        q_val <- q_values[1]
    }
}
```

| \# | File | Line | Function |
|----|----|----|----|
| 1 | `plots_violin_density.R` | 63 | [`plot_diversity_violin_density()`](https://gallardoalba.github.io/TSENAT/reference/plot_diversity_violin_density.md) |
| 2 | `plots_violin_density.R` | 134 | `.plot_diversity_density_singleq()` |
| 3 | `plots_violin_density.R` | 194 | `.plot_diversity_violin_singleq()` |
| 4 | `plots_transcript.R` | 902 | `.plot_diversity_violin_singleq()` (duplicate) |
| 5 | `plots_transcript.R` | 965 | `.plot_diversity_density_singleq()` (duplicate) |
| 6 | `plots_transcript.R` | 1061 | [`plot_diversity_violin_density()`](https://gallardoalba.github.io/TSENAT/reference/plot_diversity_violin_density.md) (duplicate) |
| 7 | `plots_transcript.R` | ~911 | `.plot_diversity_violin_singleq()` (alternate) |

**Severity**: **MEDIUM** — 126+ lines of duplicated code. Any bug fix
must be applied in 7 places.

**Fix**: Extract into a helper function:

``` r

.extract_q_value <- function(se, assay_name = "diversity") {
    q_val <- NA
    if (!is.null(S4Vectors::metadata(se)$q) && length(S4Vectors::metadata(se)$q) > 0) {
        q_vals <- unique(as.numeric(S4Vectors::metadata(se)$q))
        if (length(q_vals) > 0 && !all(is.na(q_vals))) q_val <- q_vals[1]
    }
    long <- .prepare_tsallis_long(se, assay_name = assay_name)
    if (nrow(long) == 0) stop("No data found in the long format dataframe")
    if (is.na(q_val)) {
        q_values <- unique(long$q)
        q_values <- q_values[!is.na(q_values)]
        if (length(q_values) > 0) q_val <- q_values[1]
    }
    list(q_val = q_val, long = long)
}
```

------------------------------------------------------------------------

### 3.2 — `return(NULL)` vs `return(invisible(NULL))` inconsistency

| Pattern | Count | Example location |
|----|----|----|
| `return(NULL)` | 18 | `plots_sait.R:132`, `plots_gam_helpers.R:43` |
| `return(invisible(NULL))` | 8 | `plots_heatmaps.R:205`, `plots_themes.R:1260` |

No consistent rule governs which to use. `return(invisible(NULL))` is
more appropriate for functions called for side effects (saving files,
drawing to device); `return(NULL)` is appropriate for functions that
compute a value. But the codebase doesn’t follow this distinction.

------------------------------------------------------------------------

### 3.3 — `.prepare_tsallis_long()` vs `.prepare_long_format()` inconsistency

Both perform similar data preparation but are used interchangeably:

| Function used | Locations |
|----|----|
| `.prepare_tsallis_long()` | 7 call sites (including the duplicated functions) |
| `.prepare_long_format()` | 5 call sites |

`.prepare_long_format()` is a wrapper around `.prepare_tsallis_long()`
that adds validation and condition_col auto-detection. Some callers
bypass the wrapper, losing validation.

**Severity**: **MEDIUM** — Inconsistent validation coverage.

------------------------------------------------------------------------

### 3.4 — `ggsave()` vs `.save_plot_standard()` inconsistency

| Save mechanism | Locations |
|----|----|
| `.save_plot_standard()` | `plots_violin_density.R`, `plots_tsallis_q.R` (via helper) |
| Manual `ggsave()` | `plots_transcript.R:1111`, `plots_gam_helpers.R:498` |
| [`ggplot2::ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html) via `.plot_gam_save_plot()` | `plots_sait.R:157` |

`.save_plot_standard()` was created specifically to consolidate the
ggsave pattern, yet two locations still use raw `ggsave()`.

------------------------------------------------------------------------

## 4. 🟢 Additional Architectural Issues

### 4.1 — Function ownership confusion: plots_transcript.R vs plots_violin_density.R

The `plot_diversity_violin_density` function and its two helpers
(`.plot_diversity_violin_singleq`, `.plot_diversity_density_singleq`)
exist in **both** files. The `plots_transcript.R` copies appear to be an
accidental paste during the July 2026 refactoring when code was
extracted from `plots_helpers.R`.

**Evidence**: The `plots_transcript.R` header says “Extracted from
plots_helpers.R — July 2026 refactoring (I11)” but the violin/density
functions were split into their own file `plots_violin_density.R`. The
old copies in `plots_transcript.R` were never removed.

------------------------------------------------------------------------

### 4.2 — Missing validation asymmetry in `.plot_sait()`

**File**: `R/plots_sait.R`

`.plot_sait()` validates: - `se` is a SummarizedExperiment ✓ -
`sait_res` is a data.frame ✓ - `condition_col` exists in colData ✓ -
`assay_name` exists (implicitly via `assay()` call) ✓

But does NOT validate: - `n_top` is positive integer ✗ - `sig_alpha` is
between 0 and 1 ✗ - `assay_name` actually exists before calling
`assay()` ✗ - `output_file` is a writable path ✗

Compare with `.plot_divergence_spectrum()` which validates `n_genes` and
`ncol` explicitly.

------------------------------------------------------------------------

### 4.3 — `.prepare_transcript_inputs()` global option side-effect removed from new version but old version keeps it

**File**: `R/plots_transcript.R`

The old function `.make_plot_for_geneprepare_inputs()` (line 30) has:

``` r

.cnt <- as.integer(getOption("TSENAT.plot_top_counter", 0)) + 1L
options(TSENAT.plot_top_counter = .cnt)
```

The new function `.prepare_transcript_inputs()` (line 230) does NOT have
this side-effect. This means different behavior depending on which
function is called.

------------------------------------------------------------------------

### 4.4 — `sprintf("%g", q_val)` format inconsistency

The q-extraction block uses
`sprintf("Tsallis entropy at q = %g", q_val)`. `%g` uses “general”
format which switches between `%f` and `%e` depending on magnitude. For
`q = 1`, this produces `"q = 1"` (no decimal). For `q = 0.01`, this
produces `"q = 0.01"`.

Other parts of the codebase use `sprintf("q = %.2f", q_val)` (in
`.format_q_label()` at `plots_validation.R:150`). This inconsistency
means q-value labels can appear as `"q = 1"` in one plot and
`"q = 1.00"` in another.

------------------------------------------------------------------------

## 5. 📋 Updated Prioritized Action Items

These supplement the action items from the original report.

| \# | Action | Impact | Effort | Category |
|----|----|----|----|----|
| **N1** | Delete duplicate `plot_diversity_violin_density` from `plots_transcript.R` | **CRITICAL** | 5 min | Bug |
| **N2** | Delete duplicate `.plot_diversity_violin_singleq` from `plots_transcript.R` | **CRITICAL** | 2 min | Bug |
| **N3** | Delete duplicate `.plot_diversity_density_singleq` from `plots_transcript.R` | **CRITICAL** | 2 min | Bug |
| **N4** | Remove all 5 empty `suppressPackageStartupMessages({})` blocks | Low | 5 min | Dead code |
| **N5** | Extract `q_val` detection into `.extract_q_value()` helper (7→1) | Medium | 30 min | DRY |
| **N6** | Fix `base_title` unused — pass to `.create_title_grob()` | Medium | 5 min | Bug |
| **N7** | Fix `@param q_value` documentation mismatch (3 functions) | Low | 5 min | Docs |
| **N8** | Remove or comment `gene_display_name` dead parameter | Low | 2 min | Dead code |
| **N9** | Remove or use `q_values` in `.plot_sait()` | Low | 2 min | Dead code |
| **N10** | Standardize `return(NULL)` vs `return(invisible(NULL))` | Low | 30 min | Consistency |
| **N11** | Standardize `.prepare_tsallis_long()` vs `.prepare_long_format()` usage | Medium | 1 hr | Consistency |
| **N12** | Standardize save mechanism (use `.save_plot_standard()` everywhere) | Low | 20 min | Consistency |
| **N13** | Add missing parameter validation to `.plot_sait()` | Low | 15 min | Robustness |
| **N14** | Standardize q-value label format (`%g` vs `%.2f`) | Low | 10 min | Consistency |

------------------------------------------------------------------------

## Combined Severity Summary

| Severity | Count | Key Items |
|----|----|----|
| **CRITICAL** | 3 | Duplicate exported function, duplicate internal helpers |
| **MEDIUM** | 4 | Unused `base_title` (ignores user input), 7x copy-pasted q-extraction, inconsistent data prep wrapper usage, missing validation |
| **LOW** | 13 | Dead code, doc mismatches, dead params, inconsistencies |

------------------------------------------------------------------------

*Supplementary report generated 2026-07-27.*
