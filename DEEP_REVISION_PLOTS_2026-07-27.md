# Deep Revision: TSENAT plots\_\*.R — Bugs, Inconsistencies & Unification Plan

**Date**: 2026-07-27  
**Scope**: 10 files, 7,018 total lines  
**Status**: Analysis Complete ✅

------------------------------------------------------------------------

## Table of Contents

1.  [File Inventory](#id_1-file-inventory)
2.  [Bugs](#id_2--bugs)
3.  [Formatting Inconsistencies](#id_3--formatting-inconsistencies)
4.  [Architectural Observations](#id_4--architectural-observations)
5.  [Unification Proposal](#id_5--unification-proposal)
6.  [Prioritized Action Items](#id_6--prioritized-action-items)
7.  [Summary](#id_7-summary)

------------------------------------------------------------------------

## 1. File Inventory

| \# | File | Lines | Purpose |
|----|----|----|----|
| 1 | `plots_themes.R` | 1,382 | Theme engine, palette dispatch, legend/axis helpers, grid assembly, plot saving |
| 2 | `plots_transcript.R` | 1,695 | Transcript-level expression heatmaps, divergence distribution, grid composition |
| 3 | `plots_heatmaps.R` | 1,304 | Multi-q delta influence heatmaps, per-gene expression heatmaps, layout helpers |
| 4 | `plots_tsallis_q.R` | 581 | Tsallis entropy q-curve (aggregate, gene-specific, bootstrap CI) |
| 5 | `plots_gam_helpers.R` | 578 | GAM fitting, gene selection, grid arrangement for SAIT plots |
| 6 | `plots_stats.R` | 532 | Diversity spectrum stats, gene filtering/ranking, distribution stats, CI extraction |
| 7 | `plots_spectrum.R` | 397 | Divergence spectrum (global, top genes, single gene, bootstrap CI) |
| 8 | `plots_violin_density.R` | 233 | Violin + density combo grid, single-q density, single-q violin |
| 9 | `plots_sait.R` | 160 | `.plot_sait()` — main GAM q-curve orchestration |
| 10 | `plots_validation.R` | 156 | SE validation, p-value formatting, q-value labeling |

### Dependency Graph

    plots_validation.R  ──┐
    plots_stats.R       ──┤
    plots_gam_helpers.R ──┼──► plots_sait.R
    plots_spectrum.R    ──┤
    plots_tsallis_q.R   ──┤
    plots_violin_density.R ─┤
    plots_heatmaps.R    ──┤
    plots_transcript.R  ──┘
            │
            ▼
    plots_themes.R  (central theme/palette/legend engine)
            │
            ▼
    globals.R  (.font_sizes, .theme_base, .theme_spectrum)

------------------------------------------------------------------------

## 2. 🔴 Bugs

### B2.1 — Truncated file: `plots_validation.R` ends with dangling roxygen

**File**: `R/plots_validation.R`  
**Lines**: 152–156

``` r

#' Format Label for Display
#'
#' @noRd
```

No function body follows the roxygen block. The file ends immediately
after. Meanwhile, `.format_label()` is actually defined at
`plots_transcript.R:1566`.

**Severity**: **LOW** — No runtime error, but confusing and could cause
`R CMD check` documentation warnings.

**Root cause**: The function was likely moved to `plots_transcript.R`
during the July 2026 refactoring but the roxygen block was not removed
from the original location.

**Fix**: Delete lines 152–156 from `plots_validation.R`.

------------------------------------------------------------------------

### B2.2 — Duplicated `has_valid_ci` computation in `.create_ci_ribbon_plot()`

**File**: `R/plots_themes.R`  
**Lines**: ~586–616

The `has_valid_ci` logical is computed **twice** within the same
function:

``` r

# FIRST computation (line ~586):
has_valid_ci <- FALSE
if (ci_lower_col %in% colnames(data) && ci_upper_col %in% colnames(data)) {
    has_valid_ci <- any(!is.na(data[[ci_lower_col]]) & !is.infinite(data[[ci_lower_col]]) &
        !is.na(data[[ci_upper_col]]) & !is.infinite(data[[ci_upper_col]]))
}

# ... 20+ lines of ggplot construction ...

# SECOND computation (line ~610):
has_valid_ci <- FALSE
if (ci_lower_col %in% colnames(data) && ci_upper_col %in% colnames(data)) {
    has_valid_ci <- any(!is.na(data[[ci_lower_col]]) & !is.infinite(data[[ci_lower_col]]) &
        !is.na(data[[ci_upper_col]]) & !is.infinite(data[[ci_upper_col]]))
}
```

The second computation shadows the first with an identical value. The
first computation also filters `data_ci` which is never actually used
(the ribbon uses the unfiltered `data`).

**Severity**: **LOW** — No incorrect behavior, just wasteful
recomputation.

**Fix**: Remove the second `has_valid_ci` block (lines ~610–616); reuse
the first result.

------------------------------------------------------------------------

### B2.3 — `%||%` operator dependency ordering risk

**File**: `R/orchestration.R:441` (definition)  
**Used in**: 4 plot files

The `%||%` infix operator is defined in `R/orchestration.R`:

``` r

# R/orchestration.R:441
`%||%` <- function(x, y) {
    if (is.null(x) || (length(x) == 1 && is.na(x))) y else x
}
```

But it is used in:

| File                     | Line | Context                                    |
|--------------------------|------|--------------------------------------------|
| `plots_heatmaps.R`       | 831  | `gene_info_list[[i]]$n_transcripts %||% 1` |
| `plots_themes.R`         | 301  | `color = border_color %||% "black"`        |
| `plots_transcript.R`     | 926  | `title_use <- title %||% sprintf(...)`     |
| `plots_violin_density.R` | 84   | `base_title <- title %||% sprintf(...)`    |

If R sources any of these plot files before `orchestration.R`, runtime
errors occur with `could not find function "%||%"`. In practice,
`orchestration.R` loads first due to alphabetical ordering, but this is
fragile.

**Severity**: **LOW-MEDIUM** — Works under current load order but
fragile.

**Fix**: Move `%||%` definition to `globals.R` (earliest-loaded utility
file).

------------------------------------------------------------------------

### B2.4 — Confusing sentinel value: `cellwidth=0` means “adaptive”

**File**: `R/plots_heatmaps.R`  
**Lines**: 113, 349, 928

The function signatures declare:

``` r
.plot_jis_delta <- function(..., cellwidth = 0, cellheight = 0, ...)
.plot_expression <- function(..., cellwidth = 0, cellheight = 0, ...)
.calculate_adaptive_cellsizes <- function(..., cellwidth = 0, cellheight = 0, ...)
```

A value of `0` is treated as “auto/adaptive” by internal logic. Using
`0` as a sentinel is misleading — it implies “no width” rather than
“compute automatically.”

**Severity**: **LOW** — Works correctly but API is unintuitive.

**Fix**: Change the sentinel to `NULL` or `NA` and update
`.calculate_adaptive_cellsizes()` to check `is.null(cellwidth)`.

------------------------------------------------------------------------

### B2.5 — Dead debug block in `.plot_expression()`

**File**: `R/plots_heatmaps.R`  
**Lines**: 393–396

``` r

if (FALSE) {
    # Debug mode - set to TRUE if needed
    message("[DEBUG] After resolution, genes: ", paste(gene, collapse = ", "))
    message("[DEBUG] tx2gene$Gen unique values (first 10): ", ...)
}
```

The `if (FALSE)` guard means this block will never execute. It’s dead
code.

**Severity**: **COSMETIC** — No runtime impact but clutters the file.

**Fix**: Remove the dead block or use a proper debug flag (e.g.,
`getOption("TSENAT.debug", FALSE)`).

------------------------------------------------------------------------

## 3. 🟡 Formatting Inconsistencies

### 3.1 — Font size chaos (the single biggest problem)

#### Global Constants (from `globals.R:286`)

``` r

.font_sizes <- list(
    title        = 21,
    subtitle     = 17,
    axis_title   = 17,
    axis_text    = 14,
    legend_title = 14,
    legend_text  = 13,
    heatmap_main = 14,
    heatmap_labels = 12
)
```

#### Actual Usage Across Files

| Context | `title_size` | `subtitle_size` | `base_size` | Source |
|----|----|----|----|----|
| `.font_sizes` global default | **21** | **17** | — | `globals.R:286` |
| `.apply_publication_theme` (via globals) | 21 | 17 | 11 | `plots_themes.R:195` |
| Gene-specific q-curve | **16** ✗ | — | 11 | `plots_tsallis_q.R:318` |
| Bootstrap gene CI q-curve | **14** ✗ | — | 11 | `plots_tsallis_q.R:451` |
| Basic gene q-curve fallback | **14** ✗ | — | 11 | `plots_tsallis_q.R:518` |
| GAM grid title | **20** | **16** | — | `plots_gam_helpers.R:447` |
| Violin+density grid | **19** | **15** | — | `plots_violin_density.R:99` |
| Divergence distribution | **19** | **15** | — | `plots_transcript.R:1101` |
| Cowplot gene grid | **18** | **14** | 11 | `plots_transcript.R:1482` |
| Patchwork gene grid | 21 (globals) | 17 (globals) | 11 | `plots_transcript.R:1462` |
| `.plot_jis_delta()` heatmap | 18 (fontsize) | — | — | `plots_heatmaps.R:113` |
| `.plot_expression()` heatmap | 16 (fontsize) | — | — | `plots_heatmaps.R:349` |
| `.create_tsenat_heatmap()` wrapper | 11 (row/col) | — | — | `plots_themes.R:88` |

**Summary**: - **Title sizes span 14 to 21** — a 50% range. The global
default of 21 is almost never used as-is. - **Subtitle sizes span 14 to
17** — a 21% range. - **Heatmap font sizes span 11 to 18** — a 64%
range. - Each plot type independently chose its own sizes with no
apparent system.

------------------------------------------------------------------------

### 3.2 — `base_theme` selection inconsistency

| File | Theme used | Legend position |
|----|----|----|
| `plots_tsallis_q.R` | `theme_spectrum` | right (from theme) |
| `plots_spectrum.R` | `theme_base` | bottom (from `.apply_group_aesthetics`) |
| `plots_transcript.R` | `.theme_base()` direct | bottom |
| `plots_gam_helpers.R` | `.theme_spectrum()` direct | none (per-plot) |
| `plots_violin_density.R` | `theme_base` | bottom |

**Problem**: Tsallis entropy q-curves and divergence q-spectra are
conceptually the same type of plot (both show trends across q-values),
but one uses `theme_spectrum` and the other `theme_base`. There is no
documented rationale for when to use which.

**The only differences between the two themes**:

``` r

.theme_spectrum <- function(base_size = 12) {
    .theme_base(base_size = base_size) + ggplot2::theme(
        legend.position = "right",
        panel.grid.major.y = ggplot2::element_line(color = "gray90", linewidth = 0.25),
        plot.margin = ggplot2::margin(t = 5, r = 8, b = 5, l = 5, unit = "mm")
    )
}
```

------------------------------------------------------------------------

### 3.3 — Y-axis label notation inconsistency

| File | Mode | Y-axis label |
|----|----|----|
| `plots_tsallis_q.R` | Aggregate/basic | `expression("Tsallis entropy (" * S[q] * ")")` — LaTeX-style |
| `plots_tsallis_q.R` | Gene-specific | `"Tsallis entropy"` — plain text |
| `plots_tsallis_q.R` | Bootstrap CI | `"Tsallis entropy"` — plain text |
| `plots_spectrum.R` | Global | `expression("Divergence D[q]")` — LaTeX-style |
| `plots_spectrum.R` | Single gene | `"Divergence D_q"` — plain text |
| `plots_violin_density.R` | All modes | `"Tsallis entropy"` — plain text |

**Problem**: The same plot concept (Tsallis entropy or divergence) uses
LaTeX-like math notation in aggregate mode but plain ASCII in
gene-specific mode. Users see different label styles for the same
quantity.

------------------------------------------------------------------------

### 3.4 — Title boldness double-application

**Files**: `plots_spectrum.R`, `plots_transcript.R`

Some plot functions use `expression(bold(...))` in the title:

``` r

# plots_spectrum.R:309
title = expression(bold("Global Divergence Spectrum: Average " * D[q]))

# plots_transcript.R:1673
title = expression(bold("Distribution of Tsallis Divergence (" ~ D[q] ~ ") effect sizes across genes"))
```

But `.apply_publication_theme()` already applies `face = "bold"` to all
titles:

``` r

# plots_themes.R:207
result <- result + ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5, size = title_size, face = "bold"),
    ...
)
```

**Effect**: The title gets `expression(bold(...))` (math-mode bold) AND
`element_text(face = "bold")` (theme bold). ggplot2 renders correctly
(last wins), but the code is redundant and confusing.

------------------------------------------------------------------------

### 3.5 — Legend handling: three parallel implementations

Three different approaches to legend extraction and grid assembly
coexist:

| Approach | Location | Mechanism |
|----|----|----|
| Centralized helper | `plots_themes.R` — `.assemble_grid_plot()` | `extract_legend` parameter, uses `.create_title_grob()` |
| Manual cowplot | `plots_transcript.R` — `.make_plot_for_genecombine_cowplot()` | Manual [`cowplot::get_legend()`](https://wilkelab.org/cowplot/reference/get_legend.html) + `plot_grid()` |
| Manual patchwork | `plots_transcript.R` — `.make_plot_for_genecombine_patchwork()` | [`patchwork::plot_annotation()`](https://patchwork.data-imaginist.com/reference/plot_annotation.html) |
| Manual base grid | `plots_transcript.R` — `.make_plot_for_genecombine_grid()` | [`grid::grob`](https://rdrr.io/r/grid/grid.grob.html) extraction |
| Manual GAM | `plots_gam_helpers.R` — `.plot_gam_arrange_grid()` | Manual [`cowplot::get_legend()`](https://wilkelab.org/cowplot/reference/get_legend.html) + custom spacing |

The centralized `.assemble_grid_plot()` helper exists but is only used
by `plots_tsallis_q.R`. The transcript and GAM plot code uses four
separate manual implementations.

------------------------------------------------------------------------

### 3.6 — Heatmap default parameter divergence

| Function                   | `fontsize` default  | `cellwidth` default       |
|----------------------------|---------------------|---------------------------|
| `.plot_jis_delta()`        | 18                  | 0 (adaptive)              |
| `.plot_expression()`       | 16                  | 0 (adaptive)              |
| `.create_pheatmap_grob()`  | 18                  | 0 (adaptive)              |
| `.create_tsenat_heatmap()` | 11 (row/col labels) | — (delegates to pheatmap) |

------------------------------------------------------------------------

## 4. 🟢 Architectural Observations

### 4.1 — Concatenated function names (botched refactoring artifact)

**File**: `R/plots_transcript.R`

Several functions have names that appear to be the result of
concatenating a namespace prefix with a function name during a
refactoring that went wrong:

``` r

.make_plot_for_geneprepare_inputs       # line 30   — should be .prepare_transcript_inputs()
.make_plot_for_genemake_plot_for_gene   # line 143  — should be .make_single_gene_plot()
.make_plot_for_genecombine_plots        # line 167  — should be .combine_gene_plots()
.make_plot_for_genebuild_tx_long        # line 1349 — should be .build_transcript_long()
.make_plot_for_geneaggregate_df_long    # line 1364 — should be .aggregate_transcript_data()
.make_plot_for_genebuild_plot_from_summary  # line 1371 — should be .build_plot_from_summary()
.make_plot_for_genecombine_patchwork    # line 1412
.make_plot_for_genecombine_cowplot      # line 1467
.make_plot_for_genecombine_grid         # line 1493
.make_plot_for_generead_tx2gene         # line 1313
.make_plot_for_genemake_agg             # line 1328
```

These names suggest someone attempted to namespace functions under
`make_plot_for_gene` but the prefix concatenation was applied literally
instead of using dot-notation sub-functions.

------------------------------------------------------------------------

### 4.2 — Duplicate logic: old vs new function pairs

**File**: `R/plots_transcript.R`

| Old function (concatenated name) | New function (proper name) | Lines |
|----|----|----|
| `.make_plot_for_geneprepare_inputs()` | `.prepare_transcript_inputs()` | ~109 vs ~114 |
| `.make_plot_for_genebuild_tx_long()` | `.build_transcript_long()` | ~15 vs ~23 |
| `.make_plot_for_geneaggregate_df_long()` | `.aggregate_transcript_data()` | ~7 vs ~8 |
| `.make_plot_for_generead_tx2gene()` | `.read_tx2gene()` | ~17 vs ~23 |
| `.make_plot_for_genemake_agg()` | `.create_aggregation_function()` | ~16 vs ~12 |

The old and new versions do nearly identical things. The old version has
a global option side-effect (`TSENAT.plot_top_counter`) that the new
version does not. Both remain in the codebase.

------------------------------------------------------------------------

### 4.3 — `.theme_spectrum` is a thin wrapper

**File**: `R/globals.R:314`

``` r

.theme_spectrum <- function(base_size = 12) {
    .theme_base(base_size = base_size) + ggplot2::theme(
        legend.position = "right",
        panel.grid.major.y = ggplot2::element_line(color = "gray90", linewidth = 0.25),
        plot.margin = ggplot2::margin(t = 5, r = 8, b = 5, l = 5, unit = "mm")
    )
}
```

The only differences from `.theme_base` are: 1. Legend on the right (vs
bottom) 2. Horizontal grid lines 3. Custom plot margins

This is a thin wrapper that could be a parameter on `.theme_base`.

------------------------------------------------------------------------

### 4.4 — `.plot_divergence_distribution` uses inconsistent styling

**File**: `R/plots_transcript.R:1637`

Unlike all other plot functions, `.plot_divergence_distribution()` does
NOT use `.apply_publication_theme()`. Instead, it manually applies:

``` r

.theme_base(base_size = 12) + 
ggplot2::theme(
    panel.grid.major = ggplot2::element_line(color = "gray90"),
    plot.title = ggplot2::element_text(size = 15, hjust = 0.5),
    plot.subtitle = ggplot2::element_text(size = 12, hjust = 0.5)
)
```

The hardcoded `size = 15` for title and `size = 12` for subtitle are
different from any other plot.

------------------------------------------------------------------------

### 4.5 — Redundant `.apply_tsenat_theme()` and `.set_plot_title()`

**File**: `R/plots_themes.R:23, 47`

Two functions exist that overlap with `.apply_publication_theme()`:

- `.apply_tsenat_theme()` — Applies theme + title centering (superseded
  by `.apply_publication_theme`)
- `.set_plot_title()` — Updates title/subtitle on existing plot
  (superseded by `.apply_publication_theme`’s built-in title handling)

Neither is called anywhere in the codebase. They are dead code.

------------------------------------------------------------------------

## 5. 🎯 Unification Proposal

### 5.1 — Font Size Standardization

Define a **single font hierarchy** in `globals.R` and enforce it via
`.apply_publication_theme()` everywhere:

| Token | Current | Proposed | Rationale |
|----|----|----|----|
| `title` | 21 | **18** | 21 is too large for standard layouts; 18 is the most-used override |
| `subtitle` | 17 | **14** | 17 is too large; 14 is the most-used override |
| `axis_title` | 17 | **14** | Reduce for cleaner look |
| `axis_text` | 14 | **12** | Slightly smaller |
| `legend_title` | 14 | **13** | Minor reduction |
| `legend_text` | 13 | **11** | Minor reduction |
| `heatmap_main` | 14 | **14** | Keep |
| `heatmap_labels` | 12 | **11** | Slightly smaller |

Then **remove ALL hardcoded `title_size=`, `subtitle_size=` overrides**
in plot files. Every plot should call
`.apply_publication_theme(p, title = "...")` and let the global
constants control sizing.

**Files affected**: `plots_tsallis_q.R`, `plots_gam_helpers.R`,
`plots_violin_density.R`, `plots_transcript.R`, `plots_spectrum.R`.

------------------------------------------------------------------------

### 5.2 — Theme Unification

Merge `.theme_spectrum` into `.theme_base` with a `type` parameter:

``` r

.theme_base <- function(base_size = 12, type = c("standard", "spectrum")) {
    type <- match.arg(type)
    t <- ggplot2::theme_minimal(base_size = base_size) +
        ggplot2::theme(
            plot.title    = ggplot2::element_text(hjust = 0.5, face = "bold", 
                                                   size = base_size * 1.3, margin = ggplot2::margin(b = 8)),
            plot.subtitle = ggplot2::element_text(hjust = 0.5, face = "italic", 
                                                   size = base_size * 0.95),
            axis.title    = ggplot2::element_text(size = base_size * 1.1),
            axis.text     = ggplot2::element_text(size = base_size * 0.9),
            panel.grid.minor = ggplot2::element_blank(),
            panel.border  = ggplot2::element_rect(color = "grey85", fill = NA, linewidth = 0.3)
        )
    
    if (type == "spectrum") {
        t <- t + ggplot2::theme(
            legend.position     = "right",
            panel.grid.major.y  = ggplot2::element_line(color = "gray90", linewidth = 0.25),
            plot.margin         = ggplot2::margin(t = 5, r = 8, b = 5, l = 5, unit = "mm")
        )
    }
    t
}
```

Then update `.apply_publication_theme()` to pass `type` through.

------------------------------------------------------------------------

### 5.3 — Y-Axis Label Standardization

Use LaTeX-style [`expression()`](https://rdrr.io/r/base/expression.html)
notation **consistently** across all spectrum plots:

| Plot type | Unified label |
|----|----|
| Tsallis entropy (all modes) | `expression("Tsallis entropy (" * S[q] * ")")` |
| Divergence (all modes) | `expression("Divergence " * D[q])` |

Gene-specific and global modes must use the same label format.

------------------------------------------------------------------------

### 5.4 — Remove `expression(bold(...))` Duplication

Since `.apply_publication_theme()` already sets `face = "bold"` on plot
titles, remove all `expression(bold(...))` calls. Use plain
`expression(...)` or character strings instead.

**Files affected**: `plots_spectrum.R:309, 350`,
`plots_transcript.R:1673`.

------------------------------------------------------------------------

### 5.5 — Eliminate Duplicate Functions

| Remove (old, concatenated name)          | Keep (new, proper name)          |
|------------------------------------------|----------------------------------|
| `.make_plot_for_geneprepare_inputs()`    | `.prepare_transcript_inputs()`   |
| `.make_plot_for_genebuild_tx_long()`     | `.build_transcript_long()`       |
| `.make_plot_for_geneaggregate_df_long()` | `.aggregate_transcript_data()`   |
| `.make_plot_for_generead_tx2gene()`      | `.read_tx2gene()`                |
| `.make_plot_for_genemake_agg()`          | `.create_aggregation_function()` |
| `.apply_tsenat_theme()`                  | `.apply_publication_theme()`     |
| `.set_plot_title()`                      | `.apply_publication_theme()`     |

------------------------------------------------------------------------

### 5.6 — Move `%||%` to `globals.R`

``` r

# In R/globals.R:
`%||%` <- function(x, y) {
    if (is.null(x) || (length(x) == 1 && is.na(x))) y else x
}
```

Remove from `R/orchestration.R:441`.

------------------------------------------------------------------------

### 5.7 — Standardize heatmap fontsize defaults

| Function                   | Proposed default                       |
|----------------------------|----------------------------------------|
| `.plot_jis_delta()`        | `fontsize = 12`                        |
| `.plot_expression()`       | `fontsize = 12`                        |
| `.create_pheatmap_grob()`  | `fontsize = 12`                        |
| `.create_tsenat_heatmap()` | `fontsize_row = 11, fontsize_col = 11` |

Use `.font_sizes$heatmap_main` and `.font_sizes$heatmap_labels` as the
canonical source.

------------------------------------------------------------------------

## 6. 📋 Prioritized Action Items

| \# | Action | Impact | Effort | Files |
|----|----|----|----|----|
| **P1** | Remove dangling roxygen in `plots_validation.R:152-156` | Low | 1 min | `plots_validation.R` |
| **P1** | Remove dead `if (FALSE)` debug block in `plots_heatmaps.R:393` | Low | 1 min | `plots_heatmaps.R` |
| **P1** | Fix duplicated `has_valid_ci` in `.create_ci_ribbon_plot()` | Low | 5 min | `plots_themes.R` |
| **P1** | Move `%||%` to `globals.R` | Medium | 5 min | `orchestration.R`, `globals.R` |
| **P2** | Normalize `.font_sizes` values | **High** | 10 min | `globals.R` |
| **P2** | Remove all hardcoded `title_size`/`subtitle_size` overrides | **High** | 1–2 hr | All 10 files |
| **P2** | Merge `.theme_spectrum` into `.theme_base(type=)` | Medium | 30 min | `globals.R`, `plots_themes.R`, callers |
| **P3** | Unify Y-axis labels (expression notation everywhere) | Medium | 1 hr | `plots_tsallis_q.R`, `plots_spectrum.R` |
| **P3** | Remove `expression(bold(...))` from titles | Low | 30 min | `plots_spectrum.R`, `plots_transcript.R` |
| **P3** | Delete deprecated functions (old concatenated names) | Medium | 1 hr | `plots_transcript.R` |
| **P3** | Delete `.apply_tsenat_theme()` and `.set_plot_title()` (unused) | Low | 5 min | `plots_themes.R` |
| **P4** | Standardize heatmap fontsize defaults | Low | 15 min | `plots_heatmaps.R`, `plots_themes.R` |
| **P4** | Fix `.plot_divergence_distribution()` to use `.apply_publication_theme()` | Low | 10 min | `plots_transcript.R` |
| **P4** | Change `cellwidth=0` sentinel to `NULL` | Low | 20 min | `plots_heatmaps.R` |

------------------------------------------------------------------------

## 7. Summary

### What’s Working Well ✅

The July 2026 refactoring (I11) made significant progress: -
**Centralized palette dispatch** — `.tsenat_palette()` unified 6
scattered palette functions - **`.apply_publication_theme()`** —
consolidated the 35+ line theme pattern into one call -
**`.configure_legend()`** — standardized legend across 16+ call sites -
**`.create_ci_ribbon_plot()`** — consolidated the ribbon+line+point
pattern (18x occurrences) - **`.assemble_grid_plot()`** — centralized
grid assembly with legend extraction - **`.save_plot_standard()`** —
unified plot saving (20+ occurrences → one function)

### What Needs Work ⚠️

| Issue                                 | Files affected       | Severity |
|---------------------------------------|----------------------|----------|
| Font size inconsistency (title 14–21) | All 10 files         | **HIGH** |
| Theme selection inconsistency         | 5 files              | MEDIUM   |
| Y-axis label format mismatch          | 2 files              | MEDIUM   |
| Duplicate functions (old + new)       | `plots_transcript.R` | MEDIUM   |
| Concatenated function names           | `plots_transcript.R` | LOW      |
| 5 minor bugs                          | 3 files              | LOW      |

### Bottom Line

The plotting codebase is in **good structural shape** after the July
2026 refactoring — the architecture is sound and the centralized helpers
are well-designed. The remaining issues are **cosmetic consistency
problems**, not structural flaws. The single highest-impact change would
be: **normalize `.font_sizes` values and remove all hardcoded font size
overrides**, forcing every plot through `.apply_publication_theme()`
with the global constants. This would eliminate the 50% title size
spread and ensure visual consistency across all plot types.

------------------------------------------------------------------------

*Report generated 2026-07-27. Analysis covered all 10 `plots_*.R` files
totaling 7,018 lines.*
