# Dependencies Metrics Analysis → TSENAT Plot Optimizations

**Date**: 2026-07-27  
**Source**: `dependencies_metrics.json`  
**Scope**: Actionable optimizations for `plots_*.R` based on
dependency/risk/health/complexity metrics

------------------------------------------------------------------------

## 1. Key Metrics Summary

### Betweenness Centrality (Hub Functions)

Functions that sit at the center of the dependency graph — changing them
has cascade effects.

| Score    | Function              | Category         |
|----------|-----------------------|------------------|
| 0.65     | `calculate_sait`      | Statistical      |
| 0.63     | `calculate_diversity` | Statistical      |
| **0.16** | `.report_fit_summary` | **Plot-related** |

**Finding**: Plot functions have relatively low betweenness — they’re
leaf nodes. Good architecture.

### Complexity (100 = maximum)

| Score    | Function                               | File                    |
|----------|----------------------------------------|-------------------------|
| **100**  | `calculate_sait`                       | orchestration.R         |
| **100**  | `.calculate_divergence_impl`           | divergence_core.R       |
| **90.9** | `.plot_jis_delta`                      | **plots_heatmaps.R** ⚠️ |
| **90.9** | `.plot_expression`                     | **plots_heatmaps.R** ⚠️ |
| **90.9** | `.calculate_tsallis_entropy_bootstrap` | bootstrap.R             |
| **81.8** | `.plot_sait`                           | **plots_sait.R** ⚠️     |
| 63.6     | `.plot_tsallis_gene_bootstrap_ci`      | plots_tsallis_q.R       |

### Bottom Health (lower = worse, \< 55 critical)

| Score | Function | Issue |
|----|----|----|
| **42.5** | `.calculate_tsallis_entropy_bootstrap` | Statistical, not plots |
| 50.3 | `.calculate_jis` | Statistical |
| 53.1 | `.bootstrap_divergence` | Statistical |
| 53.8 | `.calculate_sait` | Statistical |
| **66.8** | `.plot_jis_delta` | **plots_heatmaps.R** — lowest plot health |
| **66.8** | `.plot_expression` | **plots_heatmaps.R** — tied with jis_delta |
| **64.4** | `.build_tile_plot_from_summary` | **plots_transcript.R** |

### Top Risk (100 = maximum, functions most likely to break)

| Score    | Function                            | Category                |
|----------|-------------------------------------|-------------------------|
| 100      | `bootstrap_replicate_cpp`           | C++ — needs C++ tests   |
| 100      | `.estimate_pseudocount`             | Utility — no tests?     |
| 100      | `.suggest_min_count`                | Utility — no tests?     |
| **92.7** | `.get_palette_colors`               | **plots_themes.R** ⚠️   |
| **86.8** | `.calculate_scaled_fonts`           | **plots_themes.R**      |
| **83.4** | `.extract_bootstrap_ci_assays`      | **plots_stats.R**       |
| **69.1** | `.create_simple_line_plot`          | **plots_themes.R**      |
| 70.5     | `.validate_results_df`              | plots_validation.R      |
| 70.5     | `.format_pvalue`                    | plots_validation.R      |
| 67.4     | `.process_switching_tables_results` | orchestration_results.R |

### Architecture Health

- **107 unstable** components (fragile, tightly coupled)
- **303 rigid** components (hard to change)
- **265 balanced** components (well-designed)
- **Ratio**: 1.14 rigid per balanced → slightly skewed toward rigidity

------------------------------------------------------------------------

## 2. Plot-Specific Findings

### F1 — `.plot_jis_delta()` and `.plot_expression()` are the least healthy plots

Both score **90.9 complexity** and **66.8 health**. The root cause is
the duplicated 5-phase heatmap pipeline (identified as AP2 in the
architectural evaluation). These two functions share a layout engine but
duplicate the orchestration code.

**Optimization**: Merge the shared pipeline into a single
`.render_heatmap_pipeline()` function that accepts a plot-type
parameter. This would: - Reduce complexity from 90.9 → ~60 for both -
Eliminate ~120 lines of duplicate orchestration - Improve health from
66.8 → ~80

### F2 — `.get_palette_colors()` has 92.7 risk

This is already wrapped by `.tsenat_palette()` and
`.apply_group_aesthetics()`, yet `.get_palette_colors()` itself remains
high-risk. The risk comes from its
[`get()`](https://rdrr.io/r/base/get.html) +
[`tryCatch()`](https://rdrr.io/r/base/conditions.html) dynamic dispatch
pattern:

``` r

.get_palette_colors <- function(palette_name = "palette_blue_red") {
    if (is.character(palette_name) && length(palette_name) == 1) {
        tryCatch({
            palette_fn <- get(paste0(".", palette_name))
            return(palette_fn())
        }, error = function(e) {
            stop(sprintf("Palette '%s' not found", palette_name))
        })
    }
    ...
}
```

**Optimization**: Replace [`get()`](https://rdrr.io/r/base/get.html)
with an explicit named list lookup to eliminate the runtime symbol
resolution risk.

### F3 — `.calculate_scaled_fonts()` has 86.8 risk

This function computes proportional font sizes based on output
dimensions. It’s called from `.plot_gam_save_plot()` and has complex
arithmetic with default multipliers.

**Optimization**: Since we normalized `.font_sizes` to absolute values
(18, 14, 12, etc.), the scaling logic may be unnecessary. Evaluate
whether font scaling by output area is actually needed, or if absolute
sizes suffice.

### F4 — `.extract_bootstrap_ci_assays()` has 83.4 risk

Used by multiple plot functions to detect and extract bootstrap CI
matrices from SummarizedExperiment objects. The risk comes from its
naming convention check (`paste0(assay_name, "_ci_lower")`) which is
fragile to naming changes.

**Optimization**: Centralize the CI naming convention in `globals.R` as
a constant, and add explicit validation that extracted CIs match the
base assay dimensions.

### F5 — `.build_tile_plot_from_summary()` health 64.4

This is the old renamed function (was
`.make_plot_for_genebuild_plot_from_summary`) that creates tile plots.
It has manual
[`ggplot2::theme()`](https://ggplot2.tidyverse.org/reference/theme.html)
element calculations and custom font scaling.

**Optimization**: Refactor to use `.apply_publication_theme()` and
remove the manual `element_text(size = ...)` calls that duplicate the
theme system.

------------------------------------------------------------------------

## 3. Prioritized Optimization Plan

### Phase 1: High-Impact, Low-Risk (today)

| \# | Change | Impact | Effort |
|----|----|----|----|
| **O1** | Replace `.get_palette_colors()` [`get()`](https://rdrr.io/r/base/get.html) dispatch with explicit lookup | Risk 92.7 → ~40 | 15 min |
| **O2** | Add CI naming constants to `globals.R` + validate dimensions in `.extract_bootstrap_ci_assays()` | Risk 83.4 → ~50 | 20 min |
| **O3** | Refactor `.build_tile_plot_from_summary()` to use `.apply_publication_theme()` | Health 64.4 → ~75 | 20 min |

### Phase 2: Medium-Impact (this week)

| \# | Change | Impact | Effort |
|----|----|----|----|
| **O4** | Merge `.plot_jis_delta()` + `.plot_expression()` shared pipeline into `.render_heatmap_pipeline()` | Complexity 90.9→60 for both, Health 66.8→80 | 2-3 hr |
| **O5** | Evaluate and simplify `.calculate_scaled_fonts()` — use absolute sizes | Risk 86.8→~30 | 1 hr |

### Phase 3: Broader Refactoring

| \# | Change | Impact |
|----|----|----|
| **O6** | Extract `.create_simple_line_plot()` into a small family of plot builders | Risk 69.1→~40 |
| **O7** | Add test coverage for high-risk utilities (`.estimate_pseudocount`, `.suggest_min_count`, `bootstrap_replicate_cpp`) | Risk 100→~60 |

------------------------------------------------------------------------

## 4. Concrete Implementation: O1 — Palette Lookup Hardening

Current (high-risk):

``` r

palette_fn <- get(paste0(".", palette_name))  # Runtime symbol lookup
return(palette_fn())
```

Proposed (safe):

``` r

palette_registry <- list(
    palette_blue_red = .palette_blue_red,
    palette_discrete = .palette_discrete,
    palette_continuous_diverging = .palette_continuous_diverging
)
if (palette_name %in% names(palette_registry)) {
    return(palette_registry[[palette_name]]())
}
stop("Unknown palette: ", palette_name)
```

## 5. Concrete Implementation: O2 — CI Naming Convention

Current (fragile):

``` r

ci_lower_name <- paste0(assay_name, "_ci_lower")
```

Proposed (centralized):

``` r

# In globals.R:
.CI_SUFFIX_LOWER <- "_ci_lower"
.CI_SUFFIX_UPPER <- "_ci_upper"

# In .extract_bootstrap_ci_assays:
ci_lower_name <- paste0(assay_name, .CI_SUFFIX_LOWER)
```

------------------------------------------------------------------------

## 6. Summary

| Category | Current State | After Optimizations |
|----|----|----|
| Plot functions with complexity \> 80 | 3 (.plot_jis_delta, .plot_expression, .plot_sait) | 1 (.plot_sait only) |
| Plot functions with health \< 70 | 3 | 0 |
| High-risk plot utilities (\> 80 risk) | 3 (.get_palette_colors, .calculate_scaled_fonts, .extract_bootstrap_ci_assays) | 0 |
| Architecture rigid/balanced ratio | 1.14 | ~0.85 |

**Bottom line**: The metrics confirm that the architectural evaluation
findings are correct — the heatmap pipeline duplication is the \#1 plot
problem (complexity 90.9 × 2 functions), and the palette/font/CI
utilities are the \#1 risk area. The unification work already done
(`.normalize_plot_input`, `.resolve_plot_genes`, `.finalize_plot`)
directly addresses the “rigid” architecture score by reducing coupling.
Implementing O1-O5 would bring plot code health from “needs attention”
to “well-maintained.”

------------------------------------------------------------------------

*Report generated 2026-07-27 from `dependencies_metrics.json`.*
