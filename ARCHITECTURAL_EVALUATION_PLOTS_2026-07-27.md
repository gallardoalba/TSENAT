# Architectural Evaluation: TSENAT plots\_\*.R — Centralization Opportunities

**Date**: 2026-07-27  
**Scope**: 10 files, 6,279 lines  
**Status**: Evaluation Complete ✅

------------------------------------------------------------------------

## Table of Contents

1.  [Current Architecture Map](#id_1-current-architecture-map)
2.  [Identified Anti-Patterns](#id_2-identified-anti-patterns)
3.  [Centralization Opportunities](#id_3-centralization-opportunities)
4.  [Proposed Target Architecture](#id_4-proposed-target-architecture)
5.  [Implementation Roadmap](#id_5-implementation-roadmap)

------------------------------------------------------------------------

## 1. Current Architecture Map

### 1.1 — Exported (Public) Functions

| Function | File | Purpose |
|----|----|----|
| [`plot_diversity_spectrum()`](https://gallardoalba.github.io/TSENAT/reference/plot_diversity_spectrum.md) | `plots_tsallis_q.R` | Tsallis entropy q-curve (aggregate + gene-specific + bootstrap CI) |
| [`plot_diversity_violin_density()`](https://gallardoalba.github.io/TSENAT/reference/plot_diversity_violin_density.md) | `plots_violin_density.R` | Violin + density combo for single q-value |

### 1.2 — Internal Dispatch Functions

| Function | File | Lines | Purpose |
|----|----|----|----|
| `.plot_divergence_spectrum()` | `plots_spectrum.R` | 397 | Divergence q-curve (single gene / top genes / global) |
| `.plot_sait()` | `plots_sait.R` | 157 | GAM q-curve grid for top SAIT genes |
| `.plot_jis_delta()` | `plots_heatmaps.R` | ~240 | Multi-q delta influence heatmaps (5-phase pipeline) |
| `.plot_expression()` | `plots_heatmaps.R` | ~160 | Per-gene expression heatmaps (5-phase pipeline, duplicated) |
| `.plot_divergence_distribution()` | `plots_transcript.R` | ~50 | Histogram of divergence effect sizes |
| `.plot_diversity_density_singleq()` | `plots_violin_density.R` | ~45 | Single-q density plot |
| `.plot_diversity_violin_singleq()` | `plots_violin_density.R` | ~45 | Single-q violin plot |

### 1.3 — Helper Functions by Layer

    ┌─────────────────────────────────────────────────────────┐
    │  OUTPUT LAYER (plots_themes.R)                          │
    │  .save_plot_standard()    .assemble_grid_plot()         │
    │  .create_title_grob()     .create_tsenat_heatmap()      │
    ├─────────────────────────────────────────────────────────┤
    │  STYLING LAYER (plots_themes.R)                         │
    │  .apply_publication_theme()  .apply_group_aesthetics()  │
    │  .configure_legend()         .apply_facet_styling()     │
    │  .add_reference_lines()      .create_centered_theme()   │
    │  .create_ci_ribbon_plot()    .create_simple_line_plot() │
    │  .format_axis_labels()       .calculate_scaled_fonts()  │
    ├─────────────────────────────────────────────────────────┤
    │  DATA PREP LAYER (scattered across files)               │
    │  .prepare_long_format()      .prepare_tsallis_long()    │
    │  .extract_q_value()          .prepare_grouped_long()    │
    │  .prepare_transcript_inputs() .build_transcript_long()  │
    │  .aggregate_transcript_data() .compute_diversity_spectrum()│
    │  .compute_distribution_stats() .extract_bootstrap_ci_assays()│
    ├─────────────────────────────────────────────────────────┤
    │  GENE SELECTION LAYER (scattered, duplicated)           │
    │  .select_top_genes()         .filter_genes_by_pvalue()  │
    │  .heatmap_select_genes_multiq() .heatmap_select_genes_results()│
    │  .plot_select_genes()        .select_genes_from_results()│
    │  .spectrum_find_gene_column() .spectrum_find_pvalue_column()│
    ├─────────────────────────────────────────────────────────┤
    │  VALIDATION LAYER (scattered)                           │
    │  .validate_diversity_se()    .validate_results_df()     │
    │  .validate_multiq_input()    .validate_se_for_heatmaps()│
    │  .spectrum_validate_inputs() .validate_plot_data()      │
    │  .validate_control_in_samples()                         │
    ├─────────────────────────────────────────────────────────┤
    │  GRID COMPOSITION (duplicated in 3 places)              │
    │  .combine_plots_patchwork()   .combine_plots_cowplot()  │
    │  .combine_plots_grid()        .combine_gene_plots()     │
    │  .plot_gam_arrange_grid()     .render_heatmaps_to_grid()│
    └─────────────────────────────────────────────────────────┘

------------------------------------------------------------------------

## 2. Identified Anti-Patterns

### AP1 — **Every plot function reimplements mode dispatch independently**

Each main plot function has its own
`if (gene) ... else if (sait_res) ... else ...` dispatch logic with
different parameter names for the same concepts:

| Concept | `plot_diversity_spectrum` | `.plot_divergence_spectrum` | `.plot_sait` | `.plot_expression` |
|----|----|----|----|----|
| Select genes by p-value | `sait_res`, `n_top` | `sait_res`, `n_genes` | `sait_res`, `n_top`, `sig_alpha` | `res`, `top_n` |
| Specific genes | `gene` | `gene` | `genes` | `gene` |
| Grid columns | (hardcoded 2) | `ncol` | (hardcoded 2) | `layout_ncol` |
| Grid title | (hardcoded) | (hardcoded) | (hardcoded) | (hardcoded) |

**Impact**: Adding a new “select genes by effect size” mode requires
changes in 4 separate dispatch blocks.

------------------------------------------------------------------------

### AP2 — **The 5-phase heatmap pipeline is duplicated**

`.plot_jis_delta()` and `.plot_expression()` both implement:

    Phase 1: Validate input
    Phase 2: Select genes
    Phase 3: Plan layout
    Phase 4: Create heatmaps (pheatmap grobs)
    Phase 5: Render grid + finalize

These share helper functions (`.validate_multiq_input`,
`.plot_adaptive_layout`, etc.) but the **orchestration code** (the loop
that calls them) is copy-pasted with only parameter name differences.

------------------------------------------------------------------------

### AP3 — **TSENATAnalysis → SE conversion duplicated**

Both
[`plot_diversity_spectrum()`](https://gallardoalba.github.io/TSENAT/reference/plot_diversity_spectrum.md)
and
[`plot_diversity_violin_density()`](https://gallardoalba.github.io/TSENAT/reference/plot_diversity_violin_density.md)
contain nearly identical blocks:

``` r

if (methods::is(se, "TSENATAnalysis")) {
    if (length(se@diversity_results) == 0) {
        stop("No diversity results found...")
    }
    se <- se@diversity_results[[1]]
}
```

One also has condition_col extraction from `@config`, the other doesn’t.
No single entry point normalizes `TSENATAnalysis` input.

------------------------------------------------------------------------

### AP4 — **Gene selection functions scattered across 4 files**

| Function | File | Signature |
|----|----|----|
| `.select_top_genes()` | `plots_stats.R` | `(results, p_col, gene_col, n_genes)` |
| `.filter_genes_by_pvalue()` | `plots_stats.R` | `(results, p_threshold, p_col, gene_col)` |
| `.plot_select_genes()` | `plots_gam_helpers.R` | `(sait_res, genes, n_top, sig_alpha)` |
| `.select_genes_from_results()` | `plots_transcript.R` | `(res, top_n)` |
| `.heatmap_select_genes_multiq()` | `plots_heatmaps.R` | `(switching_results, n_genes, sait_results)` |
| `.heatmap_select_genes_results()` | `plots_heatmaps.R` | `(se, res, gene_col, top_n, tx2gene)` |
| `.spectrum_find_gene_column()` | `plots_spectrum.R` | `(df)` — only finds column name |
| `.spectrum_find_pvalue_column()` | `plots_spectrum.R` | `(df)` — only finds column name |

**6 different functions** do essentially the same thing: “give me the
top N genes by p-value from a results data frame.” Each has subtly
different column name auto-detection logic and return types.

------------------------------------------------------------------------

### AP5 — **Grid assembly implemented 5 different ways**

| Approach | Location | Mechanism |
|----|----|----|
| `.assemble_grid_plot()` | `plots_themes.R` | cowplot with legend extraction, title grob |
| `.plot_gam_arrange_grid()` | `plots_gam_helpers.R` | Manual cowplot with margin tweaks |
| `.combine_plots_patchwork()` | `plots_themes.R` | patchwork with row spacers |
| `.combine_plots_cowplot()` | `plots_themes.R` | cowplot with title+subtitle+spacer |
| `.render_heatmaps_to_grid()` | `plots_heatmaps.R` | grid viewports (low-level) |

`.assemble_grid_plot()` was created to consolidate this but is only used
by `plots_tsallis_q.R`. The GAM and heatmap code uses separate
implementations.

------------------------------------------------------------------------

### AP6 — **Data preparation has two parallel APIs**

| Function | What it does |
|----|----|
| `.prepare_tsallis_long()` | Basic long-format conversion |
| `.prepare_long_format()` | Wrapper: adds validation + auto-detect condition_col |

Some callers use the wrapper, some use the raw function. This means some
code paths skip validation. The `condition_col` auto-detection logic in
`.prepare_long_format()` is also duplicated inline in
[`plot_diversity_spectrum()`](https://gallardoalba.github.io/TSENAT/reference/plot_diversity_spectrum.md).

------------------------------------------------------------------------

### AP7 — **`condition_col` resolution is duplicated**

The logic “if condition_col is NULL, try config, then fall back to
‘condition’ or ‘sample_type’” appears in:

1.  [`plot_diversity_spectrum()`](https://gallardoalba.github.io/TSENAT/reference/plot_diversity_spectrum.md)
    — lines 145-163
2.  `.prepare_long_format()` — lines 795-803
3.  [`plot_diversity_violin_density()`](https://gallardoalba.github.io/TSENAT/reference/plot_diversity_violin_density.md)
    — only handles TSENATAnalysis, missing SE fallback

------------------------------------------------------------------------

## 3. Centralization Opportunities

### C1 — **Unified Plot Pipeline**

Every plot in TSENAT follows the same abstract pipeline:

    Input → Validate → Select Genes → Prepare Data → Construct Plot → Style → Assemble Grid → Output

A unified pipeline abstraction would look like:

``` r

.tsenat_plot_pipeline <- function(se, config = list()) {
    # config is a list specifying each pipeline stage:
    #   $mode: "aggregate" | "gene_specific" | "top_genes"
    #   $gene_selection: list(method = "pvalue" | "manual", ...)
    #   $plot_type: "q_curve" | "heatmap" | "violin" | "density"
    #   $grid: list(ncol = 2, title = "...")
    #   $output: list(file = NULL, width = 12)
    
    se <- .normalize_input(se)            # C3: TSENATAnalysis → SE
    se <- .validate_plot_input(se)        # C4: Unified validation
    genes <- .resolve_genes(se, config)   # C5: Unified gene selection
    data <- .prepare_plot_data(se, genes) # C6: Unified data prep
    plots <- .build_plots(data, config)   # Per-plot-type construction
    styled <- .apply_tsenat_style(plots)  # C7: Unified styling
    grid <- .assemble_output(styled)      # C8: Unified grid assembly
    .finalize_output(grid, config)        # C9: Unified save/return
}
```

**Benefit**: Every plot type shares validation, gene selection, styling,
and output. Only the “construct plot” stage differs.

------------------------------------------------------------------------

### C2 — **Unified Gene Selection Module**

Consolidate 6 gene selection functions into one:

``` r

.resolve_plot_genes <- function(se, results = NULL, genes = NULL, 
    n_top = 4, rank_by = c("pvalue", "effect_size", "variance"),
    p_col = NULL, gene_col = NULL, sig_threshold = NULL) {
    
    # Case 1: User provided specific genes
    if (!is.null(genes)) return(.validate_gene_list(genes, se))
    
    # Case 2: Select from results by ranking
    if (!is.null(results)) {
        return(.rank_and_select(results, n_top, rank_by, p_col, gene_col, sig_threshold))
    }
    
    # Case 3: Use all genes (aggregate mode)
    return(NULL)  # NULL means "all genes"
}
```

**Impact**: Eliminates `.plot_select_genes()`,
`.heatmap_select_genes_multiq()`, `.heatmap_select_genes_results()`,
`.select_genes_from_results()`, `.spectrum_find_gene_column()`,
`.spectrum_find_pvalue_column()` — replaces 6 functions with 1.

------------------------------------------------------------------------

### C3 — **Unified Input Normalizer**

``` r

.normalize_plot_input <- function(se, assay_name = "diversity", condition_col = NULL) {
    # Handle TSENATAnalysis
    if (methods::is(se, "TSENATAnalysis")) {
        condition_col <- condition_col %||% se@config$condition_col %||% "condition"
        se <- .prepare_combined_se(se)
    }
    
    # Handle SummarizedExperiment
    if (!methods::is(se, "SummarizedExperiment")) {
        stop("Input must be a SummarizedExperiment or TSENATAnalysis object")
    }
    
    # Resolve condition column
    condition_col <- .resolve_condition_col(se, condition_col)
    
    # Validate assay exists
    if (!assay_name %in% SummarizedExperiment::assayNames(se)) {
        stop("Assay '", assay_name, "' not found")
    }
    
    list(se = se, condition_col = condition_col)
}
```

**Impact**: Removes 30+ lines of duplicated TSENATAnalysis handling from
[`plot_diversity_spectrum()`](https://gallardoalba.github.io/TSENAT/reference/plot_diversity_spectrum.md)
and
[`plot_diversity_violin_density()`](https://gallardoalba.github.io/TSENAT/reference/plot_diversity_violin_density.md).

------------------------------------------------------------------------

### C4 — **Unified Validation Layer**

Currently validation functions are scattered. Centralize into a single
dispatch:

``` r

.validate_plot_input <- function(se, plot_type = c("q_curve", "divergence", "heatmap", "violin")) {
    plot_type <- match.arg(plot_type)
    
    # Common validations
    if (nrow(se) == 0) stop("No genes in input")
    if (ncol(se) == 0) stop("No samples in input")
    
    # Plot-specific validations
    switch(plot_type,
        q_curve = {
            if (length(unique(SummarizedExperiment::colData(se)$condition)) < 2)
                stop("Need at least 2 conditions for q-curve comparison")
        },
        heatmap = {
            if (!"gene_name" %in% colnames(SummarizedExperiment::rowData(se)))
                stop("rowData must have 'gene_name' column for heatmaps")
        },
        violin = {},  # No additional checks
        divergence = {}
    )
    
    invisible(se)
}
```

------------------------------------------------------------------------

### C5 — **Unified Data Preparation**

Currently there are two parallel APIs (`.prepare_tsallis_long()` and
`.prepare_long_format()`). Unify:

``` r

.prepare_plot_data <- function(se, assay_name = "diversity", condition_col = NULL,
    format = c("long", "wide", "matrix"), genes = NULL) {
    
    format <- match.arg(format)
    
    # Always go through the validated wrapper
    long <- .prepare_long_format(se, assay_name, condition_col)
    
    # Filter to genes if specified
    if (!is.null(genes)) {
        long <- long[long$Gene %in% genes, ]
    }
    
    # Convert to requested format
    switch(format,
        long = long,
        wide = tidyr::pivot_wider(long, names_from = q, values_from = tsallis),
        matrix = as.matrix(SummarizedExperiment::assay(se, assay_name))
    )
}
```

------------------------------------------------------------------------

### C6 — **Unified Grid Assembly**

Replace 5 grid assembly implementations with one:

``` r

.assemble_plot_grid <- function(plots, config = list()) {
    # config: ncol, title, subtitle, legend_position, output_file
    
    if (length(plots) == 1) {
        return(.apply_publication_aesthetics(plots[[1]], 
            title = config$title, subtitle = config$subtitle))
    }
    
    # Multi-plot: use the existing .assemble_grid_plot()
    .assemble_grid_plot(
        plots,
        ncol = config$ncol %||% 2,
        title = config$title,
        subtitle = config$subtitle,
        legend_position = config$legend_position %||% "bottom",
        extract_legend = config$extract_legend %||% TRUE
    )
}
```

------------------------------------------------------------------------

### C7 — **Unified Output Handler**

``` r

.finalize_plot <- function(plot, output_file = NULL, width = 12, aspect = "standard") {
    if (!is.null(output_file)) {
        .save_plot_standard(plot, output_file, width_inches = width, aspect_type = aspect)
    }
    invisible(plot)
}
```

------------------------------------------------------------------------

## 4. Proposed Target Architecture

### 4.1 — File Reorganization

    R/
    ├── plots_pipeline.R        # NEW: Unified pipeline + input normalization + validation
    ├── plots_gene_selection.R  # NEW: Single gene selection module
    ├── plots_data_prep.R       # NEW: Merged plots_stats.R + data prep functions
    ├── plots_themes.R          # KEEP: Styling, themes, palettes, grid assembly
    ├── plots_q_curve.R         # MERGE: plots_tsallis_q.R + plots_spectrum.R
    ├── plots_gam.R             # MERGE: plots_sait.R + plots_gam_helpers.R
    ├── plots_heatmaps.R        # KEEP: Refactored with shared pipeline
    ├── plots_distribution.R    # MERGE: plots_violin_density.R + .plot_divergence_distribution
    ├── plots_transcript.R      # KEEP: Transcript-level expression
    └── plots_validation.R      # KEEP: Validation + formatting utilities

### 4.2 — Call Flow

    User calls exported function (e.g., plot_diversity_spectrum)
      │
      ▼
    .normalize_plot_input(se)          # TSENATAnalysis → SE, resolve condition_col
      │
      ▼
    .validate_plot_input(se, "q_curve") # Plot-specific validation
      │
      ▼
    .resolve_plot_genes(se, results, genes, n_top)  # Unified gene selection
      │
      ▼
    .prepare_plot_data(se, format = "long", genes)   # Unified data prep
      │
      ▼
    .build_q_curve(data)               # Plot-type-specific construction
      │
      ▼
    .apply_publication_theme(plot)     # Unified styling
      │
      ▼
    .assemble_plot_grid(plots)         # Unified grid (single-plot: no-op)
      │
      ▼
    .finalize_plot(plot, output_file)  # Unified save/return

### 4.3 — Exported Functions Become Thin Wrappers

``` r

#' @export
plot_diversity_spectrum <- function(se, ...) {
    .tsenat_plot_pipeline(se, plot_type = "q_curve", ...)
}

#' @export
plot_diversity_violin_density <- function(se, ...) {
    .tsenat_plot_pipeline(se, plot_type = "violin_density", ...)
}

# Future additions become trivial:
#' @export
plot_divergence_heatmap <- function(se, ...) {
    .tsenat_plot_pipeline(se, plot_type = "divergence_heatmap", ...)
}
```

------------------------------------------------------------------------

## 5. Implementation Roadmap

### Phase 1: Foundation (non-breaking, additive)

| Step | Action | Files | Effort |
|----|----|----|----|
| 1.1 | Create `.normalize_plot_input()` in `plots_themes.R` | `plots_themes.R` | 1 hr |
| 1.2 | Create `.resolve_plot_genes()` in new `plots_gene_selection.R` | new file | 2 hr |
| 1.3 | Create `.finalize_plot()` in `plots_themes.R` | `plots_themes.R` | 30 min |
| 1.4 | Add deprecation warnings to old gene selection functions | 4 files | 1 hr |

### Phase 2: Migration (incremental, each plot independently)

| Step | Action | Impact |
|----|----|----|
| 2.1 | Refactor [`plot_diversity_spectrum()`](https://gallardoalba.github.io/TSENAT/reference/plot_diversity_spectrum.md) to use unified pipeline | ~60 lines removed |
| 2.2 | Refactor [`plot_diversity_violin_density()`](https://gallardoalba.github.io/TSENAT/reference/plot_diversity_violin_density.md) to use unified pipeline | ~40 lines removed |
| 2.3 | Refactor `.plot_divergence_spectrum()` to use unified pipeline | ~50 lines removed |
| 2.4 | Refactor `.plot_sait()` to use unified pipeline | ~40 lines removed |
| 2.5 | Refactor `.plot_jis_delta()` and `.plot_expression()` to share pipeline | ~100 lines removed |

### Phase 3: Consolidation (cleanup)

| Step | Action |
|----|----|
| 3.1 | Delete old gene selection functions (after tests updated) |
| 3.2 | Merge `plots_tsallis_q.R` + `plots_spectrum.R` → `plots_q_curve.R` |
| 3.3 | Merge `plots_sait.R` + `plots_gam_helpers.R` → `plots_gam.R` |
| 3.4 | Merge `plots_violin_density.R` + `.plot_divergence_distribution` → `plots_distribution.R` |

------------------------------------------------------------------------

## Summary

| Metric                        | Current         | Target        |
|-------------------------------|-----------------|---------------|
| Gene selection functions      | 6 (scattered)   | 1             |
| Grid assembly implementations | 5               | 1             |
| Input normalization blocks    | 3 (duplicated)  | 1             |
| Data preparation APIs         | 2 (parallel)    | 1             |
| Mode dispatch patterns        | 4 (independent) | 1 (pipeline)  |
| Exported function boilerplate | ~80 lines each  | ~5 lines each |

**Estimated total reduction**: 300–500 lines of duplicated orchestration
code across the 7 plot dispatch functions.

**Key principle**: The current architecture already has well-factored
**styling** and **output** layers (plots_themes.R). What’s missing is a
unified **orchestration** layer that sits above them and handles: input
normalization → validation → gene selection → data preparation → routing
to plot-specific builders. Adding this layer would make the system
dramatically easier to extend and maintain.

------------------------------------------------------------------------

*Report generated 2026-07-27.*
