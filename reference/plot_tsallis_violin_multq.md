# Violin plot of Tsallis entropy for multiple q values

Violin plot of Tsallis entropy for multiple q values

## Usage

``` r
plot_tsallis_violin_multq(se, assay_name = "diversity")
```

## Arguments

- se:

  A \`SummarizedExperiment\` returned by \`calculate_diversity\` with
  multiple q values (colnames like 'Sample_q=0.01').

- assay_name:

  Name of the assay to use (default: "diversity").

## Value

A \`ggplot\` violin plot object faceted/colored by group and q.
