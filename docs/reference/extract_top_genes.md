# Extract significant genes and prepare plotting data

This function filters differential analysis results by significance,
extracts the top genes, and prepares sample information for plotting.
It's designed as a helper for transcript-level visualization workflows.

## Usage

``` r
extract_top_genes(res, ts_se, p_threshold = 0.05, top_n = 3)
```

## Arguments

- res:

  A data frame of differential analysis results (output from
  \`calculate_difference()\`), with columns including
  \`adjusted_p_values\` and \`genes\`.

- ts_se:

  A \`SummarizedExperiment\` object with diversity data, where
  \`colData(ts_se)\$sample_type\` contains sample group assignments.

- p_threshold:

  Numeric; adjusted p-value threshold for significance (default: 0.05).

- top_n:

  Integer; number of top genes to extract (default: 3).

## Value

A list containing: - \`top_genes\`: Character vector of the top
significant gene names. - \`samples_vec\`: Character vector of sample
types from colData. - \`sample_base_names\`: Character vector of sample
base names (q-stripped).
