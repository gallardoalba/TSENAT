# Map external coldata into a SummarizedExperiment

This helper maps an external \`coldata\` table (with sample IDs and a
condition/label column) into a \`SummarizedExperiment\` produced by
\`calculate_diversity()\`. It assigns a \`sample_type\` column to
\`colData(ts_se)\` and falls back to \`infer_sample_group()\` for
unmapped samples.

## Usage

``` r
map_coldata_to_se(
  ts_se,
  coldata,
  coldata_sample_col = "Sample",
  coldata_condition_col = "Condition"
)
```

## Arguments

- ts_se:

  A \`SummarizedExperiment\` object with assay columns named possibly
  including \`\_q=\` suffixes.

- coldata:

  A data.frame with sample metadata or \`NULL\`.

- coldata_sample_col:

  Name of the column in \`coldata\` with sample IDs.

- coldata_condition_col:

  Name of the column in \`coldata\` with condition/labels.

## Value

The input \`ts_se\` with \`colData(ts_se)\$sample_type\` set when
possible.
