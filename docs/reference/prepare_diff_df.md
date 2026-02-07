# Prepare a differential data.frame from a SummarizedExperiment

This helper constructs the \`div_df\` used by downstream differential
functions and vignettes. It extracts the requested assay (default
\`diversity\`) and prepends a \`genes\` column taken from
\`rowData(se)\[\[gene_col\]\]\` if present, otherwise falls back to
\`rownames(assay)\`.

## Usage

``` r
prepare_diff_df(se, assay_name = "diversity", gene_col = "genes")
```

## Arguments

- se:

  A \`SummarizedExperiment\` produced by \`calculate_diversity()\`.

- assay_name:

  Character; name of the assay to extract (default: "diversity").

- gene_col:

  Character; column name in \`rowData(se)\` containing gene symbols
  (default: "genes").

## Value

A \`data.frame\` with a leading \`genes\` column and the assay
measurements as additional columns.
