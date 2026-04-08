# Extract Multi-Q Tabular Results with Q-Value Column

Helper to combine results from multiple q-values into a single
data.frame with q_value column.

## Usage

``` r
extract_multiq_table(
  result,
  is_multiq = NULL,
  extract_fn = NULL,
  q_value_col = "q_value"
)
```

## Arguments

- result:

  Result object (may be list with multi-q structure or single result)

- is_multiq:

  Logical. Whether result has multi-q structure. Default: auto-detect.

- extract_fn:

  Function to extract table from each result element. Signature:
  `function(result_element, q_key)`. Default: extracts 'summary_table'.

- q_value_col:

  Name of q-value column to add. Default: 'q_value'

## Value

data.frame with combined results and q_value column

## Examples

``` r
# Create sample multi-q result structure
q_result <- list(
  q_0.5 = list(summary_table = data.frame(
    gene_id = c('GENE1', 'GENE2'),
    entropy = c(1.2, 1.5),
    psi_mean = c(0.3, 0.7)
  )),
  q_1.0 = list(summary_table = data.frame(
    gene_id = c('GENE1', 'GENE2'),
    entropy = c(1.1, 1.4),
    psi_mean = c(0.32, 0.68)
  ))
)

# Extract and combine results across q-values
combined_results <- extract_multiq_table(q_result)
head(combined_results)
#>   gene_id entropy psi_mean q_value
#> 1   GENE1     1.2     0.30     0.5
#> 2   GENE2     1.5     0.70     0.5
#> 3   GENE1     1.1     0.32     1.0
#> 4   GENE2     1.4     0.68     1.0
```
