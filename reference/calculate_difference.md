# Calculate Difference Between Control and Treatment Groups (S4 Wrapper)

S4 wrapper that operates on TSENATAnalysis objects to calculate
differences between control and treatment groups. Uses diversity results
from `@diversity_results` slot (from
[`calculate_diversity()`](https://gallardoalba.github.io/TSENAT/reference/calculate_diversity.md))
and stores results in the `lm_results` slot.

## Usage

``` r
calculate_difference(
  analysis,
  control = NULL,
  q = NULL,
  condition_col = NULL,
  method = NULL,
  test = NULL,
  randomizations = NULL,
  pcorr = NULL,
  assayno = NULL,
  verbose = NULL,
  paired = FALSE,
  exact = FALSE,
  pseudocount = NULL,
  nthreads = NULL,
  robust_loss_type = NULL,
  robust_scale_method = NULL,
  pairs = NULL,
  output_file = NULL,
  ...
)
```

## Arguments

- analysis:

  A `TSENATAnalysis` object with diversity results in
  `@diversity_results`.

- control:

  Character string specifying the control group identifier. If `NULL`,
  attempts to retrieve from `analysis@config$control`.

- q:

  `numeric` or `NULL`. Q-value to use for testing. If NULL: auto-detects
  from diversity results if only ONE q-value present, or errors if
  MULTIPLE q-values present (must specify which to test). NOTE: q is NOT
  read from config - only explicit argument or auto-detected from
  diversity. If diversity contains single q-value, q is OPTIONAL. If
  diversity contains multiple q-values, q is REQUIRED.

- condition_col:

  `character` or `NULL`. Column name in colData identifying sample
  conditions. If NULL, reads from `@config$condition_col` or
  auto-detects.

- method:

  `character`. Difference calculation method. Default: 'mean'. If NULL,
  reads from `@config$method` if available.

- test:

  `character`. Statistical test type. Default: 'wilcoxon'. If NULL,
  reads from `@config$test` if available.

- randomizations:

  `numeric`. Number of randomizations. Default: 100. If NULL, reads from
  `@config$randomizations` if available.

- pcorr:

  `character`. P-value correction method. Default: 'BH'. If NULL, reads
  from `@config$pcorr` if available.

- assayno:

  `numeric`. Assay number to use. Default: 1. If NULL, reads from
  `@config$assayno` if available.

- verbose:

  `logical`. Print progress messages. Default: TRUE. If not specified,
  reads from `@config$verbose` if available.

- paired:

  `logical`. Whether data is paired. Default: FALSE. If not specified,
  reads from `@config$paired` if available.

- exact:

  `logical`. Use exact test. Default: FALSE. If not specified, reads
  from `@config$exact` if available.

- pseudocount:

  `numeric`. Pseudocount for normalization. Default: 0. If NULL, reads
  from `@config$pseudocount` if available.

- nthreads:

  `numeric` or `NULL`. Number of CPU threads for parallel processing. If
  NULL, reads from `@config$nthreads` (or defaults to 1).

- robust_loss_type:

  `character`. Robust regression loss type. Default: 'huber'. If NULL,
  reads from `@config$robust_loss_type` if available.

- robust_scale_method:

  `character`. Robust scaling method. Default: 'mad'. If NULL, reads
  from `@config$robust_scale_method` if available.

- pairs:

  `character` or `numeric` vector or `NULL`. Pairing information for
  paired designs. When `paired = TRUE`, specifies which samples are
  paired (e.g., c(1,1,2,2,3,3) for 3 pairs). Default: NULL. When NULL
  with `paired = TRUE`, auto-extracted from colData 'sample_base' column
  if available.

- output_file:

  `character` or `NULL`. Optional file path to save results. Supported
  formats: .rds (for S4 objects), .tsv, .csv, .txt (for tables).
  Default: NULL (no file output).

- ...:

  Additional arguments passed to the base function.

## Value

Returns the modified `analysis` object invisibly with results stored in
`analysis@pairwise_results$difference`.

## Details

\*\*IMPORTANT: \*\* Requires diversity results to exist first via
[`calculate_diversity()`](https://gallardoalba.github.io/TSENAT/reference/calculate_diversity.md).
This wrapper extracts the diversity SummarizedExperiment from
`@diversity_results`, not the raw input data in `@se`. This ensures
you're comparing diversity values between control and treatment groups,
not raw abundance data.

\*\*Parameter resolution priority\*\* (explicit \> @config \>
auto-detect \> error):

- `control`: Uses explicit arg, else `@config$control`, else error

- `nthreads`: Uses explicit arg, else `@config$nthreads`, else 1

- `condition_col` (sample grouping): Uses `@config$condition_col`, else
  auto-detects from colData columns: 'group', 'sample_type', 'condition'

## See also

[`calculate_diversity`](https://gallardoalba.github.io/TSENAT/reference/calculate_diversity.md)
for computing diversity.

## Examples

``` r
# Load example data (matching TSENAT.Rmd workflow)
data(readcounts)
readcounts <- as.matrix(readcounts)
mode(readcounts) <- 'numeric'
metadata_df <- read.table(
  system.file('extdata', 'metadata.tsv', package = 'TSENAT'),
  header = TRUE, sep = '\t'
)
gff3_dataset <- system.file('extdata', 'annotation.gff3.gz', package =
'TSENAT')

# Load TPM and effective length from vignette data
data(readcounts)  # Also loads tpm and effective_length

# Create config (metadata passed as explicit parameter to build_analysis)
config <- TSENAT_config(
  sample_col = 'sample',
  condition_col = 'condition',
  q = seq(0, 2, by = 0.05),
  paired = FALSE
)

# Build analysis from vignette data - metadata as explicit parameter
analysis <- build_analysis(
  readcounts = readcounts,
  metadata = metadata_df,
  tx2gene = gff3_dataset,
  config = config,
  tpm = tpm,
  effective_length = effective_length
)
analysis <- filter_analysis(analysis, min_samples = 1, subset_n_genes = 200)
analysis <- calculate_diversity(analysis, q = c(0.5, 1.0, 1.5))
#> Note: 37 genes excluded (< 75% valid values).
result <- calculate_difference(analysis, q = 1.0, control = 'normal')
```
