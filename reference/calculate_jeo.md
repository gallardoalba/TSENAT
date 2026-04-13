# Jackknife resampling with confidence intervals

Jackknife resampling with confidence intervals

## Usage

``` r
calculate_jeo(
  analysis,
  q = NULL,
  norm = NULL,
  log_base = NULL,
  top_n = NULL,
  verbose = NULL,
  nthreads = NULL,
  pseudocount = NULL,
  output_file = NULL,
  ...
)
```

## Arguments

- analysis:

  `TSENATAnalysis` object.

- q:

  `numeric`. Q-value(s) for jackknife estimation. If NULL, reads from
  `analysis@config$q`. Default: c(0, 0.5, 1, 1.5, 2) (matches
  `calculate_jis` spectrum).

- norm:

  `logical`. Normalization flag. Default: NULL (uses @config\$norm or
  TRUE).

- log_base:

  `numeric`. Logarithm base for entropy normalization. Default: NULL
  (uses e).

- top_n:

  `numeric`. Number of top outlier samples to report. Default: 5.

- verbose:

  `logical`. Print jackknife results summary. Default: FALSE.

- nthreads:

  `numeric` or `NULL`. Number of CPU threads for parallel processing. If
  NULL, reads from `@config$nthreads` (or defaults to 1). If \> 1 and
  multiple q-values provided, uses parallel PSOCK cluster.

- pseudocount:

  `numeric`. Pseudocount value for count regularization. Default: NULL
  (uses @config\$pseudocount or 0).

- output_file:

  `character` or `NULL`. Optional file path to save results. Supported
  formats: .tsv, .csv, .txt (for jackknife results table with estimates,
  influence, outliers), .rds (for entire S4 object). Default: NULL (no
  file output).

- ...:

  Additional arguments passed to the base function.

## Value

Modified TSENATAnalysis with jackknife results in @jackknife_results.

## Details

Requires diversity results to exist first. Will error if
[`calculate_diversity()`](https://gallardoalba.github.io/TSENAT/reference/calculate_diversity.md)
has not been run.

\*\*Parameter Priority Resolution:\*\*

- nthreads:

  Priority: explicit \> `@config$nthreads` \> 1

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

# TPM and effective_length REQUIRED for filter_analysis()
tpm <- matrix(runif(nrow(readcounts) * ncol(readcounts), 0.1, 10),
              nrow = nrow(readcounts), ncol = ncol(readcounts),
              dimnames = dimnames(readcounts))
effective_length <- matrix(100, nrow = nrow(readcounts), ncol = ncol(readcounts))

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

# Filter low-abundance genes (required for reliable jackknife estimates)
analysis <- filter_analysis(analysis, stringency = 'severe')

# Compute diversity first (required for jackknife)
analysis <- calculate_diversity(analysis, q = c(0.5, 1.0, 1.5))

# Run jackknife estimation
analysis <- calculate_jeo(analysis, q = c(0.5, 1.0, 1.5))
#> Warning: Total count (6.53) below recommended minimum (10-20).
#> Jackknife estimates may be unreliable (per papers S111, S114).
#> Consider aggregating samples or filtering genes with low abundance.
#> Warning: Total count (0) below recommended minimum (10-20).
#> Jackknife estimates may be unreliable (per papers S111, S114).
#> Consider aggregating samples or filtering genes with low abundance.
#> Warning: Gene unable to compute jackknife (diversity estimate is NA). Likely cause: zero or near-zero counts. Returning NA result.
#> Warning: Total count (0) below recommended minimum (10-20).
#> Jackknife estimates may be unreliable (per papers S111, S114).
#> Consider aggregating samples or filtering genes with low abundance.
#> Warning: Gene unable to compute jackknife (diversity estimate is NA). Likely cause: zero or near-zero counts. Returning NA result.
#> Warning: Total count (0) below recommended minimum (10-20).
#> Jackknife estimates may be unreliable (per papers S111, S114).
#> Consider aggregating samples or filtering genes with low abundance.
#> Warning: Gene unable to compute jackknife (diversity estimate is NA). Likely cause: zero or near-zero counts. Returning NA result.
#> Warning: Total count (0) below recommended minimum (10-20).
#> Jackknife estimates may be unreliable (per papers S111, S114).
#> Consider aggregating samples or filtering genes with low abundance.
#> Warning: Gene unable to compute jackknife (diversity estimate is NA). Likely cause: zero or near-zero counts. Returning NA result.
#> Warning: Total count (1) below recommended minimum (10-20).
#> Jackknife estimates may be unreliable (per papers S111, S114).
#> Consider aggregating samples or filtering genes with low abundance.
#> Warning: Total count (0) below recommended minimum (10-20).
#> Jackknife estimates may be unreliable (per papers S111, S114).
#> Consider aggregating samples or filtering genes with low abundance.
#> Warning: Gene unable to compute jackknife (diversity estimate is NA). Likely cause: zero or near-zero counts. Returning NA result.
#> Warning: Total count (0) below recommended minimum (10-20).
#> Jackknife estimates may be unreliable (per papers S111, S114).
#> Consider aggregating samples or filtering genes with low abundance.
#> Warning: Gene unable to compute jackknife (diversity estimate is NA). Likely cause: zero or near-zero counts. Returning NA result.
#> Warning: Total count (0) below recommended minimum (10-20).
#> Jackknife estimates may be unreliable (per papers S111, S114).
#> Consider aggregating samples or filtering genes with low abundance.
#> Warning: Gene unable to compute jackknife (diversity estimate is NA). Likely cause: zero or near-zero counts. Returning NA result.
#> Warning: Total count (0) below recommended minimum (10-20).
#> Jackknife estimates may be unreliable (per papers S111, S114).
#> Consider aggregating samples or filtering genes with low abundance.
#> Warning: Gene unable to compute jackknife (diversity estimate is NA). Likely cause: zero or near-zero counts. Returning NA result.
#> Warning: Total count (0) below recommended minimum (10-20).
#> Jackknife estimates may be unreliable (per papers S111, S114).
#> Consider aggregating samples or filtering genes with low abundance.
#> Warning: Gene unable to compute jackknife (diversity estimate is NA). Likely cause: zero or near-zero counts. Returning NA result.
#> Warning: Total count (0) below recommended minimum (10-20).
#> Jackknife estimates may be unreliable (per papers S111, S114).
#> Consider aggregating samples or filtering genes with low abundance.
#> Warning: Gene unable to compute jackknife (diversity estimate is NA). Likely cause: zero or near-zero counts. Returning NA result.
#> Warning: Total count (0) below recommended minimum (10-20).
#> Jackknife estimates may be unreliable (per papers S111, S114).
#> Consider aggregating samples or filtering genes with low abundance.
#> Warning: Gene unable to compute jackknife (diversity estimate is NA). Likely cause: zero or near-zero counts. Returning NA result.
#> Warning: Total count (8.946) below recommended minimum (10-20).
#> Jackknife estimates may be unreliable (per papers S111, S114).
#> Consider aggregating samples or filtering genes with low abundance.
#> Warning: Total count (0) below recommended minimum (10-20).
#> Jackknife estimates may be unreliable (per papers S111, S114).
#> Consider aggregating samples or filtering genes with low abundance.
#> Warning: Gene unable to compute jackknife (diversity estimate is NA). Likely cause: zero or near-zero counts. Returning NA result.
#> Warning: Total count (6.53) below recommended minimum (10-20).
#> Jackknife estimates may be unreliable (per papers S111, S114).
#> Consider aggregating samples or filtering genes with low abundance.
#> Warning: Total count (0) below recommended minimum (10-20).
#> Jackknife estimates may be unreliable (per papers S111, S114).
#> Consider aggregating samples or filtering genes with low abundance.
#> Warning: Gene unable to compute jackknife (diversity estimate is NA). Likely cause: zero or near-zero counts. Returning NA result.
#> Warning: Total count (0) below recommended minimum (10-20).
#> Jackknife estimates may be unreliable (per papers S111, S114).
#> Consider aggregating samples or filtering genes with low abundance.
#> Warning: Gene unable to compute jackknife (diversity estimate is NA). Likely cause: zero or near-zero counts. Returning NA result.
#> Warning: Total count (0) below recommended minimum (10-20).
#> Jackknife estimates may be unreliable (per papers S111, S114).
#> Consider aggregating samples or filtering genes with low abundance.
#> Warning: Gene unable to compute jackknife (diversity estimate is NA). Likely cause: zero or near-zero counts. Returning NA result.
#> Warning: Total count (0) below recommended minimum (10-20).
#> Jackknife estimates may be unreliable (per papers S111, S114).
#> Consider aggregating samples or filtering genes with low abundance.
#> Warning: Gene unable to compute jackknife (diversity estimate is NA). Likely cause: zero or near-zero counts. Returning NA result.
#> Warning: Total count (1) below recommended minimum (10-20).
#> Jackknife estimates may be unreliable (per papers S111, S114).
#> Consider aggregating samples or filtering genes with low abundance.
#> Warning: Total count (0) below recommended minimum (10-20).
#> Jackknife estimates may be unreliable (per papers S111, S114).
#> Consider aggregating samples or filtering genes with low abundance.
#> Warning: Gene unable to compute jackknife (diversity estimate is NA). Likely cause: zero or near-zero counts. Returning NA result.
#> Warning: Total count (0) below recommended minimum (10-20).
#> Jackknife estimates may be unreliable (per papers S111, S114).
#> Consider aggregating samples or filtering genes with low abundance.
#> Warning: Gene unable to compute jackknife (diversity estimate is NA). Likely cause: zero or near-zero counts. Returning NA result.
#> Warning: Total count (0) below recommended minimum (10-20).
#> Jackknife estimates may be unreliable (per papers S111, S114).
#> Consider aggregating samples or filtering genes with low abundance.
#> Warning: Gene unable to compute jackknife (diversity estimate is NA). Likely cause: zero or near-zero counts. Returning NA result.
#> Warning: Total count (0) below recommended minimum (10-20).
#> Jackknife estimates may be unreliable (per papers S111, S114).
#> Consider aggregating samples or filtering genes with low abundance.
#> Warning: Gene unable to compute jackknife (diversity estimate is NA). Likely cause: zero or near-zero counts. Returning NA result.
#> Warning: Total count (0) below recommended minimum (10-20).
#> Jackknife estimates may be unreliable (per papers S111, S114).
#> Consider aggregating samples or filtering genes with low abundance.
#> Warning: Gene unable to compute jackknife (diversity estimate is NA). Likely cause: zero or near-zero counts. Returning NA result.
#> Warning: Total count (0) below recommended minimum (10-20).
#> Jackknife estimates may be unreliable (per papers S111, S114).
#> Consider aggregating samples or filtering genes with low abundance.
#> Warning: Gene unable to compute jackknife (diversity estimate is NA). Likely cause: zero or near-zero counts. Returning NA result.
#> Warning: Total count (0) below recommended minimum (10-20).
#> Jackknife estimates may be unreliable (per papers S111, S114).
#> Consider aggregating samples or filtering genes with low abundance.
#> Warning: Gene unable to compute jackknife (diversity estimate is NA). Likely cause: zero or near-zero counts. Returning NA result.
#> Warning: Total count (8.946) below recommended minimum (10-20).
#> Jackknife estimates may be unreliable (per papers S111, S114).
#> Consider aggregating samples or filtering genes with low abundance.
#> Warning: Total count (0) below recommended minimum (10-20).
#> Jackknife estimates may be unreliable (per papers S111, S114).
#> Consider aggregating samples or filtering genes with low abundance.
#> Warning: Gene unable to compute jackknife (diversity estimate is NA). Likely cause: zero or near-zero counts. Returning NA result.
#> Warning: Total count (6.53) below recommended minimum (10-20).
#> Jackknife estimates may be unreliable (per papers S111, S114).
#> Consider aggregating samples or filtering genes with low abundance.
#> Warning: Total count (0) below recommended minimum (10-20).
#> Jackknife estimates may be unreliable (per papers S111, S114).
#> Consider aggregating samples or filtering genes with low abundance.
#> Warning: Gene unable to compute jackknife (diversity estimate is NA). Likely cause: zero or near-zero counts. Returning NA result.
#> Warning: Total count (0) below recommended minimum (10-20).
#> Jackknife estimates may be unreliable (per papers S111, S114).
#> Consider aggregating samples or filtering genes with low abundance.
#> Warning: Gene unable to compute jackknife (diversity estimate is NA). Likely cause: zero or near-zero counts. Returning NA result.
#> Warning: Total count (0) below recommended minimum (10-20).
#> Jackknife estimates may be unreliable (per papers S111, S114).
#> Consider aggregating samples or filtering genes with low abundance.
#> Warning: Gene unable to compute jackknife (diversity estimate is NA). Likely cause: zero or near-zero counts. Returning NA result.
#> Warning: Total count (0) below recommended minimum (10-20).
#> Jackknife estimates may be unreliable (per papers S111, S114).
#> Consider aggregating samples or filtering genes with low abundance.
#> Warning: Gene unable to compute jackknife (diversity estimate is NA). Likely cause: zero or near-zero counts. Returning NA result.
#> Warning: Total count (1) below recommended minimum (10-20).
#> Jackknife estimates may be unreliable (per papers S111, S114).
#> Consider aggregating samples or filtering genes with low abundance.
#> Warning: Total count (0) below recommended minimum (10-20).
#> Jackknife estimates may be unreliable (per papers S111, S114).
#> Consider aggregating samples or filtering genes with low abundance.
#> Warning: Gene unable to compute jackknife (diversity estimate is NA). Likely cause: zero or near-zero counts. Returning NA result.
#> Warning: Total count (0) below recommended minimum (10-20).
#> Jackknife estimates may be unreliable (per papers S111, S114).
#> Consider aggregating samples or filtering genes with low abundance.
#> Warning: Gene unable to compute jackknife (diversity estimate is NA). Likely cause: zero or near-zero counts. Returning NA result.
#> Warning: Total count (0) below recommended minimum (10-20).
#> Jackknife estimates may be unreliable (per papers S111, S114).
#> Consider aggregating samples or filtering genes with low abundance.
#> Warning: Gene unable to compute jackknife (diversity estimate is NA). Likely cause: zero or near-zero counts. Returning NA result.
#> Warning: Total count (0) below recommended minimum (10-20).
#> Jackknife estimates may be unreliable (per papers S111, S114).
#> Consider aggregating samples or filtering genes with low abundance.
#> Warning: Gene unable to compute jackknife (diversity estimate is NA). Likely cause: zero or near-zero counts. Returning NA result.
#> Warning: Total count (0) below recommended minimum (10-20).
#> Jackknife estimates may be unreliable (per papers S111, S114).
#> Consider aggregating samples or filtering genes with low abundance.
#> Warning: Gene unable to compute jackknife (diversity estimate is NA). Likely cause: zero or near-zero counts. Returning NA result.
#> Warning: Total count (0) below recommended minimum (10-20).
#> Jackknife estimates may be unreliable (per papers S111, S114).
#> Consider aggregating samples or filtering genes with low abundance.
#> Warning: Gene unable to compute jackknife (diversity estimate is NA). Likely cause: zero or near-zero counts. Returning NA result.
#> Warning: Total count (0) below recommended minimum (10-20).
#> Jackknife estimates may be unreliable (per papers S111, S114).
#> Consider aggregating samples or filtering genes with low abundance.
#> Warning: Gene unable to compute jackknife (diversity estimate is NA). Likely cause: zero or near-zero counts. Returning NA result.
#> Warning: Total count (8.946) below recommended minimum (10-20).
#> Jackknife estimates may be unreliable (per papers S111, S114).
#> Consider aggregating samples or filtering genes with low abundance.
#> Warning: Total count (0) below recommended minimum (10-20).
#> Jackknife estimates may be unreliable (per papers S111, S114).
#> Consider aggregating samples or filtering genes with low abundance.
#> Warning: Gene unable to compute jackknife (diversity estimate is NA). Likely cause: zero or near-zero counts. Returning NA result.
# Check jackknife results using unified accessor
jackknife_res <- results(analysis, type = 'jackknife')
if (!is.null(jackknife_res)) names(jackknife_res)
```
