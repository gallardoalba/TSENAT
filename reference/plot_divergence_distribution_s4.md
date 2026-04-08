# Plot Tsallis Divergence Effect Size Distribution (S4 Wrapper)

S4 wrapper that extracts effect size results from a TSENATAnalysis
object and generates a histogram visualization of Tsallis divergence
effect sizes across genes.

## Usage

``` r
plot_divergence_distribution_s4(
  analysis,
  threshold = 0.1,
  output_file = NULL,
  width = 12,
  height = 6,
  verbose = FALSE,
  ...
)
```

## Arguments

- analysis:

  `TSENATAnalysis` object with effect sizes computed (typically via
  [`effect_sizes_divergence_s4`](https://gallardoalba.github.io/TSENAT/reference/effect_sizes_divergence_s4.md)).

- threshold:

  `numeric`. Effect size threshold for visual marking in the plot.
  Default is 0.1 (information-theoretic significance level).

- output_file:

  `character`. Optional file path to save the plot. If NULL, plot is
  returned but not saved.

- width:

  `numeric`. Plot width in inches. Default is 10.

- height:

  `numeric`. Plot height in inches. Default is 6.

- verbose:

  `logical`. Print status messages. Default is TRUE.

- ...:

  Additional arguments passed to the underlying plotting function.

## Value

Invisibly returns the file path if saved, otherwise the ggplot object.
If the plot cannot be created (missing data, ggplot2 not available),
returns NULL invisibly with an informative message.

## Details

This wrapper extracts the interaction results (with effect size columns)
from `analysis@metadata$effect_sizes_divergence$interaction_results` and
passes them to the base `.plot_divergence_distribution()` function.

The function visualizes the distribution of effect sizes using the
median q-value's divergence (typically around q=1.0, close to Shannon
entropy). A red dashed line marks the information-theoretic significance
threshold.

\*\*Data Requirements:\*\*

- Effect sizes must be computed via
  [`effect_sizes_divergence_s4()`](https://gallardoalba.github.io/TSENAT/reference/effect_sizes_divergence_s4.md)

- `@metadata$effect_sizes_divergence$interaction_results` must contain
  columns matching pattern `effect_size_D_q*`

## See also

[`effect_sizes_divergence_s4`](https://gallardoalba.github.io/TSENAT/reference/effect_sizes_divergence_s4.md)
for computing effect sizes.

## Examples

``` r
# Plot 2: Distribution of effect sizes across genes
data(readcounts)
readcounts <- as.matrix(readcounts)
mode(readcounts) <- 'numeric'
metadata_df <- read.table(
  system.file('extdata', 'metadata.tsv', package = 'TSENAT'),
  header = TRUE, sep = '\t'
)
gff3_dataset <- system.file('extdata', 'annotation.gff3.gz', package =
'TSENAT')
# Configure analysis parameters first
config <- tsenat_config(
  sample_col = 'sample',
  condition_col = 'condition',
  control = 'normal'
)

# Build analysis with configured parameters
analysis <- build_analysis_s4(
  readcounts = readcounts,
  tx2gene = gff3_dataset,
  metadata = metadata_df,
  config = config,
  tpm = tpm,
  effective_length = effective_length
)

analysis <- filter_analysis_s4(analysis, stringency = 'severe')
analysis <- calculate_diversity_s4(analysis, q = c(0.5, 1.0, 1.5, 2.0, 2.5), verbose
= FALSE)
analysis <- calculate_divergence_s4(analysis, q = c(0.5, 1.0, 1.5, 2.0, 2.5), verbose
= FALSE)
analysis <- calculate_lm_interaction_s4(analysis, method = 'gam')
analysis <- effect_sizes_divergence_s4(analysis)
p_dist <- plot_divergence_distribution_s4(analysis)
print(p_dist)

```
