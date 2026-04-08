# Calculate diversity and store in TSENATAnalysis

Wrapper around .calculate_diversity() that manages TSENATAnalysis
object. Calculates Tsallis entropy (diversity) across multiple q-values
for each gene to quantify isoform complexity and transcript
heterogeneity.

## Usage

``` r
calculate_diversity_s4(
  analysis,
  q = NULL,
  norm = TRUE,
  norm_method = NULL,
  reference_group = NULL,
  tpm = FALSE,
  assayno = NULL,
  verbose = NULL,
  show_messages = FALSE,
  what = NULL,
  nthreads = NULL,
  pseudocount = NULL,
  min_valid_frac = NULL,
  shrinkage = NULL,
  genes = NULL,
  effective_length = NULL,
  metadata = NULL,
  bootstrap = NULL,
  nboot = NULL,
  bootstrap_method = NULL,
  bootstrap_ci = NULL,
  bootstrap_include_diagnostics = NULL,
  output_file = NULL,
  ...
)
```

## Arguments

- analysis:

  `TSENATAnalysis` object.

- q:

  `numeric`. Q-value(s) for Tsallis entropy. If NULL, uses q_values from
  `analysis@config$q_values` if available, else defaults to seq(0. 01,
  2, by = 0. 05).

- norm:

  `logical` or `character`. Normalization method: TRUE, FALSE, 'none',
  'range', 'zscore', 'log_odds_ratio', 'relative_reference'. If NULL,
  reads from `@config$norm` or defaults to TRUE.

- norm_method:

  `character`. Post-hoc normalization method applied after diversity
  computation. Options:

  - `'default'` - Simple normalization by theoretical maximum (current
    behavior)

  - `'zscore'` - Z-score normalization per q-value

  - `'log_odds_ratio'` - Log-odds ratio relative to max entropy (q and
    isoform-aware)

  - `'relative_reference'` - Divide by reference group mean (requires
    reference_group)

  - `NULL` - No post-hoc normalization (default)

  If NULL, reads from `@config$norm_method` if available.

- reference_group:

  `character`. For `norm_method = 'relative_reference'`, the reference
  group column name (e.g., from colData). If NULL, uses first group in
  colData. If NULL, reads from `@config$reference_group` if available.

- tpm:

  `logical`. TPM normalization. Default: FALSE. If not specified, reads
  from `@config$tpm` if available.

- assayno:

  `numeric`. Assay number to use. Default: 1. If NULL, reads from
  `@config$assayno` if available.

- verbose:

  `logical`. Print progress messages. Default: TRUE. If not specified,
  reads from `@config$verbose` if available.

- show_messages:

  `logical`. Display verbose messages during computation. Default:
  FALSE.

- what:

  `character`. Output type: 'S' (entropy) or 'D' (diversity). Default:
  'S'. If NULL, reads from `@config$what` if available.

- nthreads:

  `numeric` or `NULL`. Number of CPU threads for parallel processing. If
  NULL, reads from `@config$nthreads` (or defaults to 1).

- pseudocount:

  `numeric` or `character`. Pseudocount value or 'auto'. Default: 0. If
  NULL, reads from `@config$pseudocount` if available.

- min_valid_frac:

  `numeric`. Minimum fraction of valid samples for gene filtering.
  Default: NULL. If NULL, reads from `@config$min_valid_frac` if
  available.

- shrinkage:

  `character`. Shrinkage method: 'none' or 'empirical_bayes'. Default:
  'none'. If NULL, reads from `@config$shrinkage` if available.

- genes:

  `character` or `NULL`. Gene set specification. Default: NULL (use all
  genes). If NULL, reads from `@config$genes` if available.

- effective_length:

  `numeric` or `NULL`. Effective gene lengths. Default: NULL. If NULL,
  reads from `@config$effective_length` if available.

- metadata:

  `list` or `NULL`. Additional metadata. Default: NULL.

- bootstrap:

  `logical`. Compute bootstrap confidence intervals. Default: FALSE. If
  not specified, reads from `@config$bootstrap` if available.

- nboot:

  `numeric` or `NULL`. Number of bootstrap replicates. Default: NULL. If
  NULL, reads from `@config$nboot` if available.

- bootstrap_method:

  `character`. Bootstrap method: 'percentile' or others. Default:
  'percentile'. If NULL, reads from `@config$bootstrap_method` if
  available.

- bootstrap_ci:

  `numeric`. Bootstrap confidence interval level (0-1). Default: 0. 95.
  If NULL, reads from `@config$bootstrap_ci` if available.

- bootstrap_include_diagnostics:

  `logical`. Include bootstrap diagnostic information. Default: FALSE.
  If NULL, reads from `@config$bootstrap_include_diagnostics` if
  available.

- output_file:

  `character` or `NULL`. Optional file path to save results. When
  provided, generates TWO files:

  1.  Primary output: Analysis object (.rds) or table (.tsv/.csv/.txt)

  2.  Secondary output: Diversity spectrum statistics (TSV format) with
      suffix `_spectrum.tsv`

  Example: output_file = 'analysis.rds' generates:

  - `analysis.rds` - TSENATAnalysis object

  - `analysis_spectrum.tsv` - Spectrum statistics

  The spectrum file contains columns: q, central (median diversity),
  spread (IQR), count, and group (if grouping variable available).
  Default: NULL (no file output).

- ...:

  Additional arguments passed to underlying functions for extensibility.

## Value

Modified TSENATAnalysis object with diversity results stored in
`@diversity_results`, keyed by 'q_X. XX. . . ' format (e. g. , 'q_1.
000'). When `output_file` is provided, also generates:

- Primary file: Analysis object or table export

- Spectrum file: `*_spectrum.tsv` containing aggregated diversity
  statistics across q-values and groups

## Details

Key Features:

- Multi-q analysis: Entropy computed for q = 0.01 to 2.00 (by default)

- Normalized entropy: Scale to \[0, 1\] using theoretical maximum log(m)

- Bootstrap confidence intervals: Quantify uncertainty in estimates

- TPM normalization: Optional SALMON TPM-based weighting

- Multiple normalization methods: Range, Z-score, log-odds-ratio,
  relative

- Spectrum computation: Aggregate diversity statistics across all
  q-values

\*\*Mathematical Background:\*\* Tsallis entropy H_q for q-parameter:


      H_q(X) = (1/(q-1)) * (1 - Sum p_i^q)  [for q != 1]
      H_1(X) = -Sum p_i * log(p_i)  [Shannon entropy, limit q→1]

where p_i = relative abundance of isoform i for a gene. Larger q
emphasizes dominant isoforms; smaller q emphasizes rare ones.

\*\*Example Interpretation:\*\*

- Gene with 1 isoform: H_q = 0 for all q (no diversity)

- Gene with 2 equal isoforms: H_q ~ 0.5-1.0 depending on q

- Gene with m equally abundant isoforms: H_q = 1.0 (maximum diversity)

This wrapper calls `.calculate_diversity()` once per q-value, storing
results as SummarizedExperiment objects. It extracts key parameters from
`analysis@config` with priority resolution (explicit \> `@config` \>
default).

\*\*Diversity Spectrum Computation:\*\* By default, this function
computes and saves a diversity spectrum (aggregated statistics across
all q-values and groups) when `output_file` is provided. The spectrum
contains:

- `q`: Diversity parameter value

- `central`: Median (or mean) diversity across all genes

- `spread`: IQR (or SD) around central value

- `count`: Number of valid measurements

- `group`: Condition group (if applicable)

This provides a quick summary of how diversity changes across q-values,
useful for q-curve visualization and statistical comparisons. Spectrum
is saved as: `*_spectrum.tsv`

\*\*Parameter Priority Resolution:\*\*

- q:

  Priority 1 (explicit) \> Priority 2 (`@config$q_values`) \> Priority 3
  (default: seq(0. 01, 2, by = 0. 05))  
  \*\*Note: \*\* If explicit q AND `@config$q_values` both provided,
  explicit wins.

- nthreads:

  Priority: explicit \> `@config$nthreads` \> 1

- verbose:

  Priority: explicit \> `@config$verbose` \> TRUE

- bootstrap:

  Priority: explicit \> `@config$bootstrap` \> FALSE

- pseudocount:

  Priority: explicit \> `@config$pseudocount` \> 0

- norm:

  Priority: explicit \> `@config$norm` \> TRUE

- what:

  Priority: explicit \> `@config$what` \> 'S' (Tsallis entropy)

\*\*Audit Trail:\*\* After execution, check:

- `analysis@config$last_diversity_run$parameters_used`: Actual
  parameters (not original `@config`)

- `attr(diversity(analysis, q=X), 'computed_with')`: Per-q-value
  metadata (timestamp, bootstrap setting, nthreads, etc.)

For additional details on diversity spectrum calculations and
normalization methods, see the package vignettes.

## Examples

``` r
# Load vignette data and build analysis
data(readcounts)
metadata_df <- read.table(
  system.file('extdata', 'metadata.tsv', package = 'TSENAT'),
  header = TRUE, sep = '\t'
)
gff3_dataset <- system.file('extdata', 'annotation.gff3.gz', package =
'TSENAT')
readcounts <- as.matrix(readcounts)
mode(readcounts) <- 'numeric'

config <- tsenat_config(sample_col = 'sample', condition_col = 'condition')
analysis <- build_analysis_s4(readcounts = readcounts, tx2gene =
gff3_dataset, metadata = metadata_df, config = config,
  tpm = tpm, effective_length = effective_length)

# Filter to manageable size (use 200+ genes to survive diversity filtering)
analysis <- filter_analysis_s4(analysis, min_samples = 1, subset_n_genes
= 200)

# Compute diversity and access results
analysis <- calculate_diversity_s4(analysis, q = c(0.5, 1.0), verbose =
FALSE)
head(diversity(analysis, q = 1.0))
#> class: SummarizedExperiment 
#> dim: 6 16 
#> metadata(11): q what ... sample_base_names samples
#> assays(2): diversity counts
#> rownames(6): ENSG00000065970.10 ENSG00000072954.8 ...
#>   ENSG00000110429.15 ENSG00000122912.16
#> rowData names(1): gene_id
#> colnames(16): SRR14800481 SRR14800480 ... SRR14800483 SRR14800482
#> colData names(5): condition sample_type sample_base paired_samples
#>   sample_id
```
