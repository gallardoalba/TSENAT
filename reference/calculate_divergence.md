# Calculate divergence metrics and store in TSENATAnalysis

Wrapper around .calculate_divergence() that manages TSENATAnalysis
object. Calculates Tsallis divergence between experimental conditions
for transcripts across multiple q-values to detect condition-specific
transcript remodeling.

## Usage

``` r
calculate_divergence(
  analysis,
  q = NULL,
  verbose = FALSE,
  nthreads = NULL,
  output_file = NULL,
  control_group = NULL,
  paired = NULL,
  method = NULL,
  bootstrap = NULL,
  nboot = NULL,
  progress = FALSE,
  ...
)
```

## Arguments

- analysis:

  `TSENATAnalysis` object.

- q:

  `numeric`. Q-value(s) for divergence (single value or vector). If
  NULL, reads from `analysis@config$q`. If not in config, defaults to
  seq(0.01, 2, by = 0.05) for full spectrum computation.

- verbose:

  `logical`. Print progress messages. Default: TRUE.

- nthreads:

  `numeric` or `NULL`. Number of CPU threads for parallel processing. If
  NULL, reads from `@config$nthreads` (or defaults to 1).

- output_file:

  `character` or `NULL`. Optional file path to save results. Supported
  formats: .rds (for S4 objects), .tsv, .csv, .txt (for tables).
  Default: NULL (no file output).

- control_group:

  `character` or `NULL`. Control group identifier for divergence
  comparison. If NULL, reads from `@config$control_group` if available.

- paired:

  `logical`. Whether to use paired design. Default: FALSE. If not
  specified, reads from `@config$paired` if available.

- method:

  `character` or `NULL`. Statistical method for divergence calculation.
  If NULL, reads from `@config$method` if available.

- bootstrap:

  `logical`. Whether to compute bootstrap confidence intervals. Default:
  FALSE. If not specified, reads from `@config$bootstrap` if available.

- nboot:

  `numeric` or `NULL`. Number of bootstrap replicates. Default: NULL. If
  NULL, reads from `@config$nboot` if available.

- progress:

  `logical`. Show progress bar during computation. Default: FALSE.

- ...:

  Additional arguments passed to the base divergence function.

## Value

Modified TSENATAnalysis with divergence metrics in `@divergence_results`
(stored as list of data.frames or matrices).

## Details

Key Features:

- Multi-q analysis: Divergence computed across full q-spectrum
  simultaneously

- Bootstrap confidence intervals: Quantify uncertainty in divergence
  estimates

- Multiple testing correction: Hochberg, Benjamini-Yekutieli, or no
  correction

- Paired designs: Supports paired/repeated measures via subject_col
  parameter

- Effect size reporting: Log-fold-change and confidence intervals per
  gene

- Flexible control group: Compare any condition vs. any other condition

\*\*Mathematical Background:\*\* Tsallis divergence between the
per-condition isoform-usage distributions \\P\\ (control) and \\Q\\
(treatment), built by summing transcript counts across samples within
each condition:


      D_q(P||Q) = (sum_i P_i^q * Q_i^(1-q) - 1) / (q - 1)   for q > 0, q != 1
      D_1(P||Q) = sum_i P_i * log(P_i / Q_i)                 (KL limit)

q = 0 convention: D_0(P\|\|Q) = 1 - sum_i: P_i \> 0 Q_i (the Q-mass on
P's zero support, with 0^0 = 0). Under pseudocount regularization all
bins are positive and D_0 = 0. Measures how much isoform composition
changes from control to condition. Values near 0: Similar isoform
composition; Large positive values: Major change. Note: the default is
norm = 'none', so reported values are raw divergences on the count scale
(not rescaled across genes).

\*\*Example Use Case:\*\* Control sample: All reads from dominant
isoform (low entropy)  
Tumor sample: Reads spread across multiple isoforms (high entropy)  
Result: Large divergence indicating condition-specific isoform
switching.  

\*\*Parameter Resolution:\*\* Parameters are resolved in priority
order: 1. Explicit arguments passed to function 2. Values from
analysis@config (if present) 3. Function defaults

Requires diversity results from calculate_diversity() as prerequisite.

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
gff3_dataset <- system.file('extdata', 'annotation.gff3.gz', package = 'TSENAT')

# Build analysis from vignette data and create manageable subset
config <- TSENAT_config(sample_col = 'sample', condition_col = 'condition')
analysis <- build_analysis(
  readcounts = readcounts,
  tx2gene = gff3_dataset,
  metadata = metadata_df,
  config = config,
  tpm = tpm,
  effective_length = effective_length
)
# Use 200+ genes to ensure diversity filtering doesn't remove all genes
analysis <- filter_analysis(analysis, min_samples = 1, subset_n_genes = 200)

# Compute diversity first (required for divergence)
analysis <- calculate_diversity(analysis, q = c(0.5, 1.0, 1.5))

# Calculate divergence across q-values
analysis <- calculate_divergence(analysis, q = c(0.5, 1.0, 1.5))

# Check divergence results using unified accessor
head(results(analysis, type = 'divergence'))
#>                 q_0.5          q_1        q_1.5
#> DES      0.000000e+00 0.000000e+00 0.000000e+00
#> SYNM     9.296544e-06 1.856638e-05 2.780964e-05
#> HIPK2    0.000000e+00 0.000000e+00 0.000000e+00
#> REG3A    2.960352e-03 5.910282e-03 8.854202e-03
#> SH3PXD2A 1.492930e-04 2.919431e-04 4.284469e-04
#> EFNB2    0.000000e+00 0.000000e+00 0.000000e+00
```
