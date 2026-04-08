# Calculate divergence metrics and store in TSENATAnalysis

Wrapper around .calculate_divergence() that manages TSENATAnalysis
object. Calculates Tsallis divergence between experimental conditions
for transcripts across multiple q-values to detect condition-specific
transcript remodeling.

## Usage

``` r
calculate_divergence_s4(
  analysis,
  q = NULL,
  verbose = FALSE,
  nthreads = NULL,
  output_file = NULL,
  control_group = NULL,
  paired = FALSE,
  method = NULL,
  bootstrap = FALSE,
  nboot = NULL,
  progress = FALSE,
  ...
)
```

## Arguments

- analysis:

  `TSENATAnalysis` object.

- q:

  `numeric`. Q-value for divergence. If NULL, uses first q_value from
  @config\$q_values if available, else defaults to 1.0.

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

\*\*Mathematical Background:\*\* Tsallis divergence D_q between two
probability distributions:


      D_q(P||Q) = (log_2(N) - entropy_q(P) + entropy_q(Q)) / (q - 1)

Measures how much transcript composition changes from control to
condition. Values near 0: Similar isoform composition; Large positive
values: Major change.

\*\*Example Use Case:\*\* Control sample: All reads from dominant
isoform (low entropy)  
Tumor sample: Reads spread across multiple isoforms (high entropy)  
Result: Large divergence indicating condition-specific isoform
switching.  

\*\*Parameter Resolution:\*\* Parameters are resolved in priority
order: 1. Explicit arguments passed to function 2. Values from
analysis@config (if present) 3. Function defaults

Requires diversity results from calculate_diversity_s4() as
prerequisite.

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
config <- tsenat_config(sample_col = 'sample', condition_col = 'condition')
analysis <- build_analysis_s4(
  readcounts = readcounts,
  tx2gene = gff3_dataset,
  metadata = metadata_df,
  config = config,
  tpm = tpm,
  effective_length = effective_length
)
# Use 200+ genes to ensure diversity filtering doesn't remove all genes
analysis <- filter_analysis_s4(analysis, min_samples = 1, subset_n_genes = 200)

# Compute diversity first (required for divergence)
analysis <- calculate_diversity_s4(analysis, q = c(0.5, 1.0, 1.5))
#> Note: 12 genes excluded (< 75% valid values).

# Calculate divergence across q-values
analysis <- calculate_divergence_s4(analysis, q = c(0.5, 1.0, 1.5))

# Check divergence results  
head(divergence(analysis))
#> $divergence_se
#> class: SummarizedExperiment 
#> dim: 86 3 
#> metadata(7): summary_stats elapsed_time_sec ... normalization
#>   computation_mode
#> assays(1): divergence
#> rownames(86): DES SYNM ... ENSG00000291072 GIPR
#> rowData names(26): gene_name error ... ci_width per_q_pattern
#> colnames(3): q_0.5 q_1 q_1.5
#> colData names(3): q_value sample_type computation_mode
#> 
```
