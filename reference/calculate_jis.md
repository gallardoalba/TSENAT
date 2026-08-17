# Jackknife isoform switching analysis on TSENATAnalysis object

Wrapper around .calculate_jis() that manages TSENATAnalysis object.
Identifies transcripts with significant isoform switching patterns using
jackknife resampling across samples to detect influential isoforms.

## Usage

``` r
calculate_jis(
  analysis,
  condition_col = NULL,
  subject_col = NULL,
  gene_col = NULL,
  isoform_col = NULL,
  q = c(0, 0.5, 1, 1.5, 2),
  norm = NULL,
  log_base = NULL,
  threshold = 90,
  nboot = 1000,
  pseudocount = NULL,
  sait_results = NULL,
  sait_p_threshold = 0.05,
  use_sait_fdr = TRUE,
  nthreads = NULL,
  output_file = NULL,
  verbose = FALSE,
  ...
)
```

## Arguments

- analysis:

  `TSENATAnalysis` object containing:

  - `@se`: SummarizedExperiment with count data

  - `@config`: Configuration metadata

- condition_col:

  `character`. Column name in colData(se) specifying group assignments
  (default: 'sample_type'). If NULL, attempts auto-detection.

- subject_col:

  `character`. Optional column for paired/repeated measures design. If
  provided, enables paired analysis. Default: NULL (unpaired).

- gene_col:

  `character`. Column name in rowData(se) or metadata identifying genes.
  Default: 'gene'.

- isoform_col:

  `character`. Column name in rowData(se) or metadata identifying
  isoforms/transcripts. Default: 'transcript' or 'isoform'.

- q:

  `numeric`. Tsallis entropy parameter(s) to analyze. Can be single
  value or vector for multi-q analysis (default: c(0, 0.5, 1, 1.5, 2)).
  If NULL, uses @config\$q.

- norm:

  `logical`. Whether to use normalized diversity values (default: TRUE).

- log_base:

  `numeric`. Logarithm base for entropy calculations. Default: NULL
  (uses e).

- threshold:

  `numeric`. Percentile threshold for detecting transcript switching
  (default: 90). Transcripts with delta_influence \>= threshold
  percentile are classified as 'switching'.

- nboot:

  `integer`. Number of bootstrap resamples for confidence intervals
  (default: 1000).

- pseudocount:

  `numeric`. Pseudocount value for count regularization. Default: NULL
  (uses @config\$pseudocount or 0).

- sait_results:

  `data. frame`. Optional SAIT interaction results to filter genes. If
  provided, only genes in sait_results are analyzed.

- sait_p_threshold:

  `numeric`. P-value threshold for filtering genes from sait_results
  (default: 0.05).

- use_sait_fdr:

  `logical`. If TRUE, uses adjusted p-values from sait_results (default:
  TRUE).

- nthreads:

  `numeric` or `NULL`. Number of CPU threads for parallel processing. If
  NULL, reads from `@config$nthreads` (or defaults to 1). If \> 1 and
  multiple q-values provided, uses parallel PSOCK cluster.

- output_file:

  `character` or `NULL`. Optional file path to save results. Supported
  formats: .tsv, .csv, .txt (for tables), or .rds (for S4 objects). When
  text format is specified, generates TWO files:

  - `output_file`: Gene-level summary (one row per gene with switching
    statistics)

  - `output_file_transcripts. ext`: Transcript-level details (one row
    per transcript with p-values and FDR)

  For .rds format, saves only the full analysis object. Default: NULL
  (no file output).

- verbose:

  `logical`. Print progress messages (default: FALSE).

- ...:

  Additional arguments for future extensibility.

## Value

`TSENATAnalysis` object with jackknife results stored in
`@jackknife_results` slot. Results are keyed by q-value (e.g.,
'q_1.00'). For multi-q analysis, multiple calls will accumulate results
in the slot.

The analysis object is returned visibly to support method chaining:


        analysis <- calculate_jis(analysis, q = 0.5)
        analysis <- calculate_jis(analysis, q = 1.0)
      

## Details

Key Features:

- Jackknife resampling: Robust outlier detection across all samples

- Delta-influence metric: Measures how much each isoform drives
  phenotype

- Confidence intervals: Bootstrap-based uncertainty quantification

- Multi-q analysis: Tests across full q-spectrum simultaneously

- LM filtering: Optional restriction to genes with significant
  interactions

- Paired designs: Supports repeated measures/longitudinal data

\*\*Parameter Resolution from Config\*\*

The following parameters are resolved using a three-level priority
system:

1.  User-provided argument (if not NULL)

2.  Value from `analysis@config` (if key exists)

3.  Function default value

Affected parameters:

- `q`: Multi-q vector c(0, 0.5, 1, 1.5, 2) if not provided, or
  `@config$q` if available

- `nboot`: Uses `@config$nboot` if available, else 1000

- `threshold`: Uses `@config$threshold` if available, else 90

- `sait_p_threshold`: Uses `@config$sait_p_threshold` if available, else
  0.05

This allows setting defaults once in the config and reusing across
multiple analyses.

\*\*Mathematical Background:\*\* Delta-influence measures how much
removing each sample changes entropy:


      Delta = H_q(leave-one-out) - H_q(original)

High \|Delta\| for specific isoforms indicates those isoforms drive
differences. Identifies 'outlier samples' where isoforms contribute
unusually much.

\*\*Example Use Case:\*\* Sample shows high Delta for isoform X → X has
outsized importance in that sample  
Classifying transcripts as 'switching' if top percentile (e.g., 90th)
Delta  
Reveals condition-specific isoforms crucial for phenotype
determination.  

\*\*Automatic Setup:\*\* This wrapper automatically: 1. Extracts
SummarizedExperiment from `@se` slot 2. Detects condition_col, gene_col,
isoform_col from colData/rowData or `@config` 3. Calls
`.calculate_jis()` with extracted parameters

\*\*Parameter Auto-Detection:\*\*

1.  `condition_col`: Uses explicit parameter, then @config, then
    'sample_type'

2.  `gene_col`: Uses explicit parameter, then looks for 'gene' or 'Gene'

3.  `isoform_col`: Uses explicit parameter, then looks for 'transcript',
    'isoform', or 'Isoform'

## See also

`jackknife_isoform_switching` for the underlying implementation,
[`TSENATAnalysis`](https://gallardoalba.github.io/TSENAT/reference/TSENATAnalysis.md)
for object structure.

## Examples

``` r
data(readcounts)
metadata_df <- read.table(system.file('extdata', 'metadata.tsv', package
= 'TSENAT'),
                          header = TRUE, sep = '\t')
gff3_file <- system.file('extdata', 'annotation.gff3.gz', package = 'TSENAT')
config <- TSENAT_config(sample_col = 'sample', condition_col = 'condition')
analysis <- build_analysis(readcounts = readcounts, tx2gene =
gff3_file, metadata = metadata_df, config = config,
                             tpm = tpm,
 effective_length = effective_length)
analysis <- filter_analysis(analysis, min_samples = 1, subset_n_genes
= 20, subset_n_samples = 8)
analysis <- calculate_diversity(analysis, q = 1)
```
