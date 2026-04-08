# Prepare Gene Switching Tables from TSENATAnalysis Object

S4 wrapper for `.prepare_gene_switching_tables()` that extracts results
directly from a TSENATAnalysis object. Automatically retrieves LM
results and jackknife switching results from the analysis object slots.

## Usage

``` r
prepare_gene_switching_tables_s4(
  analysis,
  n_top_genes = NULL,
  n_transcripts_per_gene = 10,
  verbose = FALSE,
  output_file = NULL,
  ...
)
```

## Arguments

- analysis:

  `TSENATAnalysis`. An S4 object containing completed LM interaction and
  jackknife isoform switching analyses.

- n_top_genes:

  `numeric` or `NULL`. Number of top genes (by adjusted p-value) to
  include in output tables. If `NULL`, all genes with significant LM
  results are included.

- n_transcripts_per_gene:

  `numeric`. Maximum number of transcripts to display per gene in the
  output tables (default: 10).

- verbose:

  `logical`. If `TRUE`, print diagnostic messages during table
  preparation.

- output_file:

  `character` or `NULL`. Optional file path to save results. Supported
  formats: .rds (for S4 objects), .tsv, .csv, .txt (for tables).
  Default: NULL (no file output).

- ...:

  Additional arguments passed to the base function.

## Value

A list containing:

- `$summary_table`:

  Gene-level summary with LM p-values and significant q-values

- `$transcript_tables`:

  Named list of data.frames, one per gene, showing transcript-level
  switching metrics

- `$q_vector`:

  Vector of q-values analyzed

## Details

This function extracts the following from `analysis`:

- LM results:

  From `analysis@lm_results$lm_interaction$results`

- Jackknife results:

  From `analysis@jackknife_results` or extracted from the switching
  analysis metadata

The wrapper automatically handles column detection and parameter
extraction, providing a simplified interface compared to the base
function.

## Examples

``` r
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
  subject_col = 'paired_samples',
  q_values = seq(0, 2, by = 0.1)
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
analysis <- calculate_diversity_s4(analysis, q = c(0.5, 1.0, 1.5, 2.0, 2.5))
#> Note: 7 genes excluded (< 75% valid values).
analysis <- calculate_lm_interaction_s4(analysis, method = 'gam')
#> Warning: nlminb problem, convergence error code = 1
#>   message = singular convergence (7)
#> Warning: nlminb problem, convergence error code = 1
#>   message = iteration limit reached without convergence (10)
#> Warning: nlminb problem, convergence error code = 1
#>   message = iteration limit reached without convergence (10)
analysis <- jackknife_isoform_switching_s4(analysis, n_bootstrap = 50)
tables <- prepare_gene_switching_tables_s4(analysis)
head(tables$summary_df)
#>       gene gene_name p_interaction adj_p_interaction
#> 1   ZNF493    ZNF493  3.075886e-07      9.842836e-06
#> 2    CENPV     CENPV  5.046008e-06      1.564262e-04
#> 3     VRK2      VRK2  1.229501e-05      3.688503e-04
#> 4     TLN2      TLN2  7.474492e-05      2.167603e-03
#> 5 TMEM183A  TMEM183A  2.528512e-04      7.079833e-03
#> 6   CXCL12    CXCL12  2.425630e-03      6.549200e-02
```
