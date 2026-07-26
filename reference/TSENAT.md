# Run complete TSENAT analysis pipeline

Coordinates the full TSENAT workflow: diversity -\> jackknife -\> LM
interactions -\> divergence -\> gene interactions -\> rank-based tests
-\> concordance -\> visualizations.

## Usage

``` r
TSENAT(
  analysis,
  output_dir = "tsenat_outputs",
  save_output = TRUE,
  output_format = c("tsv", "csv", "txt", "rds"),
  verbose = TRUE,
  tips = FALSE
)
```

## Arguments

- analysis:

  `TSENATAnalysis` object created by
  [`build_analysis`](https://gallardoalba.github.io/TSENAT/reference/build_analysis.md).

- output_dir:

  `character`. Directory to save results and plots. Default:
  'tsenat_outputs'. Set to NULL to disable automatic output saving.

- save_output:

  `logical`. Whether to save output files (results tables). Default:
  TRUE. If FALSE, no TSV/CSV output files are written to disk.

- output_format:

  `character`. Format for output files: 'tsv' (tab-separated), 'csv'
  (comma-separated), 'txt' (text), or 'rds' (R serialized). Default:
  'tsv'.

- verbose:

  `logical`. Print progress messages. Default: TRUE.

- tips:

  `logical`. Show tips and usage examples at the end of the pipeline
  run. Default: FALSE. Set to TRUE to see the \`\[TIPS\]\` section
  displaying common result-extraction commands.

## Value

`TSENATAnalysis` object containing complete analysis results, plots, and
metadata.

## Details

Pipeline execution order (enforced, follows TSENAT.Rmd vignette):

1.  [`filter_analysis()`](https://gallardoalba.github.io/TSENAT/reference/filter_analysis.md) -
    Filter low-abundance transcripts

2.  [`calculate_diversity()`](https://gallardoalba.github.io/TSENAT/reference/calculate_diversity.md) -
    Tsallis entropy per q-value

3.  [`plot_diversity_spectrum()`](https://gallardoalba.github.io/TSENAT/reference/plot_diversity_spectrum.md) -
    Visualize q-spectrum

4.  [`calculate_m_estimator()`](https://gallardoalba.github.io/TSENAT/reference/calculate_m_estimator.md) -
    Sample influence QC analysis

5.  [`calculate_sait()`](https://gallardoalba.github.io/TSENAT/reference/calculate_sait.md) -
    LM/SAIT interaction testing

6.  [`plot_sait()`](https://gallardoalba.github.io/TSENAT/reference/plot_sait.md) -
    SAIT results visualization

7.  [`calculate_jis()`](https://gallardoalba.github.io/TSENAT/reference/calculate_jis.md) -
    Transcript switching detection

8.  [`plot_jis_delta()`](https://gallardoalba.github.io/TSENAT/reference/plot_jis_delta.md) -
    Multi-q influence heatmap (gene switching tables computed lazily via
    results())

9.  [`plot_expression()`](https://gallardoalba.github.io/TSENAT/reference/plot_expression.md) -
    Top transcript visualization

10. [`calculate_divergence()`](https://gallardoalba.github.io/TSENAT/reference/calculate_divergence.md) -
    Pairwise divergence metrics

11. [`calculate_effect_sizes()`](https://gallardoalba.github.io/TSENAT/reference/calculate_effect_sizes.md) -
    Effect size computation

12. [`plot_divergence_distribution()`](https://gallardoalba.github.io/TSENAT/reference/plot_divergence_distribution.md) -
    Divergence distribution plot

13. [`plot_divergence_spectrum()`](https://gallardoalba.github.io/TSENAT/reference/plot_divergence_spectrum.md) -
    Divergence spectrum plot

14. [`calculate_assumptions()`](https://gallardoalba.github.io/TSENAT/reference/calculate_assumptions.md) -
    Validate rank-based test assumptions

15. [`calculate_rank_transform()`](https://gallardoalba.github.io/TSENAT/reference/calculate_rank_transform.md) -
    ART (Aligned Rank Transform) interaction test

16. [`calculate_concordance()`](https://gallardoalba.github.io/TSENAT/reference/calculate_concordance.md) -
    Compare LM and rank test results

## Examples

``` r
# \donttest{
data(readcounts, package = 'TSENAT')
metadata_df <- read.table(
  system.file('extdata', 'metadata.tsv', package = 'TSENAT'),
  header = TRUE, sep = '\t'
)
gff3_file <- system.file('extdata', 'annotation.gff3.gz', package = 'TSENAT')

config <- TSENAT_config(
  sample_col = 'sample',
  condition_col = 'condition',
  q = seq(0, 2, length.out = 10),
  generate_plots = FALSE
)
analysis <- build_analysis(
  readcounts = as.matrix(readcounts),
  tx2gene = gff3_file,
  metadata = metadata_df,
  config = config,
  tpm = tpm,
  effective_length = effective_length
)

result <- TSENAT(analysis)
#> Created output directory: tsenat_outputs
#> 
#> +============================================================+
#> |          TSENAT: Tsallis Entropy Analysis Toolbox          |
#> +============================================================+
#>                                                               
#>       Science is an essentially anarchic enterprise.          
#>                                                               
#>                        -- Paul Feyerabend, Against Method     
#>                                                               
#> [DATA] Data Summary
#>   Transcripts .......... 3,089
#>   Samples .............. 16
#>   Conditions ........... 2
#>   Q-spectrum range ..... 0 to 2 (10 values)
#> 
#> [CONFIG] Analysis Configuration
#>   Design ............... unpaired
#>   Filter stringency .... medium
#>   Normalization ........ enabled [0-1]
#>   Normalization method . range (default)
#>   Pseudocount .......... disabled
#>   Shrinkage ............ disabled
#>   Significance ......... p < 0.050 | FDR < 0.050
#>   SAIT method .......... GAM
#>   SAIT p-corr method ... BH
#>   Jackknife use_sait_fdr . TRUE
#> 
#> =============================================================
#> [>] [ 1/16] Filtering low-abundance transcripts
#>           [OK] Complete - 341 transcripts remaining
#> [>] [ 2/16] Computing Tsallis entropy
#>           [OK] 10 q-values processed
#> [>] [ 3/16] Plotting diversity q-spectrum
#>           [OK] Plot generated
#> [>] [ 4/16] Computing M-estimator influence
#>           [OK] M-estimate QC complete
#> [>] [ 5/16] Fitting SAIT interaction models
#>           [OK] SAIT interaction analysis complete
#> [>] [ 6/16] Plotting SAIT results
#>           [OK] SAIT interaction plot generated
#> [>] [ 7/16] Computing jackknife isoform switching
#>           [OK] Jackknife isoform switching complete
#> [>] [ 8/16] Plotting influence heatmap
#>           [OK] Influence heatmap generated
#> [>] [ 9/16] Plotting top transcripts
#>           [OK] Top transcripts plot generated
#> [>] [10/16] Computing divergence metrics
#>           [OK] Divergence computed
#> [>] [11/16] Computing effect sizes
#>           [OK] Effect sizes computed
#> [>] [12/16] Plotting divergence distributions
#>           [OK] Divergence distribution plot generated
#> [>] [13/16] Plotting divergence spectrum
#>           [OK] Global divergence spectrum plot generated
#>           [OK] Multi-gene divergence spectrum plot generated
#> [>] [14/16] Checking statistical assumptions
#>           [OK] Assumptions validated
#> [>] [15/16] Performing ART (Aligned Rank Transform) interaction test
#>           [OK] Aligned Rank Transform (ART) completed
#> [>] [16/16] Computing SAIT-ART concordance
#>           [OK] Concordance analysis completed
#> =============================================================
#> 
#> +============================================================+
#> |               [OK] ANALYSIS COMPLETE                       |
#> +============================================================+
#> 
#> [RESULTS] Results Summary
#> 
#> [PERF] Performance
#>   Total time ........... 32.3s
#>   Slowest steps:
#>     1. sait_interaction     14.3s (44.4%)
#>     2. rank_transform_test  5.8s (18.1%)
#>     3. jackknife            4.3s (13.2%)
#> 
#> [OUTPUT] Output
#>   Directory ........... tsenat_outputs
#>   Files saved ......... 18
# }
```
