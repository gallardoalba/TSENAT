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
  output_format = "tsv",
  verbose = TRUE
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

5.  [`calculate_lm()`](https://gallardoalba.github.io/TSENAT/reference/calculate_lm.md) -
    LM interaction testing

6.  [`plot_lm()`](https://gallardoalba.github.io/TSENAT/reference/plot_lm.md) -
    LM results visualization

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

15. [`calculate_srh()`](https://gallardoalba.github.io/TSENAT/reference/calculate_srh.md) -
    Scheirer-Ray-Hare rank-based interaction test

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
#>   LM method ............ GAM
#>   LM p-corr method ..... BH
#>   Jackknife use_lm_fdr . TRUE
#> 
#> =============================================================
#> [>] [%2d/16] Filtering low-abundance transcripts
#>           [OK] Complete
#> [>] [%2d/16] Computing Tsallis entropy2
#> [>] [ 2/14] Computing Tsallis diversity
#>           [OK] 10 q-values processed
#> [>] [%2d/16] Plotting diversity q-spectrum3
#> [>] [ 3/14] Plotting q-spectrum curve
#>           [OK] Plot generated
#> [>] [%2d/16] Computing M-estimator influence4
#> [>] [ 4/14] Running sample influence QC analysis (m-estimator)
#>           [OK] M-estimate QC complete
#> [>] [%2d/16] Fitting linear models5
#> [>] [ 5/14] Testing LM interactions with GAM
#> Warning: nlminb problem, convergence error code = 1
#>   message = singular convergence (7)
#>           [OK] LM interaction analysis complete
#> [>] [%2d/16] Plotting LM results6
#> [>] [ 6/14] Plotting LM interaction GAM smoother
#>           [OK] LM interaction plot generated
#> [>] [%2d/16] Computing jackknife isoform switching7
#> [>] [ 7/14] Computing jackknife isoform switching analysis
#>           [OK] Jackknife isoform switching complete
#> [>] [%2d/16] Plotting influence heatmap8
#> [>] [ 9/14] Plotting multi-q influence heatmap
#>           [OK] Influence heatmap generated
#> [>] [%2d/16] Plotting top transcripts9
#> [>] [10/14] Plotting top transcript counts
#>           [OK] Top transcripts plot generated
#> [>] [%2d/16] Computing divergence metrics10
#> [>] [11/14] Computing divergence metrics
#>           [OK] Divergence computed
#> [>] [%2d/16] Computing effect sizes11
#> [>] [12/17] Computing effect sizes for divergence
#>           [OK] Effect sizes computed
#> [>] [%2d/16] Plotting divergence distributions12
#> [>] [13/17] Plotting divergence distribution
#>           [OK] Divergence distribution plot generated
#> [>] [%2d/16] Plotting divergence spectrum13
#> [>] [14/17] Plotting divergence spectrum
#>           [OK] Global divergence spectrum plot generated
#>           [OK] Multi-gene divergence spectrum plot generated
#> [>] [%2d/16] Checking statistical assumptions14
#> [>] [15/17] Validating rank-based test assumptions
#>           [OK] Assumptions validated
#> [>] [%2d/16] Performing Scheirer-Ray-Hare test15
#> [>] [16/17] Running Scheirer-Ray-Hare rank-based test
#> Using condition_col='condition' from @config
#>           [OK] Scheirer-Ray-Hare test completed
#> [>] [%2d/16] Computing LM-rank test concordance16
#> [>] [17/17] Computing concordance between LM and rank test results
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
#>   Total time ........... 25.4s
#>   Slowest steps:
#>     1. lm_interaction       12.6s (49.8%)
#>     2. jackknife            4.5s (17.8%)
#>     3. lm_plot              2.4s (9.5%)
#> 
#> [OUTPUT] Output
#>   Directory ........... tsenat_outputs
#>   Files saved ......... 18
#> 
#> [TIPS] Extract Results - Common Examples:
#>   # View object structure
#>   show(result)
#> 
#>   # View detailed statistics summary
#>   summary(result)
#> 
#>   # Tsallis Entropy Diversity
#>   # Get results for specific sample at q=1.0
#>   div <- results(result, type = 'diversity',
#>                  n_genes = 4, sample = 'SRR14800481')
#> 
#>   # Linear Model Interaction Results
#>   # Top 10 genes by p-value
#>   lm <- results(result, type = 'lm',
#>                 rankBy = 'pvalue', n = 10)
#> 
#>   # Visualizations
#>   plot_diversity <- results(result, type = 'diversity', plot = TRUE)
#>   plot_lm <- results(result, type = 'lm', plot = TRUE)
#> 
# }
```
