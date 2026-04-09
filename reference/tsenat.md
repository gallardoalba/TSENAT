# Run complete TSENAT analysis pipeline

Coordinates the full TSENAT workflow: diversity -\> jackknife -\> LM
interactions -\> divergence -\> gene interactions -\> visualizations.

## Usage

``` r
tsenat(
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
  [`build_analysis_s4`](https://gallardoalba.github.io/TSENAT/reference/build_analysis_s4.md).

- output_dir:

  `character`. Directory to save results and plots. Default:
  "tsenat_outputs". Set to NULL to disable automatic output saving.

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

1.  [`filter_analysis_s4()`](https://gallardoalba.github.io/TSENAT/reference/filter_analysis_s4.md) -
    Filter low-abundance transcripts

2.  [`calculate_diversity_s4()`](https://gallardoalba.github.io/TSENAT/reference/calculate_diversity_s4.md) -
    Tsallis entropy per q-value

3.  [`plot_tsallis_q_curve_s4()`](https://gallardoalba.github.io/TSENAT/reference/plot_tsallis_q_curve_s4.md) -
    Visualize q-spectrum

4.  [`m_estimate_s4()`](https://gallardoalba.github.io/TSENAT/reference/m_estimate_s4.md) -
    Sample influence QC analysis

5.  [`calculate_lm_interaction_s4()`](https://gallardoalba.github.io/TSENAT/reference/calculate_lm_interaction_s4.md) -
    LM interaction testing

6.  [`plot_lm_interaction_gam_s4()`](https://gallardoalba.github.io/TSENAT/reference/plot_lm_interaction_gam_s4.md) -
    GAM visualization of LM results

7.  [`jackknife_isoform_switching_s4()`](https://gallardoalba.github.io/TSENAT/reference/jackknife_isoform_switching_s4.md) -
    Transcript switching detection

8.  [`prepare_gene_switching_tables_s4()`](https://gallardoalba.github.io/TSENAT/reference/prepare_gene_switching_tables_s4.md) -
    Prepare gene switching summary tables

9.  [`plot_multiq_delta_influence_heatmaps_s4()`](https://gallardoalba.github.io/TSENAT/reference/plot_multiq_delta_influence_heatmaps_s4.md) -
    Multi-q influence heatmap

10. [`plot_top_transcripts_s4()`](https://gallardoalba.github.io/TSENAT/reference/plot_top_transcripts_s4.md) -
    Top transcript visualization

11. [`calculate_divergence_s4()`](https://gallardoalba.github.io/TSENAT/reference/calculate_divergence_s4.md) -
    Pairwise divergence metrics

12. [`effect_sizes_divergence_s4()`](https://gallardoalba.github.io/TSENAT/reference/effect_sizes_divergence_s4.md) -
    Effect size computation

13. [`plot_divergence_distribution_s4()`](https://gallardoalba.github.io/TSENAT/reference/plot_divergence_distribution_s4.md) -
    Divergence distribution plot

14. [`plot_divergence_spectrum_s4()`](https://gallardoalba.github.io/TSENAT/reference/plot_divergence_spectrum_s4.md) -
    Divergence spectrum plot

## Examples

``` r
data(readcounts, package = "TSENAT")
metadata_df <- read.table(
  system.file("extdata", "metadata.tsv", package = "TSENAT"),
  header = TRUE, sep = "\t"
)
gff3_file <- system.file("extdata", "annotation.gff3.gz", package = "TSENAT")

config <- tsenat_config(
  sample_col = "sample",
  condition_col = "condition",
  q_values = c(0.5, 1.0, 1.5, 2.0, 2.5),
  generate_plots = FALSE
)
analysis <- build_analysis_s4(
  readcounts = as.matrix(readcounts),
  tx2gene = gff3_file,
  metadata = metadata_df,
  config = config,
  tpm = tpm,
  effective_length = effective_length
)

result <- tsenat(analysis)
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
#>   Q-spectrum range ..... 0.5 to 2.5 (5 values)
#> 
#> [CONFIG] Analysis Configuration
#>   Design ............... unpaired
#>   Filter stringency .... medium
#>   Normalization ........ enabled [0-1]
#>   Normalization method . NONE
#>   Pseudocount .......... disabled
#>   Shrinkage ............ disabled
#>   Significance ......... p < 0.050 | FDR < 0.050
#>   LM method ............ GAM
#>   LM p-corr method ..... BH
#>   Jackknife use_lm_fdr . TRUE
#> 
#> =============================================================
#> [>] [ 1/14] Filtering low-abundance transcripts
#>           [OK] Complete
#> [>] [ 2/14] Computing Tsallis diversity
#>           [OK] 5 q-values processed
#> [>] [ 3/14] Plotting q-spectrum curve
#>           [OK] Plot generated
#> [>] [ 4/14] Running sample influence QC analysis (m-estimator)
#>           [OK] M-estimate QC complete
#> [>] [ 5/14] Testing LM interactions with GAM smoother
#>           [OK] LM interaction analysis complete
#> [>] [ 6/14] Plotting LM interaction GAM smoother
#> Warning: No valid plots generated
#> [>] [ 7/14] Computing jackknife isoform switching analysis
#>           [OK] Jackknife isoform switching complete
#> [>] [ 8/14] Preparing gene switching tables
#>           [OK] Gene switching tables prepared
#> [>] [ 9/14] Plotting multi-q influence heatmap
#>           [OK] Influence heatmap generated
#> [>] [10/14] Plotting top transcript counts
#>           [OK] Top transcripts plot generated
#> [>] [11/14] Computing divergence metrics
#>           [OK] Divergence computed
#> [>] [12/14] Computing effect sizes for divergence
#>           [OK] Effect sizes computed
#> [>] [13/14] Plotting divergence distribution
#>           [OK] Divergence distribution plot generated
#> [>] [14/14] Plotting divergence spectrum
#>           [OK] Global divergence spectrum plot generated
#>           [OK] Multi-gene divergence spectrum plot generated
#> =============================================================
#> 
#> +============================================================+
#> |               [OK] ANALYSIS COMPLETE                       |
#> +============================================================+
#> 
#> [RESULTS] Results Summary
#> 
#> [PERF] Performance
#>   Total time ........... 17.0s
#>   Slowest steps:
#>     1. lm_interaction       11.7s (69.1%)
#>     2. jackknife            1.3s (7.6%)
#>     3. div_spectrum_plot    1.0s (6.1%)
#> 
#> [OUTPUT] Output
#>   Directory ........... tsenat_outputs
#>   Files saved ......... 15
#> 
#> [TIPS] Next steps:
#>   show(result)             - View object structure and slots
#>   summary(result)          - Print detailed statistics summary
#>   getResults(result)       - Extract numerical results (diversity, divergence, etc.)
#>   getPlot(result, type)    - Retrieve specific visualization (e.g., 'diversity', 'volcano')
#>   getMeta(result)          - Access metadata and workflow parameters
#>   getSE(result)            - Get SummarizedExperiment object for downstream analysis
#> 
```
