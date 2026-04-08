# Add or cache a plot in TSENATAnalysis object

Cache plots for later retrieval, maintaining analysis visualization
history.

## Usage

``` r
addPlot(object, type, plot, replace = FALSE)
```

## Arguments

- object:

  TSENATAnalysis object

- type:

  character. Plot type identifier

- plot:

  ggplot or list. The plot object to cache

- replace:

  logical. If TRUE, replace existing plot of same type. If FALSE
  (default), warn if plot already exists and do not overwrite.

## Value

invisible(object) for method chaining

## Examples

``` r
# Create a simple plot and add to analysis
data(readcounts, package = 'TSENAT')
metadata <- read.table(
  system.file('extdata', 'metadata.tsv', package = 'TSENAT'),
  header = TRUE, sep = '\t'
)
gff3_file <- system.file('extdata', 'annotation.gff3.gz', package = 'TSENAT')

# TPM and effective_length REQUIRED for filter_analysis_s4()
tpm <- matrix(runif(nrow(readcounts) * ncol(readcounts), 0.1, 10),
              nrow = nrow(readcounts), ncol = ncol(readcounts),
              dimnames = dimnames(readcounts))
effective_length <- matrix(100, nrow = nrow(readcounts), ncol = ncol(readcounts))

config <- tsenat_config(q_values = c(0.5, 1.0), generate_plots = FALSE)
analysis <- build_analysis_s4(readcounts, tx2gene = gff3_file,
    metadata = metadata, tpm = tpm, effective_length = effective_length,
    config = config)
analysis <- filter_analysis_s4(analysis, stringency = 'severe')
analysis <- calculate_diversity_s4(analysis, norm = TRUE)
#> Note: 3 genes excluded (< 75% valid values).

# Create and cache a plot
p <- plot_tsallis_q_curve_s4(analysis)
analysis <- addPlot(analysis, type = 'tsallis_q_curve', plot = p)
```
