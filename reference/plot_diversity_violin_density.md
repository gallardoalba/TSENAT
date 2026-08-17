# Combined Violin and Density Plot Grid for Single q Value

Creates a side-by-side grid layout with a violin plot on the left and a
density plot on the right, both showing Tsallis entropy distribution for
the q value in the provided SummarizedExperiment (which should contain a
single q value).

## Usage

``` r
plot_diversity_violin_density(
  se,
  assay_name = "diversity",
  title = NULL,
  output_file = NULL,
  q = NULL
)
```

## Arguments

- se:

  A \`SummarizedExperiment\` returned by \`calculate_diversity\`
  containing entropy values at a single q value.

- assay_name:

  Name of the assay to use (default: 'diversity').

- title:

  Optional base title. If NULL, auto-generated based on q value.

- output_file:

  Character or NULL. Optional file path to save the plot as an image. If
  provided, the plot will be saved with appropriate dimensions. Default:
  NULL (no file output, only return object).

## Value

A \`ggplot2\` object showing a 1x2 grid with violin plot on the left and
density plot on the right.

## Examples

``` r
# Plot 8: Violin and density plots of Tsallis entropy distribution
data(readcounts)
metadata_df <- read.table(
  system.file('extdata', 'metadata.tsv', package = 'TSENAT'),
  header = TRUE, sep = '\t'
)
gff3_dataset <- system.file('extdata', 'annotation.gff3.gz', package =
'TSENAT')
readcounts <- as.matrix(readcounts)
mode(readcounts) <- 'numeric'

# Create configuration (required when metadata is provided)
config <- TSENAT_config(sample_col = 'sample', condition_col = 'condition')
analysis <- build_analysis(readcounts = readcounts, tx2gene =
gff3_dataset, metadata = metadata_df, config = config,
  tpm = tpm, effective_length = effective_length)
analysis <- filter_analysis(analysis, min_samples = 1, subset_n_genes
= 200)
analysis <- calculate_diversity(analysis, q = 1.0)
p <- plot_diversity_violin_density(analysis)
# print(p)
```
