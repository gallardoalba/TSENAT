# Plot top transcripts for a gene

For a given gene, find transcripts using a tx-\>gene mapping, compute
per-transcript statistics between two sample groups, select the top N
transcripts by p-value and plot their expression across groups.

## Arguments

- counts:

  A numeric matrix of transcript-level expression (rows = transcripts,
  columns = samples).

- gene:

  Character; gene symbol to inspect.

- samples:

  Character vector of sample group labels (length = ncol(counts)).

- tx2gene:

  Path to a two-column tab-delimited file with columns \`Transcript\`
  and \`Gen\`, or a data.frame with those columns. Required.

- top_n:

  Integer; number of transcripts to show (default = 3). If NULL, all
  transcripts for the gene are plotted.

- pseudocount:

  Numeric value added when computing log2 fold-change to avoid division
  by zero (default = 1e-6).

- output_file:

  Optional path to save the plot (ggsave will be used). If NULL the
  ggplot object is returned.

## Value

A \`ggplot\` object (or invisibly saved file if \`output_file\`
provided).
