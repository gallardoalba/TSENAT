# Build a SummarizedExperiment from transcript readcounts and tx-\>gene map

This helper creates a \`SummarizedExperiment\` with an assay named
\`counts\`, stores the \`tx2gene\` table and the raw \`readcounts\` in
\`metadata()\`, and extracts gene assignments from the tx2gene mapping
to populate \`rowData(se)\$genes\`.

## Usage

``` r
build_se(readcounts, tx2gene, assay_name = "counts")
```

## Arguments

- readcounts:

  Numeric matrix or data.frame of transcript-level counts (rows =
  transcripts, columns = samples). Rownames should contain transcript
  IDs that match the first column of tx2gene.

- tx2gene:

  Transcript-to-gene mapping. Can be one of:

  **TSV File Path**

  :   A path to a tab-separated file with at least two columns:
      Transcript ID (first column) and Gene ID (second column). The file
      should have a header row with column names. File extension should
      be .tsv, .txt, or similar. Example:

          Transcript\tGene
              ENST00000001\tENSG00000101
              ENST00000002\tENSG00000102

  **GFF3 File Path**

  :   A path to a GFF3 annotation file (.gff3 or .gff3.gz format). The
      file should contain transcript/mRNA features with ID and Parent
      attributes. Parent attributes reference gene IDs. Example GFF3
      line:

          chr1\tgencode\ttranscript\t100\t2100\t.\t+\t.\t
              ID=ENST00000001;Parent=ENSG00000101

      The function automatically detects GFF3 format by file extension
      (.gff3 or .gff3.gz) and parses accordingly.

  **Data Frame**

  :   An in-memory data.frame with Transcript and Gene columns. Useful
      when mapping data is already loaded in R. Example:

            Transcript        Gene
              1 ENST00000001 ENSG00000101
              2 ENST00000002 ENSG00000102

- assay_name:

  Name for the assay to store readcounts (default: 'counts').

## Value

A \`SummarizedExperiment\` with assay, \`metadata()\$tx2gene\`,
\`metadata()\$readcounts\` and \`rowData(se)\$genes\` populated.

## Details

The function accepts three different input types for the tx2gene
mapping: a file path (TSV or GFF3 format) or an in-memory data.frame.

**Input Format Detection:**

- If tx2gene is a character string ending in .gff3 or .gff3.gz, it is
  parsed as a GFF3 file.

- If tx2gene is a character string with any other extension or no
  extension, it is parsed as a tab-separated file.

- If tx2gene is a data.frame, it is used directly.

**Transcript Matching:** All transcript IDs in the readcounts object
(rownames) must have a corresponding entry in the tx2gene mapping. If
any transcript is missing, an error is raised.

**Performance:**

- GFF3 files are processed efficiently even for large annotations (e.g.,
  full GENCODE with 100k+ transcripts)

- TSV files are standard tab-separated format for fast parsing

- Data.frame inputs have no I/O overhead

## Examples

``` r
# Example 1: Using data.frame (in-memory mapping)
tx2gene <- data.frame(
    Transcript = c('ENST00000001', 'ENST00000002'),
    Gene = c('ENSG00000101', 'ENSG00000102')
)
readcounts <- matrix(c(10, 5, 2, 3),
    nrow = 2,
    dimnames = list(c('ENST00000001', 'ENST00000002'), c('s1', 's2'))
)
se <- build_se(readcounts, tx2gene)

# Example 2: Using TSV file path
# Assuming you have a file 'tx2gene.tsv' with Transcript and Gene columns
# se <- build_se(readcounts, 'path/to/tx2gene.tsv')

# Example 3: Using GFF3.gz file path
# Assuming you have a file 'annotation.gff3.gz' with transcript features
# se <- build_se(readcounts, 'path/to/annotation.gff3.gz')
```
