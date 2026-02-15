# Extract Transcript-to-Gene Mapping from GFF3 File

Internal function that parses GFF3 or GFF3.gz files to extract
transcript-to-gene mappings. This function is called internally by
[`build_se()`](https://gallardoalba.github.io/TSENAT/reference/build_se.md)
when a GFF3 file path is provided. Users should not call this function
directly.

## Usage

``` r
extract_tx2gene_from_gff3(gff3_file)
```

## Arguments

- gff3_file:

  Path to a GFF3 or GFF3.gz file containing transcript/mRNA features
  with ID and Parent attributes.

## Value

A data.frame with two columns (Transcript, Gene) containing the
transcript-to-gene mapping extracted from the GFF3 file.

## Details

The function expects GFF3 format with the following structure:

- 9 tab-separated columns: seqname, source, feature, start, end, score,
  strand, phase, attributes

- Feature type column (3rd column) should contain 'transcript' or 'mRNA'

- Attributes column (9th column) should contain ID and Parent fields

- ID field: unique identifier for the transcript

- Parent field: references the gene ID that this transcript belongs to

Example GFF3 line:

    chr1\tgencode\ttranscript\t1000\t3000\t.\t+\t.\t
    ID=ENST00000001;Parent=ENSG00000101;Name=BRCA1-001

**Performance:** This function is optimized for large GFF3 files:

- Reads files in 10,000-line chunks (not line-by-line)

- Uses fast pre-filtering (feature type check before regex)

- Employs efficient string operations instead of heavy regex on every
  line

- Handles both compressed (.gz) and uncompressed files seamlessly

## Examples

``` r
if (FALSE) { # \dontrun{
# Create a temporary GFF3 file with transcript features
gff3_lines <- c(
  'chr1\tgencode\ttranscript\t1000\t3000\t.\t+\t.\tID=ENST001;Parent=ENSG001',
  'chr1\tgencode\ttranscript\t1500\t3500\t.\t+\t.\tID=ENST002;Parent=ENSG001',
  'chr2\tgencode\ttranscript\t5000\t8000\t.\t-\t.\tID=ENST003;Parent=ENSG002'
)
tf <- tempfile(fileext = '.gff3')
writeLines(gff3_lines, tf)

# Extract transcript-to-gene mapping
tx2gene <- extract_tx2gene_from_gff3(tf)
head(tx2gene)
} # }
```
