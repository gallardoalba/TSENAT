## Helper: Extract tx2gene mapping from GFF3 file
#' @title Extract Transcript-to-Gene Mapping from GFF3 File
#' @description Internal function that parses GFF3 or GFF3.gz files to extract
#' transcript-to-gene mappings. This function is called internally by \code{build_se()}
#' when a GFF3 file path is provided. Users should not call this function directly.
#' @param gff3_file Path to a GFF3 or GFF3.gz file containing transcript/mRNA
#'   features with ID and Parent attributes.
#' @return A data.frame with two columns (Transcript, Gene) containing the
#'   transcript-to-gene mapping extracted from the GFF3 file.
#' @details
#' The function expects GFF3 format with the following structure:
#' \itemize{
#'   \item 9 tab-separated columns: seqname, source, feature, start, end,
#'         score, strand, phase, attributes
#'   \item Feature type column (3rd column) should contain 'transcript' or 'mRNA'
#'   \item Attributes column (9th column) should contain ID and Parent fields
#'   \item ID field: unique identifier for the transcript
#'   \item Parent field: references the gene ID that this transcript belongs to
#' }
#' Example GFF3 line:
#' \preformatted{chr1\tgencode\ttranscript\t1000\t3000\t.\t+\t.\t
#' ID=ENST00000001;Parent=ENSG00000101;Name=BRCA1-001}
#'
#' \strong{Performance:} This function is optimized for large GFF3 files:
#' \itemize{
#'   \item Reads files in 10,000-line chunks (not line-by-line)
#'   \item Uses fast pre-filtering (feature type check before regex)
#'   \item Employs efficient string operations instead of heavy regex on every line
#'   \item Handles both compressed (.gz) and uncompressed files seamlessly
#' }
extract_tx2gene_from_gff3 <- function(gff3_file) {
    # Handle both .gff3 and .gff3.gz files
    if (grepl("\\.gff3\\.gz$", gff3_file)) {
        con <- gzfile(gff3_file, "rt")
    } else {
        con <- file(gff3_file, "r")
    }

    on.exit(close(con))

    # Pre-allocate vectors for efficiency (start with capacity for 10k
    # transcripts) This avoids repeatedly growing lists/data.frames which is
    # slow
    tx2gene_transcripts <- character(10000)
    tx2gene_genes <- character(10000)
    idx <- 0
    chunk_size <- 10000  # Read in chunks for better performance

    while (TRUE) {
        # Read lines in chunks instead of one-by-one for better performance
        lines <- readLines(con, n = chunk_size)
        if (length(lines) == 0)
            break

        # Process each line in the chunk
        for (line in lines) {
            # Skip comments and empty lines (fast pre-filter)
            if (startsWith(line, "#") || line == "")
                next

            # Parse GFF3 line format: seqname source feature start end score
            # strand phase attributes Use strsplit only once per line
            fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
            if (length(fields) < 9)
                next

            # Extract feature type (3rd column) - do this check first to skip
            # early
            feature_type <- fields[3]
            if (!feature_type %in% c("transcript", "mRNA"))
                next

            # Only now extract attributes from the (potentially long) 9th
            # column
            attributes <- fields[9]

            # Extract ID and Parent using efficient substring operations Find
            # positions of ID= and Parent= patterns
            id_start <- regexpr("ID=", attributes, fixed = TRUE) + 3
            if (id_start > 3) {
                # ID found, extract until semicolon or end of string
                id_end <- regexpr(";", substr(attributes, id_start, nchar(attributes)),
                  fixed = TRUE)
                if (id_end > 0) {
                  transcript_id <- substr(attributes, id_start, id_start + id_end -
                    2)
                } else {
                  transcript_id <- substr(attributes, id_start, nchar(attributes))
                }
            } else {
                transcript_id <- NA_character_
            }

            parent_start <- regexpr("Parent=", attributes, fixed = TRUE) + 7
            if (parent_start > 7) {
                # Parent found, extract until semicolon or end of string
                parent_end <- regexpr(";", substr(attributes, parent_start, nchar(attributes)),
                  fixed = TRUE)
                if (parent_end > 0) {
                  gene_id <- substr(attributes, parent_start, parent_start + parent_end -
                    2)
                } else {
                  gene_id <- substr(attributes, parent_start, nchar(attributes))
                }
            } else {
                gene_id <- NA_character_
            }

            # Store mapping if both IDs are present
            if (!is.na(transcript_id) && !is.na(gene_id)) {
                idx <- idx + 1
                # Grow vectors if needed
                if (idx > length(tx2gene_transcripts)) {
                  tx2gene_transcripts <- c(tx2gene_transcripts, rep(NA_character_,
                    10000))
                  tx2gene_genes <- c(tx2gene_genes, rep(NA_character_, 10000))
                }
                tx2gene_transcripts[idx] <- transcript_id
                tx2gene_genes[idx] <- gene_id
            }
        }
    }

    if (idx == 0) {
        stop("No transcript-gene mappings found in GFF3 file. ", "Ensure the file contains 'transcript' or 'mRNA' features with ID and Parent attributes.",
            call. = FALSE)
    }

    # Trim to actual size and create data frame
    tx2gene_df <- data.frame(Transcript = tx2gene_transcripts[1:idx], Gene = tx2gene_genes[1:idx],
        stringsAsFactors = FALSE)
    rownames(tx2gene_df) <- NULL

    return(tx2gene_df)
}

## Helper: build SummarizedExperiment from readcounts + tx2gene
#' Build a SummarizedExperiment from transcript readcounts and tx->gene map
#'
#' This helper creates a `SummarizedExperiment` with an assay named
#' `counts`, stores the `tx2gene` table and the raw `readcounts` in
#' `metadata()`, and extracts gene assignments from the tx2gene mapping
#' to populate `rowData(se)$genes`.
#'
#' The function accepts three different input types for the tx2gene mapping:
#' a file path (TSV or GFF3 format) or an in-memory data.frame.
#'
#' @param readcounts Numeric matrix or data.frame of transcript-level counts
#'   (rows = transcripts, columns = samples). Rownames should contain
#'   transcript IDs that match the first column of tx2gene.
#'
#' @param tx2gene Transcript-to-gene mapping. Can be one of:
#'   \describe{
#'     \item{\strong{TSV File Path}}{A path to a tab-separated file with at least
#'       two columns: Transcript ID (first column) and Gene ID (second column).
#'       The file should have a header row with column names.
#'       File extension should be .tsv, .txt, or similar.
#'       Example:
#'       \preformatted{Transcript\tGene
#'       ENST00000001\tENSG00000101
#'       ENST00000002\tENSG00000102}
#'     }
#'     \item{\strong{GFF3 File Path}}{A path to a GFF3 annotation file (.gff3 or
#'       .gff3.gz format). The file should contain transcript/mRNA features with
#'       ID and Parent attributes. Parent attributes reference gene IDs.
#'       Example GFF3 line:
#'       \preformatted{chr1\tgencode\ttranscript\t100\t2100\t.\t+\t.\t
#'       ID=ENST00000001;Parent=ENSG00000101}
#'       The function automatically detects GFF3 format by file extension
#'       (.gff3 or .gff3.gz) and parses accordingly.
#'     }
#'     \item{\strong{Data Frame}}{An in-memory data.frame with Transcript
#'       and Gene columns. Useful when mapping data is already loaded in R.
#'       Example:
#'       \preformatted{  Transcript        Gene
#'       1 ENST00000001 ENSG00000101
#'       2 ENST00000002 ENSG00000102}
#'     }
#'   }
#'
#' @param assay_name Name for the assay to store readcounts (default: 'counts').
#'
#' @return A `SummarizedExperiment` with assay, `metadata()$tx2gene`,
#'   `metadata()$readcounts` and `rowData(se)$genes` populated.
#'
#' @details
#' \strong{Input Format Detection:}
#' \itemize{
#'   \item If tx2gene is a character string ending in .gff3 or .gff3.gz,
#'         it is parsed as a GFF3 file.
#'   \item If tx2gene is a character string with any other extension or
#'         no extension, it is parsed as a tab-separated file.
#'   \item If tx2gene is a data.frame, it is used directly.
#' }
#'
#' \strong{Transcript Matching:}
#' All transcript IDs in the readcounts object (rownames) must have a
#' corresponding entry in the tx2gene mapping. If any transcript is missing,
#' an error is raised.
#'
#' \strong{Performance:}
#' \itemize{
#'   \item GFF3 files are processed efficiently even for large annotations
#'         (e.g., full GENCODE with 100k+ transcripts)
#'   \item TSV files are standard tab-separated format for fast parsing
#'   \item Data.frame inputs have no I/O overhead
#' }
#'
#' @export
#'
#' @examples
#' # Example 1: Using data.frame (in-memory mapping)
#' tx2gene <- data.frame(
#'     Transcript = c('ENST00000001', 'ENST00000002'),
#'     Gene = c('ENSG00000101', 'ENSG00000102')
#' )
#' readcounts <- matrix(c(10, 5, 2, 3),
#'     nrow = 2,
#'     dimnames = list(c('ENST00000001', 'ENST00000002'), c('s1', 's2'))
#' )
#' se <- build_se(readcounts, tx2gene)
#'
#' # Example 2: Using TSV file path
#' # Assuming you have a file 'tx2gene.tsv' with Transcript and Gene columns
#' # se <- build_se(readcounts, 'path/to/tx2gene.tsv')
#'
#' # Example 3: Using GFF3.gz file path
#' # Assuming you have a file 'annotation.gff3.gz' with transcript features
#' # se <- build_se(readcounts, 'path/to/annotation.gff3.gz')
build_se <- function(readcounts, tx2gene, assay_name = "counts") {
    if (is.character(tx2gene) && length(tx2gene) == 1) {
        if (!file.exists(tx2gene)) {
            stop("tx2gene file not found: ", tx2gene, call. = FALSE)
        }

        # Detect file type and parse accordingly
        if (grepl("\\.gff3(\\.gz)?$", tx2gene, ignore.case = TRUE)) {
            message("Detected GFF3 format. Extracting transcript-to-gene mapping...")
            tx2gene_df <- extract_tx2gene_from_gff3(tx2gene)
        } else {
            # Assume TSV format (backward compatible)
            tx2gene_df <- utils::read.table(tx2gene, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
        }
    } else if (is.data.frame(tx2gene)) {
        tx2gene_df <- tx2gene
    } else {
        stop("'tx2gene' must be a path (TSV or GFF3) or a data.frame.", call. = FALSE)
    }

    if (is.data.frame(readcounts)) {
        readcounts <- as.matrix(readcounts)
    }
    if (!is.matrix(readcounts) || !is.numeric(readcounts)) {
        stop("'readcounts' must be a numeric matrix or numeric data.frame.", call. = FALSE)
    }

    # Extract transcript IDs from rownames or use row indices
    tx_ids <- rownames(readcounts)
    if (is.null(tx_ids)) {
        stop("'readcounts' must have transcript IDs as rownames.", call. = FALSE)
    }

    # Extract gene names from tx2gene based on transcript IDs
    tx_col <- if ("Transcript" %in% colnames(tx2gene_df))
        "Transcript" else colnames(tx2gene_df)[1]
    genes <- tx2gene_df$Gene[match(tx_ids, tx2gene_df[[tx_col]])]

    if (any(is.na(genes))) {
        stop("Some transcript IDs in readcounts were not found in tx2gene mapping.",
            call. = FALSE)
    }

    assays_list <- S4Vectors::SimpleList()
    assays_list[[assay_name]] <- readcounts

    se <- SummarizedExperiment::SummarizedExperiment(assays = assays_list)
    S4Vectors::metadata(se)$tx2gene <- tx2gene_df
    S4Vectors::metadata(se)$readcounts <- readcounts

    # Populate rowData with gene assignments and transcript IDs
    SummarizedExperiment::rowData(se) <- S4Vectors::DataFrame(genes = genes, row.names = tx_ids)

    return(se)
}
