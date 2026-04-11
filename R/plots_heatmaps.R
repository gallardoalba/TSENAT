#' Plot Divergence Spectrum Heatmaps (Multi-q Transcript Switching)
#'
#' Creates combined heatmap panels showing delta influence (transcript
#' switching magnitude)
#' across multiple q-values (diversity scales) for selected genes. Each
#' heatmap shows
#' how transcript importance differs between conditions (delta_influence)
#' across the
#' q-spectrum from 0.01 to 2.0.
#'
#' @param switching_results A multi-q switching analysis result from
#'   \code{\link{jackknife_isoform_switching}(q = c(...))}. Must include
#'   gene_ids, gene_name_map, and per-gene results for each q-value.
#' @param n_genes Numeric: number of top genes to visualize (default 4).
#'   When \code{lm_results} is provided,  selects the n genes with 
#' lowest p-values.
#'   Otherwise, selects the first n genes from results.
#' @param lm_results Optional data.
#' frame from \code{\link{calculate_lm_interaction}()}
#' containing gene interaction statistics. Should have columns for gene
#' identifiers
#' ('gene_name' or 'gene_id') and p-values ('p_interaction' or
#' 'adj_p_interaction').
#' If provided, genes are ranked by p-value significance for selection of
#' top genes.
#' @param verbose Logical; if TRUE, print detailed validation report of
#' heatmap data
#' including which genes were included and any skipped due to insufficient
#' data.
#'   Default: FALSE (no validation output).
#' @param cellwidth Numeric; width of heatmap cells in pixels (default: 35).
#'   Following pheatmap best practices for publication-quality heatmaps.
#' Larger values (50+) make cells more visible but reduce number of visible
#' transcripts.
#' @param cellheight Numeric; height of heatmap cells in pixels (default:
#' 10.25).
#' Following pheatmap best practices. Smaller values allow more q-values to
#' be visible.
#' @param fontsize Numeric; font size in points for heatmap labels (default:
#' 11).
#' Following pheatmap best practices for publication-quality figures.
#' Applies to
#'   row labels (q-values) and column labels (transcript IDs).
#' @param width Output image width in inches. If NULL, automatically
#' calculated (12 inches).
#' @param height Output image height in inches. If NULL, automatically
#' calculated based on number of layout rows.
#'
#' @return Character path to saved PNG file containing the combined heatmaps.
#'   The plot is automatically saved to a temporary file and can be displayed
#'   in R Markdown with \code{knitr::include_graphics()}.
#'
#' @details
#' **Heatmap interpretation:**
#'
#' - \bold{Rows}: Different q-values from 0.01 (rare isoforms) to 2.0
#' (dominant isoforms)
#' - \bold{Columns}: Individual transcripts of each gene
#' - \bold{Color scale}: Blue (negative delta_influence) = transcript more
#' important in second condition;
#' Red (positive delta_influence) = transcript more important in first
#' condition;
#'   White = no switching effect
#' - \bold{Intensity}: Darker colors indicate stronger switching magnitude
#'
#' **Layout:**
#' - Multiple genes displayed in separate panels (up to 2 per row)
#' - Panels combined into single PNG for reproducible visualization
#' - Outliers (>95th percentile) are capped for better color contrast
#' - NaN and Inf values treated as missing (light gray)
#'
#' **Use case:**
#' Identify whether transcript switching is consistent across diversity scales
#' or scale-dependent. Genes with similar patterns across q-values show robust
#' isoform shifts; genes with varying patterns across q indicate q-dependent
#' switching driven by rare vs. abundant isoforms.
#'
#' @examples
#' # Example: Create synthetic multi-q switching results
#' # For real analysis, use .calculate_jis() output
#' set.seed(123)
#' gene_names <- paste0('gene_', 1:4)
#' names(gene_names) <- 1:4
#' 
#' # Create multi-q results structure
#' q_values <- c(0.5, 1.0, 1.5)
#' switching_results <- structure(
#'   list(
#'     q_0_50 = list(
#'       gene_ids = 1:4,
#'       gene_name_map = gene_names,
#'       switching_results = list(
#'         list(gene_id = 1, delta_influence = matrix(rnorm(30, 0, 0.2), 5, 6)),
#'         list(gene_id = 2, delta_influence = matrix(rnorm(30, 0, 0.2), 5, 6)),
#'         list(gene_id = 3, delta_influence = matrix(rnorm(30, 0, 0.2), 5, 6)),
#'         list(gene_id = 4, delta_influence = matrix(rnorm(30, 0, 0.2), 5, 6))
#'       )
#'     )
#'   ),
#'   class = 'tsenat_isoform_switching_multiq'
#' )
#' 
#' # Create heatmap visualization
#' .plot_jis_delta(switching_results, n_genes = 2)
#'
#' @import grid
#' @import pheatmap
#' @importFrom grDevices png dev.off colorRampPalette

#' @noRd

.plot_jis_delta <- function(switching_results, n_genes = 4,
    lm_results = NULL, verbose = FALSE, cellwidth = 0, cellheight = 0, fontsize = 18,
    layout_ncol = 2, output_file = NULL, width = NULL, height = NULL) {
    # Phase 1: Validate input
    result_data <- .validate_multiq_input(switching_results)
    q_result_keys <- result_data$q_result_keys
    gene_ids <- result_data$gene_ids
    gene_name_map <- result_data$gene_name_map

    # Phase 2: Select genes
    top_genes <- .heatmap_select_genes_multiq(switching_results, n_genes, lm_results)

    # Phase 3: Collect gene info for layout planning
    gene_info_list <- lapply(seq_along(top_genes), function(i) {
        mat <- .heatmap_prepare_multiq_data(switching_results, top_genes[i], q_result_keys)
        if (is.null(mat)) {
            list(n_transcripts = 0)
        } else {
            list(n_transcripts = ncol(mat))
        }
    })

    # Phase 4: Plan layout and calculate dimensions
    layout_result <- .plot_adaptive_layout(gene_info_list, use_fixed_layout = !is.null(layout_ncol) &&
        layout_ncol > 0, layout_ncol = layout_ncol)
    gene_layout <- layout_result$layout
    n_layout_rows <- layout_result$n_layout_rows
    dims <- .calculate_heatmap_dimensions(n_layout_rows, length(q_result_keys), width_in = if (is.null(width)) 12 else width, height_in = height)

    # Phase 5: Create heatmaps (using refactored loop)
    all_gene_matrices <- list()
    all_gene_info <- list()

    for (gene_idx in seq_along(top_genes)) {
        gene_id <- top_genes[gene_idx]

        # Try to find the gene in any of the q-value results Build list of all
        # available gene names from all q-values
        available_genes <- unique(unlist(lapply(switching_results[q_result_keys],
            function(qres) {
                if (!is.null(qres$results_per_gene))
                  names(qres$results_per_gene) else NULL
            })))

        # Check if gene_id matches - if not, try to find a match by partial
        # match or case-insensitive match
        actual_gene_id <- gene_id
        if (!(gene_id %in% available_genes) && length(available_genes) > 0) {
            # Try case-insensitive match
            matches <- grep(paste0("^", tolower(gene_id), "$"), tolower(available_genes))
            if (length(matches) > 0) {
                actual_gene_id <- available_genes[matches[1]]
            }
        }

        mat <- .heatmap_prepare_multiq_data(switching_results, actual_gene_id, q_result_keys)

        if (is.null(mat) || nrow(mat) == 0 || ncol(mat) == 0) {
            all_gene_matrices[[gene_idx]] <- NULL
            all_gene_info[[gene_idx]] <- NULL
            next
        }

        # Get gene name
        gene_name_idx <- which(gene_ids == actual_gene_id)[1]
        gene_name <- if (!is.na(gene_name_idx) && !is.na(gene_name_map[gene_name_idx])) {
            gene_name_map[gene_name_idx]
        } else {
            actual_gene_id
        }

        all_gene_matrices[[gene_idx]] <- mat
        all_gene_info[[gene_idx]] <- list(gene_id = actual_gene_id, gene_name = gene_name,
            n_transcripts = ncol(mat))
    }

    # Security check: if all matrices are NULL, we can't proceed BUT: Allow
    # partial data - if at least SOME genes have data, continue
    non_null_count <- sum(!vapply(all_gene_matrices, is.null, logical(1)))
    if (non_null_count == 0) {
        available_genes <- unique(unlist(lapply(switching_results[q_result_keys],
            function(qres) {
                if (!is.null(qres$results_per_gene))
                  names(qres$results_per_gene) else NULL
            })))

        error_msg <- paste0("No valid heatmap data generated for any genes. ", "Requested genes: ",
            paste(top_genes, collapse = ", "), ". ", "Available genes: ", paste(head(available_genes,
                5), collapse = ", "), if (length(available_genes) > 5)
                "..." else "")

        warning(error_msg, call. = FALSE)
        return(invisible(NULL))
    }

    # Phase 6: Create pheatmap objects
    if (!requireNamespace("pheatmap", quietly = TRUE)) {
        stop("pheatmap package required for this function. Install with: install.packages('pheatmap')",
            call. = FALSE)
    }

    heatmap_plots <- list()
    for (gene_idx in seq_along(top_genes)) {
        if (is.null(all_gene_matrices[[gene_idx]]) || is.null(all_gene_info[[gene_idx]])) {
            heatmap_plots[[gene_idx]] <- NULL
            next
        }

        mat <- all_gene_matrices[[gene_idx]]
        gene_info <- all_gene_info[[gene_idx]]
        layout_info <- gene_layout[[gene_idx]]
        width_frac <- if (!is.null(layout_info))
            layout_info$width else 1

        # Calculate cell sizes
        cells <- .calculate_adaptive_cellsizes(ncol(mat), nrow(mat), width_frac,
            cellwidth, cellheight, fontsize)

        # Create pheatmap grob
        heatmap_plots[[gene_idx]] <- .create_pheatmap_grob(mat, title = gene_info$gene_name,
            cellw = cells$cellwidth, cellh = cells$cellheight, fontsize = cells$fontsize_adj,
            cluster_rows = FALSE)
    }

    # Phase 7: Render grid Filter to only non-NULL heatmaps for rendering
    non_null_idx <- !vapply(heatmap_plots, is.null, logical(1))
    heatmap_plots_filtered <- heatmap_plots[non_null_idx]
    gene_layout_filtered <- if (!is.null(gene_layout) && length(gene_layout) > 0) {
        gene_layout[non_null_idx]
    } else {
        NULL
    }

    tryCatch({
        .plot_grid_setup(n_layout_rows, output_file, dims$png_width, dims$png_height,
            title = "Delta Influence Across Diversity Scales", subtitle = "Jackknife weights across q-spectrum for selected genes")

        .render_heatmaps_to_grid(heatmap_plots_filtered, gene_layout_filtered, layout_ncol)

        .plot_grid_finalize(output_file, verbose = verbose)
    }, error = function(e) {
        if (!is.null(output_file)) {
            tryCatch(grDevices::dev.off(), silent = TRUE)
        }
        stop("Heatmap creation failed: ", e$message, call. = FALSE)
    })
}

#' Plot top transcripts for a gene using pheatmap
#' @param se A `SummarizedExperiment` with transcript counts as assay and
#' gene information in rowData.
#' Must have a 'genes' column in rowData specifying which gene each
#' transcript belongs to.
#' If `use_tpm = TRUE`, requires TPM data in metadata (provided to
#' `build_analysis()` or `.build_se()`).
#' @param gene Character vector; gene symbol(s) to inspect. If NULL and
#' `res` is provided,
#'   top genes are selected by p-value.
#' @param condition_col Character; column name in colData(se) to use for
#' sample grouping
#'   (default: 'sample_type').
#' @param res Optional result data.frame from differential/interaction
#' analysis with gene identifiers and p-values.
#'   Supported sources:
#' - `.calculate_lm(..., return_model_data = TRUE)` returns a
#' list with $results and $model_data
#' - `.calculate_lm(..., return_model_data = FALSE)` returns a
#' data.frame with adj_p_interaction column
#' - `.calculate_rank_test()` returns a data.frame with adj_p_value column
#' (for Scheirer-Ray-Hare rank tests)
#'   If provided and `gene` is NULL, top genes are selected by adjusted p-value.
#' @param top_n Integer number of transcripts to show (default = 3). Use
#' NULL to plot all transcripts for the gene.
#' @param output_file Optional file path to save the plot. If `NULL`,
#' renders to active graphics device.
#' @param metric Aggregation metric: 'median', 'mean', 'variance', or 'iqr'
#' (default: 'median').
#' @param use_tpm Logical; if TRUE, uses TPM (Transcripts Per Million) from
#' metadata instead of raw counts
#' (default: FALSE). TPM is normalized for sequencing depth and is
#' recommended for comparing
#' expression across samples. Requires TPM data in `metadata(se)$tpm`
#' from `build_analysis()` or `.build_se()`
#' with `tpm` parameter. Raises error if TPM not available and `use_tpm =
#' TRUE`.
#' @param width Output image width in inches. If NULL, automatically
#' calculated (12 inches).
#' @param height Output image height in inches. If NULL, automatically
#' calculated based on number of genes.
#' @param fontsize Base font size for heatmap titles and labels (default:
#' 16pt). Automatically scaled for readability.
#' @param cellwidth Width of individual heatmap cells in pixels (default: 0
#' = adaptive). Set > 0 to use fixed sizing.
#' @param cellheight Height of individual heatmap cells in pixels (default:
#' 0 = adaptive). Set > 0 to use fixed sizing.
#' @param layout_ncol Number of heatmaps per row in fixed layout (default:
#' 2). If NULL, uses adaptive layout based on transcript counts.
#' @return Invisibly returns the output file path (if `output_file`
#' provided), or invisible(NULL) if rendering to active device.
#' Graphics are rendered to the active grid device for capture during
#' vignette compilation.
#' @details
#' Visualizes transcript abundances across conditions as pheatmap heatmaps
#' (one per gene with conditions as columns).
#' Aggregates expression by condition using the specified metric before
#' log-normalization with pseudocount.
#' Uses hierarchical clustering of transcripts and condition-based samples.
#' Following pheatmap best practices:
#' publication-quality colors, dynamic cell sizing, and no artificial gaps
#' between cells.
#'
#' Architecture follows the pattern established by
#' `.plot_jis_delta()`:
#' - Phase 1: Input validation and extraction
#' - Phase 2: Gene/condition selection
#' - Phase 3: Layout planning (before creating heatmaps)
#' - Phase 4: Data preparation and heatmap creation
#' - Phase 5: Grid layout and rendering
#' @examples
#' library(SummarizedExperiment)
#' library(S4Vectors)
#' # Create example SummarizedExperiment
#' counts <- matrix(sample(1:100, 36, replace = TRUE), nrow = 6, ncol = 6)
#' rownames(counts) <- paste0('tx', 1:6)
#' rowData_df <- DataFrame(genes = rep(paste0('G', 1:3), each = 2))
#' colData_df <- DataFrame(sample_type = rep(c('Normal', 'Tumor'), 3))
#' se <- SummarizedExperiment(assays = list(counts = counts), 
#'                           rowData = rowData_df, colData = colData_df)
#' # Plot top transcripts
#' .plot_expression(se, gene = 'G1', top_n = 2, output_file =
#' '/tmp/heatmap.png')

#' @noRd

.plot_expression <- function(se, gene = NULL, condition_col = "condition", res = NULL,
    top_n = 3, output_file = NULL, metric = c("median", "mean", "variance", "iqr"),
    use_tpm = TRUE, width = NULL, height = NULL, fontsize = 16, cellwidth = 0, cellheight = 0,
    layout_ncol = 2) {
    if (!requireNamespace("pheatmap", quietly = TRUE)) {
        stop("pheatmap package required", call. = FALSE)
    }

    # Phase 1: Validate input and extract components
    se_data <- .validate_se_for_heatmaps(se, condition_col = condition_col)
    counts <- se_data$counts
    rd <- se_data$rowdata
    cd <- se_data$coldata
    gene_col <- se_data$gene_col

    # Handle TPM override if available
    if (use_tpm) {
        md <- S4Vectors::metadata(se)
        if (!is.null(md$tpm)) {
            tpm_data <- as.matrix(md$tpm)
            if (nrow(tpm_data) == nrow(counts) && ncol(tpm_data) == ncol(counts)) {
                counts <- tpm_data
            }
        }
    }

    # Build tx2gene mapping
    tx2gene <- data.frame(Transcript = rownames(counts), Gen = as.character(rd[[gene_col]]),
        stringsAsFactors = FALSE)

    conditions <- as.character(cd[[condition_col]])
    unique_conditions <- unique(conditions)

    # Phase 2: Select genes
    metric_choice <- match.arg(metric)
    if (is.null(gene) && !is.null(res)) {
        gene <- .heatmap_select_genes_results(se, res, gene_col, top_n, tx2gene)
    }
    if (is.null(gene)) {
        stop("gene must be provided or derivable from res", call. = FALSE)
    }

    # Resolve gene identifiers: convert gene names/transcript IDs to gene IDs
    gene <- .resolve_gene_identifiers(gene, tx2gene, rd, gene_col)
    
    if (FALSE) {  # Debug mode - set to TRUE if needed
        message("[DEBUG] After resolution, genes: ", paste(gene, collapse=", "))
        message("[DEBUG] tx2gene$Gen unique values (first 10): ", paste(head(unique(tx2gene$Gen), 10), collapse=", "))
    }

    # Phase 3: Plan layout
    gene_info_list <- lapply(seq_along(gene), function(i) {
        tx_idx <- which(tx2gene$Gen == gene[i])
        list(n_transcripts = length(tx_idx))
    })

    layout_result <- .plot_adaptive_layout(gene_info_list, use_fixed_layout = !is.null(layout_ncol) &&
        layout_ncol > 0, layout_ncol = layout_ncol)
    n_layout_rows <- layout_result$n_layout_rows
    gene_layout <- layout_result$layout

    dims <- list(png_width = if (is.null(width)) 12 else width, png_height = if (is.null(height)) 3.5 *
        n_layout_rows + 2.5 else height)

    # Phase 4: Create heatmaps
    heatmap_plots <- list()
    for (gene_idx in seq_along(gene)) {
        gene_name <- gene[gene_idx]
        tx_indices <- which(tx2gene$Gen == gene_name)

        if (length(tx_indices) == 0) {
            heatmap_plots[[gene_idx]] <- NULL
            next
        }

        gene_counts <- counts[tx_indices, , drop = FALSE]
        mat <- .heatmap_prepare_condition_data(gene_counts, seq_along(tx_indices),
            conditions, metric_choice)

        if (is.null(mat) || nrow(mat) == 0) {
            heatmap_plots[[gene_idx]] <- NULL
            next
        }

        # Calculate cell sizes
        layout_info <- gene_layout[[gene_idx]]
        width_frac <- if (!is.null(layout_info))
            layout_info$width else 1
        cells <- .calculate_adaptive_cellsizes(ncol(mat), nrow(mat), width_frac,
            cellwidth, cellheight, fontsize)

        # Create pheatmap
        heatmap_plots[[gene_idx]] <- tryCatch({
            .create_pheatmap_grob(mat, title = gene_name,
                cellw = cells$cellwidth, cellh = cells$cellheight, fontsize = cells$fontsize_adj,
                cluster_rows = FALSE)
        }, error = function(e) {
            NULL
        })
    }

    if (all(vapply(heatmap_plots, is.null, logical(1)))) {
        stop("No valid heatmaps created", call. = FALSE)
    }

    # Phase 5: Render grid
    tryCatch({
        metric_label <- if (metric_choice == "iqr")
            "IQR" else metric_choice

        .plot_grid_setup(n_layout_rows, output_file, dims$png_width, dims$png_height,
            title = "Isoform Expression Profiles", subtitle = paste("Log2-normalized",
                metric_label, "by condition"))

        .render_heatmaps_to_grid(heatmap_plots, gene_layout, layout_ncol)
        .plot_grid_finalize(output_file, verbose = FALSE)
    }, error = function(e) {
        if (!is.null(output_file)) {
            tryCatch(grDevices::dev.off(), silent = TRUE)
        }
        stop("Heatmap rendering failed: ", e$message, call. = FALSE)
    })
}

# ============================================================================
# Internal Helper Functions for Heatmap Refactoring
# ============================================================================
# This file contains shared helper functions extracted to support
# .plot_jis_delta() and .plot_expression()
# refactoring to meet Bioconductor's 50-line function guideline.  All functions
# marked @keywords internal @noRd are NOT exported.
# ============================================================================

# ============================================================================
# SECTION 1: INPUT VALIDATION HELPERS
# ============================================================================

#' Validate Multi-Q Isoform Switching Results
#'
#' Checks that switching_results is a valid multi-q result object and
#' extracts necessary components (q-values, gene information).
#'
#' @param switching_results Object to validate
#'
#' @return List with:
#' - q_result_keys: character vector of q-value result keys (q_0_01, q_0_50,
#' etc.)
#'   - first_result: the results for the first q-value
#'   - gene_ids: unique gene IDs from first result
#'   - gene_name_map: gene ID to name mapping
#'

#' @noRd
.validate_multiq_input <- function(switching_results) {
    if (!inherits(switching_results, "tsenat_isoform_switching_multiq")) {
        stop("switching_results must be a multi-q result from .calculate_jis()",
            call. = FALSE)
    }

    q_result_keys <- names(switching_results)[grepl("^q_", names(switching_results))]
    if (length(q_result_keys) == 0) {
        stop("No multi-q results found in switching_results", call. = FALSE)
    }

    first_result <- switching_results[[q_result_keys[1]]]
    gene_ids <- first_result$gene_ids
    gene_name_map <- first_result$gene_name_map

    if (length(gene_ids) == 0) {
        stop("No genes found in switching_results", call. = FALSE)
    }

    list(q_result_keys = q_result_keys, first_result = first_result, gene_ids = gene_ids,
        gene_name_map = gene_name_map)
}

#' Validate SummarizedExperiment Input
#'
#' Checks that se is valid and has required columns for condition and gene info.
#'
#' @param se SummarizedExperiment object to validate
#' @param gene_col Character: name of gene column in rowData
#' @param condition_col Character: name of condition column in colData
#'
#' @return List with:
#'   - counts: assay matrix
#'   - rowdata: rowData as data.frame
#'   - coldata: colData as data.frame
#'   - gene_col: validated gene column name
#'   - condition_col: validated condition column name
#'

#' @keywords internal
#' @noRd
.resolve_gene_identifiers <- function(genes, tx2gene, rd, gene_col) {
    # Flexible gene identifier resolution
    # Accepts: gene IDs, gene names, or transcript IDs
    # Returns: vector of gene IDs for lookup in tx2gene
    
    if (is.null(genes) || length(genes) == 0) {
        return(genes)
    }
    
    genes <- as.character(genes)
    
    # Get available identifiers from rowData
    gene_ids <- if ("gene_id" %in% colnames(rd)) {
        as.character(rd$gene_id)
    } else {
        NULL
    }
    gene_names <- if ("gene_name" %in% colnames(rd)) {
        as.character(rd$gene_name)
    } else {
        NULL
    }
    transcript_ids <- tx2gene$Transcript
    
    # Key: tx2gene$Gen is built from rd[[gene_col]], so we need to resolve TO that column
    # If gene_col="gene_name", we need to convert gene_ids to gene_names
    target_col <- as.character(rd[[gene_col]])
    
    # Try to resolve each gene
    resolved_genes <- character(length(genes))
    
    for (i in seq_along(genes)) {
        gene_input <- genes[i]
        
        # Direct match: already in target column
        if (gene_input %in% target_col) {
            resolved_genes[i] <- gene_input
            next
        }
        
        # Is it a gene_id that needs mapping to target_col?
        if (!is.null(gene_ids) && gene_input %in% gene_ids) {
            idx <- which(gene_ids == gene_input)[1]
            resolved_genes[i] <- target_col[idx]
            next
        }
        
        # Is it a transcript ID?
        if (gene_input %in% transcript_ids) {
            idx <- which(transcript_ids == gene_input)[1]
            resolved_genes[i] <- tx2gene$Gen[idx]
            next
        }
        
        # If not found, keep original (will fail downstream with informative error)
        resolved_genes[i] <- gene_input
    }
    
    resolved_genes
}

#' @noRd
.validate_se_for_heatmaps <- function(se, gene_col = NULL, condition_col = NULL) {
    if (!inherits(se, "SummarizedExperiment")) {
        stop("se must be a SummarizedExperiment object", call. = FALSE)
    }

    # Validate/auto-detect gene column (prefer gene_name for human readability)
    if (is.null(gene_col)) {
        gene_col <- if ("gene_name" %in% colnames(rowData(se))) {
            "gene_name"
        } else if ("genes" %in% colnames(rowData(se))) {
            "genes"
        } else if ("gene_id" %in% colnames(rowData(se))) {
            "gene_id"
        } else {
            stop("rowData(se) must contain 'gene_name', 'genes', or 'gene_id' column",
                call. = FALSE)
        }
    } else {
        if (!gene_col %in% colnames(rowData(se))) {
            stop("gene_col '", gene_col, "' not found in rowData(se)", call. = FALSE)
        }
    }

    # Validate condition column (optional)
    if (!is.null(condition_col)) {
        if (!condition_col %in% colnames(colData(se))) {
            stop("condition_col '", condition_col, "' not found in colData(se)",
                call. = FALSE)
        }
    }

    list(counts = as.matrix(assay(se)), rowdata = as.data.frame(rowData(se)), coldata = as.data.frame(colData(se)),
        gene_col = gene_col, condition_col = condition_col)
}

# ============================================================================
# SECTION 2: GENE SELECTION HELPERS
# ============================================================================

#' Select Top Genes from Multi-Q Results
#'
#' Ranks genes using LM results (if provided) and selects top N genes
#' for visualization from multi-q switching analysis results.
#'
#' @param switching_results Multi-q result list (from
#' jackknife_isoform_switching)
#' @param n_genes Integer: number of top genes to select
#' @param lm_results Optional data.frame with LM interaction results
#'
#' @return Character vector of selected gene IDs
#'

#' @noRd
.heatmap_select_genes_multiq <- function(switching_results, n_genes = 4, lm_results = NULL) {
    q_key <- names(switching_results)[grepl("^q_", names(switching_results))][1]
    first_result <- switching_results[[q_key]]
    gene_ids <- first_result$gene_ids

    if (is.null(lm_results)) {
        # No ranking: use first N genes
        return(gene_ids[seq_len(min(n_genes, length(gene_ids)))])
    }

    # Find p-value column
    p_col <- if ("adj_p_interaction" %in% colnames(lm_results)) {
        "adj_p_interaction"
    } else if ("p_interaction" %in% colnames(lm_results)) {
        "p_interaction"
    } else {
        NULL
    }

    if (is.null(p_col)) {
        return(gene_ids[seq_len(min(n_genes, length(gene_ids)))])
    }

    # Find gene identifier column in lm_results
    lm_gene_col <- if ("gene_id" %in% colnames(lm_results)) {
        "gene_id"
    } else if ("gene" %in% colnames(lm_results)) {
        "gene"
    } else if ("gene_name" %in% colnames(lm_results)) {
        "gene_name"
    } else {
        NULL
    }

    if (is.null(lm_gene_col)) {
        return(gene_ids[seq_len(min(n_genes, length(gene_ids)))])
    }

    # Rank genes by p-value (lowest = most significant)
    matches <- match(gene_ids, lm_results[[lm_gene_col]])
    p_values <- rep(Inf, length(gene_ids))
    matched_idx <- !is.na(matches)
    p_values[matched_idx] <- lm_results[[p_col]][matches[matched_idx]]

    gene_order <- order(p_values)
    genes_sorted <- gene_ids[gene_order]
    genes_sorted[seq_len(min(n_genes, length(genes_sorted)))]
}

#' Select Top Genes from Results DataFrame
#'
#' Extracts and ranks genes from a results data.frame (typically from
#' calculate_lm_interaction or similar statistical results), and returns
#' top N genes found in the SummarizedExperiment.
#'
#' @param se SummarizedExperiment object
#' @param res Results data.frame with gene and p-value columns
#' @param gene_col Character: gene column name in rowData(se)
#' @param top_n Integer: number of top genes to select
#' @param tx2gene Data.frame mapping transcript IDs to gene IDs
#'
#' @return Character vector of selected gene IDs
#'

#' @noRd
.heatmap_select_genes_results <- function(se, res, gene_col = "genes", top_n = 3,
    tx2gene = NULL) {
    if (!is.data.frame(res)) {
        stop("res must be a data.frame", call. = FALSE)
    }

    # Build tx2gene if not provided
    if (is.null(tx2gene)) {
        rd <- as.data.frame(rowData(se))
        tx2gene <- data.frame(Transcript = rownames(se), Gen = as.character(rd[[gene_col]]),
            stringsAsFactors = FALSE)
    }

    # Find gene column in results
    res_gene_col <- if ("gene" %in% colnames(res)) {
        "gene"
    } else if ("genes" %in% colnames(res)) {
        "genes"
    } else if ("gene_id" %in% colnames(res)) {
        "gene_id"
    } else if ("gene_name" %in% colnames(res)) {
        "gene_name"
    } else {
        NULL
    }

    if (is.null(res_gene_col)) {
        stop("res must contain 'gene', 'genes', 'gene_id', or 'gene_name' column",
            call. = FALSE)
    }

    # Find p-value column
    p_col <- if ("adj_p_value" %in% colnames(res)) {
        "adj_p_value"
    } else if ("adj_p_interaction" %in% colnames(res)) {
        "adj_p_interaction"
    } else if ("padj" %in% colnames(res)) {
        "padj"
    } else if ("adjusted_p_values" %in% colnames(res)) {
        "adjusted_p_values"
    } else if ("p_value" %in% colnames(res)) {
        "p_value"
    } else if ("p_interaction" %in% colnames(res)) {
        "p_interaction"
    } else {
        NULL
    }

    if (!is.null(p_col)) {
        res_sorted <- res[order(res[[p_col]], na.last = NA), ]
    } else {
        # No p-value column: use order as-is
        res_sorted <- res
    }

    # Extract genes and filter to those in SE
    se_genes <- unique(tx2gene$Gen)
    selected_genes <- as.character(res_sorted[[res_gene_col]])
    selected_genes <- selected_genes[selected_genes %in% se_genes]
    selected_genes[seq_len(min(top_n, length(selected_genes)))]
}

# ============================================================================
# SECTION 3: LAYOUT PLANNING HELPERS
# ============================================================================

#' Plan Adaptive or Fixed Grid Layout for Heatmaps
#'
#' Determines the position of each heatmap in a grid layout, supporting
#' both fixed layouts (N columns per row) and adaptive layouts (based on
#' data dimensions). Returns layout information for each gene/heatmap.
#'
#' @param gene_info_list List where each element is a gene with info including
#'   n_transcripts (or similar data row count). Each element should have
#'   a $n_transcripts or $n_cols field.
#' @param use_fixed_layout Logical: if TRUE, use layout_ncol columns per row
#' @param layout_ncol Integer: number of columns per row (for fixed layout)
#'
#' @return List with:
#'   - $layout: named list where each element has $row, $col, $width
#'   - $n_layout_rows: total number of rows in layout
#'   - $row_heights: relative heights for each row (if needed)
#'

#' @noRd
.plot_adaptive_layout <- function(gene_info_list, use_fixed_layout = TRUE, layout_ncol = 2) {
    n_genes <- length(gene_info_list)

    if (use_fixed_layout && layout_ncol > 0) {
        # FIXED LAYOUT: layout_ncol genes per row
        n_cols <- as.integer(layout_ncol)
        gene_layout <- list()
        grid_row <- 0

        for (i in seq_len(n_genes)) {
            col_pos <- ((i - 1)%%n_cols) + 1
            if (col_pos == 1)
                grid_row <- grid_row + 1
            gene_layout[[i]] <- list(row = grid_row, col = col_pos, width = 1/n_cols)
        }
        n_layout_rows <- grid_row
    } else {
        # ADAPTIVE LAYOUT: genes with >5 transcripts get full width
        gene_layout <- list()
        n_layout_rows <- 0
        i <- 1

        while (i <= n_genes) {
            n_tx_i <- gene_info_list[[i]]$n_transcripts %||% 1
            has_next <- i < n_genes
            n_tx_next <- if (has_next) {
                gene_info_list[[i + 1]]$n_transcripts %||% 1
            } else {
                0
            }

            if (n_tx_i > 5) {
                # Full-width row
                gene_layout[[i]] <- list(row = n_layout_rows + 1, col = 1, width = 1)
                n_layout_rows <- n_layout_rows + 1
                i <- i + 1
            } else if (n_tx_i <= 5 && has_next && n_tx_next <= 5) {
                # Pair two small genes (half-width each)
                gene_layout[[i]] <- list(row = n_layout_rows + 1, col = 1, width = 0.5)
                gene_layout[[i + 1]] <- list(row = n_layout_rows + 1, col = 2, width = 0.5)
                n_layout_rows <- n_layout_rows + 1
                i <- i + 2
            } else {
                # Single full-width row
                gene_layout[[i]] <- list(row = n_layout_rows + 1, col = 1, width = 1)
                n_layout_rows <- n_layout_rows + 1
                i <- i + 1
            }
        }
    }

    list(layout = gene_layout, n_layout_rows = n_layout_rows, row_heights = rep(c(1,
        0.15), n_layout_rows)[seq_len(n_layout_rows * 2 - 1)])
}

#' Calculate Heatmap Output Dimensions
#'
#' Calculates height for PNG output based on number of layout rows
#' and data rows (q-values, conditions, etc.), following pheatmap
#' and publication best practices.
#'
#' @param n_layout_rows Integer: number of rows in grid layout
#' @param n_data_rows Integer: number of data rows per heatmap (q-values,
#' conditions)
#' @param width_in Numeric: desired width in inches (default 12)
#' @param height_in Numeric: optional fixed height in inches
#'
#' @return List with:
#'   - $png_width: width in inches (default 12, ~1200px @ 100 DPI)
#'   - $png_height: height in inches
#'   - $heatmap_height: height allocated for heatmaps only
#'

#' @noRd
.calculate_heatmap_dimensions <- function(n_layout_rows, n_data_rows, width_in = 12,
    height_in = NULL) {
    png_width <- width_in

    # Scale height: 3 inches per layout row + gaps
    height_per_layout_row <- 3 * (n_data_rows/5)
    gap_between_rows <- 1.8
    heatmap_height <- height_per_layout_row * n_layout_rows + gap_between_rows *
        (n_layout_rows - 1)

    # Total PNG height includes title/subtitle space
    total_height <- if (is.null(height_in)) {
        heatmap_height + 2.5
    } else {
        height_in
    }

    list(png_width = png_width, png_height = total_height, heatmap_height = heatmap_height)
}

# ============================================================================
# SECTION 4: CELL SIZING HELPER
# ============================================================================

#' Calculate Adaptive Cell Dimensions for Pheatmap
#'
#' Calculates pheatmap cellwidth and cellheight based on matrix dimensions,
#' available grid space, and layout context. Supports both full-width and
#' half-width (2-per-row) layouts.
#'
#' @param n_cols_mat Integer: number of columns in data matrix (transcripts)
#' @param n_rows_mat Integer: number of rows in data matrix (q-values,
#' conditions)
#' @param width_frac Numeric: fraction of plot width allocated (1 or 0.5)
#' @param cellwidth Numeric: override cell width (0 = adaptive)
#' @param cellheight Numeric: override cell height (0 = adaptive)
#' @param fontsize Numeric: base font size in points
#'
#' @return List with:
#'   - $cellwidth: width for pheatmap cellwidth parameter
#'   - $cellheight: height for pheatmap cellheight parameter
#'   - $fontsize_adj: adjusted font size
#'

#' @noRd
.calculate_adaptive_cellsizes <- function(n_cols_mat, n_rows_mat, width_frac = 1,
    cellwidth = 0, cellheight = 0, fontsize = 18) {
    # Base sizes (in pixels, for ~1200px wide plots)
    base_cellwidth <- 35
    base_cellheight <- 29

    if (cellwidth > 0 && cellheight > 0) {
        # Use explicit sizes provided
        final_cellwidth <- cellwidth
        final_cellheight <- cellheight
    } else {
        # Adaptive sizing
        if (width_frac < 1) {
            # Half-width: ~32% of plot width (1200px × 0.65 × 0.5)
            available_width_px <- 1200 * 0.65 * width_frac - 40 - 30
            cellwidth_calc <- available_width_px/max(1, n_cols_mat)
            final_cellwidth <- max(15, cellwidth_calc)  # Minimum 15px
        } else {
            # Full-width
            scale_factor_width <- if (n_cols_mat > 8)
                2.3 else 2.6
            final_cellwidth <- base_cellwidth * scale_factor_width
        }

        # Height scaling based on number of rows
        scale_factor_height <- max(0.7, 1.15 - n_rows_mat * 0.03)
        final_cellheight <- base_cellheight * scale_factor_height
    }

    final_cellwidth <- final_cellwidth * 0.95
    final_cellheight <- final_cellheight * 0.591

    list(cellwidth = final_cellwidth, cellheight = final_cellheight, fontsize_adj = fontsize *
        0.7)
}

# ============================================================================
# SECTION 5: PHEATMAP CREATION WRAPPER
# ============================================================================

#' Create Pheatmap Grob with Standard Settings
#'
#' Wrapper around pheatmap::pheatmap() that applies standard TSENAT
#' plotting settings (colors, sizing, clustering, etc.) and returns
#' a grob object for grid-based layout.
#'
#' @param matrix_data Numeric matrix to plot (rows by columns)
#' @param title Character: title for heatmap
#' @param cellw Numeric: cell width for pheatmap
#' @param cellh Numeric: cell height for pheatmap
#' @param fontsize Numeric: font size in points
#' @param cluster_rows Logical: cluster rows with dendrogram
#' @param color_palette Character vector: color palette (optional)
#'
#' @return Pheatmap grob object (from pheatmap::pheatmap)
#'

#' @noRd
.create_pheatmap_grob <- function(matrix_data, title = "", cellw = 35, cellh = 29,
    fontsize = 18, cluster_rows = FALSE, color_palette = NULL) {
    if (is.null(color_palette)) {
        color_palette <- (grDevices::colorRampPalette(c("#4575B4", "#FFFFFF", "#D73027")))(70)
    }

    pheatmap::pheatmap(matrix_data, main = title, cluster_rows = cluster_rows, cluster_cols = (ncol(matrix_data) >
        1), display_numbers = FALSE, na_col = "lightgray", border_color = "black",
        color = color_palette, cellwidth = cellw, cellheight = cellh, fontsize = fontsize,
        fontsize_row = fontsize, fontsize_col = fontsize, fontsize_number = fontsize *
            0.8, margins = c(8, 10), show_rownames = TRUE, show_colnames = TRUE,
        silent = TRUE)
}

# ============================================================================
# SECTION 6: GRID RENDERING HELPERS
# ============================================================================

#' Setup Grid for Rendering Heatmaps
#'
#' Opens PNG file if needed, initializes grid page, and adds title/subtitle.
#' Pushes main viewport for grid layout of heatmaps.
#'
#' @param n_layout_rows Integer: number of rows in grid layout
#' @param output_file Character: path to output PNG file (NULL = use active
#' device)
#' @param png_width Numeric: width in inches for PNG
#' @param png_height Numeric: height in inches for PNG
#' @param title Character: main title for figure
#' @param subtitle Character: subtitle for figure
#'
#' @return Invisibly returns NULL. Side effects: opens PNG, initializes grid.
#'

#' @noRd
.plot_grid_setup <- function(n_layout_rows, output_file = NULL, png_width = 12, png_height = 8,
    title = "Heatmap Analysis", subtitle = "") {
    # Open PNG if specified
    if (!is.null(output_file)) {
        # Create parent directories if they don't exist
        output_dir <- dirname(output_file)
        if (!dir.exists(output_dir)) {
            dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
        }
        grDevices::png(output_file, width = png_width, height = png_height, units = "in",
            res = 100)
    }

    # Initialize grid page
    grid::grid.newpage()

    # Calculate title sizes
    title_fontsize <- 16 * (1 + 0.15 * n_layout_rows)
    subtitle_fontsize <- 12 * (1 + 0.15 * n_layout_rows)

    # Add title
    grid::grid.text(title, x = 0.5, y = 0.97, just = "top", gp = grid::gpar(fontsize = title_fontsize,
        fontface = "bold"))

    # Add subtitle if provided
    if (nzchar(subtitle)) {
        grid::grid.text(subtitle, x = 0.5, y = 0.94, just = "top", gp = grid::gpar(fontsize = subtitle_fontsize,
            fontface = "italic", col = "gray40"))
    }

    # Push main viewport for grid layout
    n_grid_rows <- n_layout_rows * 2 - 1
    row_heights <- rep(c(1, 0.15), n_layout_rows)[seq_len(n_grid_rows)]

    grid::pushViewport(grid::viewport(x = 0.5, y = 0.47, width = 0.99, height = 0.85,
        layout = grid::grid.layout(n_grid_rows, 3, heights = grid::unit(row_heights,
            "null"), widths = c(1, 0.08, 1), respect = FALSE)))

    invisible(NULL)
}

#' Render Heatmaps into Grid Layout
#'
#' Positions and draws heatmap grob objects into the grid layout
#' established by .plot_grid_setup().
#'
#' @param heatmap_plots List of pheatmap grob objects
#' @param gene_layout List of layout positions (from .plot_adaptive_layout)
#' @param layout_ncol Integer: columns in fixed layout (or NULL for adaptive)
#'
#' @return Invisibly returns NULL. Side effect: draws heatmaps in grid.
#'

#' @noRd
.render_heatmaps_to_grid <- function(heatmap_plots, gene_layout, layout_ncol = 2) {
    if (is.null(gene_layout)) {
        # Simple rendering: assume layout_ncol columns per row
        n_cols <- as.integer(layout_ncol)
        for (i in seq_along(heatmap_plots)) {
            if (is.null(heatmap_plots[[i]]))
                next

            grid_row <- ((i - 1)%/%n_cols) * 2 + 1
            col_pos <- ((i - 1)%%n_cols) + 1

            if (n_cols == 1) {
                grid_col_start <- 1
                grid_col_end <- 3
            } else if (col_pos == 1) {
                grid_col_start <- 1
                grid_col_end <- 1
            } else {
                grid_col_start <- 3
                grid_col_end <- 3
            }

            grid::pushViewport(grid::viewport(layout.pos.row = grid_row, layout.pos.col = grid_col_start:grid_col_end))
            grid::grid.draw(heatmap_plots[[i]])
            grid::popViewport()
        }
    } else {
        # Use gene_layout positions
        for (i in seq_along(heatmap_plots)) {
            if (is.null(heatmap_plots[[i]]) || is.null(gene_layout[[i]]))
                next

            layout_info <- gene_layout[[i]]
            grid_row <- layout_info$row * 2 - 1

            if (layout_info$width == 1) {
                grid_col_start <- 1
                grid_col_end <- 3
            } else if (layout_info$col == 1) {
                grid_col_start <- 1
                grid_col_end <- 1
            } else {
                grid_col_start <- 3
                grid_col_end <- 3
            }

            grid::pushViewport(grid::viewport(layout.pos.row = grid_row, layout.pos.col = grid_col_start:grid_col_end))
            grid::grid.draw(heatmap_plots[[i]])
            grid::popViewport()
        }
    }

    invisible(NULL)
}

#' Finalize Grid Rendering and Output
#'
#' Pops viewport and closes PNG device if one was opened.
#'
#' @param output_file Character: output file path (if PNG was opened)
#' @param verbose Logical: print message when file saved
#'
#' @return Invisibly returns output_file (or NULL if rendered to device)
#'

#' @noRd
.plot_grid_finalize <- function(output_file = NULL, verbose = FALSE) {
    grid::popViewport()

    if (!is.null(output_file)) {
        grDevices::dev.off()
        if (verbose) {
            message("Heatmap saved to: ", output_file)
        }
        return(invisible(output_file))
    }

    invisible(NULL)
}

# ============================================================================
# SECTION 7: DATA PREPARATION HELPERS
# ============================================================================

#' Prepare Heatmap Data from Multi-Q Results
#'
#' Extracts delta_influence values across all q-values for a single gene
#' and formats into a matrix suitable for pheatmap visualization.
#' Rows = q-values, Columns = transcripts.
#'
#' @param switching_results Multi-q result list from jackknife_isoform_switching
#' @param gene_id Character: gene ID to extract
#' @param q_result_keys Character vector: q-value result keys (q_0_01, etc.)
#' @param cap_outliers_pctl Numeric: percentile at which to cap outliers
#' (default 95)
#'
#' @return Numeric matrix with q-values as rows and transcripts as columns.
#'   NAs indicate missing/infinite values. Rownames are q-value labels.
#'

#' @noRd
.heatmap_prepare_multiq_data <- function(switching_results, gene_id, q_result_keys,
    cap_outliers_pctl = 0.95) {
    heatmap_data <- NULL

    for (q_key in q_result_keys) {
        if (is.null(switching_results[[q_key]]) || is.null(switching_results[[q_key]]$results_per_gene) ||
            !gene_id %in% names(switching_results[[q_key]]$results_per_gene)) {
            next
        }

        gene_res <- switching_results[[q_key]]$results_per_gene[[gene_id]]
        if (is.null(gene_res$delta_influence))
            next

        # Extract q value from key (q_0_01 -> 0.01)
        q_str_cleaned <- gsub("_", ".", gsub("^q_", "", q_key))
        q_num <- as.numeric(q_str_cleaned)
        col_name <- paste0("q_", sprintf("%.2f", q_num))

        delta_vals <- as.numeric(gene_res$delta_influence)
        delta_vals[!is.finite(delta_vals)] <- NA

        if (is.null(heatmap_data)) {
            heatmap_data <- data.frame(transcript = as.character(gene_res$transcript_ids),
                stringsAsFactors = FALSE)
        }

        # Ensure row alignment
        n_rows <- nrow(heatmap_data)
        if (length(delta_vals) < n_rows) {
            delta_vals <- c(delta_vals, rep(NA_real_, n_rows - length(delta_vals)))
        } else if (length(delta_vals) > n_rows) {
            delta_vals <- delta_vals[seq_len(n_rows)]
        }

        heatmap_data[[col_name]] <- as.numeric(delta_vals)
    }

    if (is.null(heatmap_data) || nrow(heatmap_data) == 0) {
        return(NULL)
    }

    # Convert to matrix
    heatmap_matrix <- as.matrix(heatmap_data[, -1, drop = FALSE])
    rownames(heatmap_matrix) <- heatmap_data$transcript
    colnames(heatmap_matrix) <- colnames(heatmap_data)[-1]

    # Remove all-NA rows
    valid_rows <- rowSums(!is.na(heatmap_matrix)) > 0
    heatmap_matrix <- heatmap_matrix[valid_rows, , drop = FALSE]

    if (nrow(heatmap_matrix) == 0 || ncol(heatmap_matrix) == 0) {
        return(NULL)
    }

    # Cap outliers
    finite_vals <- heatmap_matrix[is.finite(heatmap_matrix)]
    if (length(finite_vals) > 0) {
        cap_val <- as.numeric(quantile(abs(as.numeric(finite_vals)), probs = cap_outliers_pctl))
        abs_hm <- abs(heatmap_matrix)
        mask <- which(is.finite(abs_hm) & abs_hm > cap_val)
        if (length(mask) > 0) {
            heatmap_matrix[mask] <- sign(heatmap_matrix[mask]) * cap_val
        }
    }

    # Transpose: q-values as rows, transcripts as columns
    t(heatmap_matrix)
}

#' Prepare Heatmap Data from Condition Samples
#'
#' Aggregates transcript counts across samples by condition, applies
#' log-normalization, and returns matrix for pheatmap.
#' Rows = conditions, Columns = transcripts.
#'
#' @param counts Numeric matrix: transcript counts (transcripts × samples)
#' @param gene_transcripts Integer vector: indices of transcripts for this gene
#' @param conditions Character vector: condition labels for each sample
#' @param metric Character: aggregation metric (mean, median, variance, iqr)
#' @param pseudocount Numeric: added before log transform (default 1e-6)
#'
#' @return Numeric matrix with conditions as rows and transcripts as columns,
#' log2-transformed. Rownames are condition labels, colnames are transcript
#' IDs.
#'

#' @noRd
.heatmap_prepare_condition_data <- function(counts, gene_transcripts, conditions,
    metric = "median", pseudocount = 1e-06) {
    unique_conditions <- unique(conditions)

    # Get expression for this gene's transcripts
    gene_counts <- counts[gene_transcripts, , drop = FALSE]

    # Aggregate by condition
    condition_matrix <- matrix(0, nrow = nrow(gene_counts), ncol = length(unique_conditions))
    rownames(condition_matrix) <- rownames(gene_counts)
    colnames(condition_matrix) <- unique_conditions

    for (cond in unique_conditions) {
        cond_mask <- conditions == cond
        if (sum(cond_mask) == 0)
            next

        if (metric == "mean") {
            condition_matrix[, cond] <- rowMeans(gene_counts[, cond_mask, drop = FALSE])
        } else if (metric == "median") {
            condition_matrix[, cond] <- apply(gene_counts[, cond_mask, drop = FALSE],
                1, median)
        } else if (metric == "variance") {
            condition_matrix[, cond] <- apply(gene_counts[, cond_mask, drop = FALSE],
                1, var)
        } else if (metric == "iqr") {
            condition_matrix[, cond] <- apply(gene_counts[, cond_mask, drop = FALSE],
                1, IQR)
        }
    }

    # Log-normalize
    condition_matrix_log <- log2(condition_matrix + pseudocount)

    # Transpose: conditions as rows, transcripts as columns
    t(condition_matrix_log)
}

# ============================================================================
# END OF HELPER FUNCTIONS
# ============================================================================

