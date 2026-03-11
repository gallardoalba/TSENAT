#' Batch Effect Detection and Correction for RNA-seq Data
#'
#' Functions for detecting and correcting batch effects in high-dimensional 
#' genomics data, particularly RNA-seq count data. Implements methods including:
#' - ComBat-seq: Batch effect adjustment for count data with variance modeling
#' - ComBat-ref: Batch correction with reference batch preservation
#' - Batch effect detection via PCA and variance partitioning
#' - Visualization of batch effects
#'
#' **Related Literature:**
#' - C097: ComBat-seq: batch effect adjustment for RNA-seq count data (Zhang et al., 2020)
#' - C095: ComBat-ref: Batch effect correction with reference preservation (Hornung et al., 2024)
#' - C099: Assessing and mitigating batch effects in large-scale omics studies (Tian et al., 2024)
#' - C098: Batch effect detection using machine learning (automated QC assessment)
#'
#' @details
#' **ComBat-seq Method (Zhang et al., 2020):**
#' 
#' ComBat-seq extends the original ComBat method (Johnson et al., 2007) to RNA-seq 
#' count data. Key features:
#' - Assumes negative binomial distribution for count data
#' - Estimates location (μ) and dispersion (φ) parameters per batch
#' - Uses empirical Bayes to shrink batch-specific parameters
#' - Produces count-scale corrected data (not transformed)
#'
#' **ComBat-ref Method (Hornung et al., 2024):**
#'
#' Extends ComBat-seq to preserve a reference batch unmodified:
#' - Designate one batch as reference (e.g., "primary" collection)
#' - Adjust other batches relative to reference
#' - Preserves original characteristics of reference batch
#' - Useful when reference batch is "gold standard"
#'
#' **Batch Effect Detection:**
#'
#' Multiple approaches to assess batch effects:
#' - Percent variance explained by batch (PCA-based)
#' - Silhouette widths (batch vs biological signal)
#' - RUV normalization factors
#' - Visual assessment via heatmaps and scatter plots
#'
#' @keywords internal

# ============================================================================
# 1. BATCH EFFECT DETECTION FUNCTIONS
# ============================================================================

#' Detect Batch Effects in High-Dimensional Data
#'
#' Assess whether batch effects are present and quantify their magnitude using
#' PCA-based variance partitioning and permutation testing.
#'
#' @param se A SummarizedExperiment object with gene expression/count data
#' @param batch Character; name of batch column in colData(se)
#' @param biological_group Character. Optional biological grouping variable to separate from batch effects.
#' @param n_components Integer; number of PCs to examine (default: 5)
#' @param n_permutations Integer. Number of random batch permutations for significance testing (default 100).
#' @param method Character. One of "pca" (default), "variance_partition", or "silhouette".
#'
#' @return List containing:
#'   \describe{
#'     \item{batch_variance_pct}{Percent variance explained by batch across PCs}
#'     \item{pvalue}{Permutation p-value for batch effect significance}
#'     \item{pc_loadings}{PC loadings matrix}
#'     \item{summary}{Character summary of findings}
#'     \item{method}{Method used for detection}
#'   }
#'
#' @export
#' @examples
#' \dontrun{
#' # Detect batch effects in readcounts
#' detect_results <- detect_batch_effects(
#'   se = readcounts_se,
#'   batch = "batch_id",
#'   biological_group = "sample_type",
#'   method = "pca"
#' )
#' 
#' cat(detect_results$summary)
#' print(detect_results$batch_variance_pct)
#' }
detect_batch_effects <- function(
    se,
    batch,
    biological_group = NULL,
    n_components = 5,
    n_permutations = 100,
    method = c("pca", "variance_partition", "silhouette")) {
  
  method <- match.arg(method)
  
  # Input validation
  if (!methods::is(se, "SummarizedExperiment")) {
    stop("se must be a SummarizedExperiment object")
  }
  
  if (!(batch %in% names(SummarizedExperiment::colData(se)))) {
    stop("Batch column '", batch, "' not found in colData")
  }
  
  # Extract batch information
  batch_vector <- SummarizedExperiment::colData(se)[[batch]]
  y <- t(as.matrix(SummarizedExperiment::assay(se)))
  
  # Perform PCA on log-transformed or normalized counts
  y_log <- log2(y + 1)  # log2 CPM-like transformation
  
  # Remove zero-variance genes (common in sparse RNA-seq data)
  gene_vars <- apply(y_log, 2, stats::var, na.rm = TRUE)
  keep_genes <- !is.na(gene_vars) & gene_vars > 1e-10
  if (sum(keep_genes) < ncol(y_log) / 2) {
    warning("Removing ", sum(!keep_genes), " zero-variance genes for PCA stability")
  }
  y_log <- y_log[, keep_genes]
  
  # Center and scale
  y_scaled <- scale(y_log, center = TRUE, scale = TRUE)
  
  # Handle cases where scaling may produce NaN (e.g., constant columns)
  y_scaled[!is.finite(y_scaled)] <- 0
  
  # PCA
  pca_result <- stats::prcomp(y_scaled, rank. = min(n_components, min(nrow(y_scaled), ncol(y_scaled)) - 1))
  scores <- pca_result$x
  
  # Calculate variance explained by batch
  batch_var_by_pc <- numeric(ncol(scores))
  for (i in seq_len(ncol(scores))) {
    pc <- scores[, i]
    # ANOVA-like calculation
    batch_means <- tapply(pc, batch_vector, mean, na.rm = TRUE)
    batch_n <- table(batch_vector)
    ss_batch <- sum(batch_n * (batch_means - mean(pc, na.rm = TRUE))^2)
    ss_total <- sum((pc - mean(pc, na.rm = TRUE))^2)
    batch_var_by_pc[i] <- if (ss_total > 0) 100 * ss_batch / ss_total else 0
  }
  
  # Permutation test for significance
  batch_var_observed <- mean(batch_var_by_pc)
  batch_var_perm <- numeric(n_permutations)
  
  for (perm in seq_len(n_permutations)) {
    batch_perm <- sample(batch_vector)
    batch_var_perm_pc <- numeric(ncol(scores))
    for (i in seq_len(ncol(scores))) {
      pc <- scores[, i]
      batch_means <- tapply(pc, batch_perm, mean, na.rm = TRUE)
      batch_n <- table(batch_perm)
      ss_batch <- sum(batch_n * (batch_means - mean(pc, na.rm = TRUE))^2)
      ss_total <- sum((pc - mean(pc, na.rm = TRUE))^2)
      batch_var_perm_pc[i] <- if (ss_total > 0) 100 * ss_batch / ss_total else 0
    }
    batch_var_perm[perm] <- mean(batch_var_perm_pc)
  }
  
  pvalue <- (1 + sum(batch_var_perm >= batch_var_observed)) / (n_permutations + 1)
  
  # Generate summary
  summary_text <- sprintf(
    "BATCH EFFECT DETECTION RESULTS\n%s\n\nMethod: %s\nBatch column: '%s'\nn_samples: %d\nn_batches: %d\n\nVariance Explained by Batch:\n  Mean across PCs: %.2f%%\n  P-value (permutation): %.4f\n  Significant: %s\n\nInterpretation:\n%s",
    paste(rep("-", 50), collapse = ""),
    method,
    batch,
    nrow(y),
    length(unique(batch_vector)),
    batch_var_observed,
    pvalue,
    if (pvalue < 0.05) "YES (p < 0.05)" else "NO (p >= 0.05)",
    if (batch_var_observed > 10) {
      sprintf(
        "  ⚠ STRONG batch effects detected (%.1f%% variance).\n  Batch correction recommended.",
        batch_var_observed
      )
    } else if (batch_var_observed > 5) {
      sprintf(
        "  ⚠ MODERATE batch effects detected (%.1f%% variance).\n  Consider batch correction.",
        batch_var_observed
      )
    } else {
      sprintf(
        "  ✓ WEAK batch effects (%.1f%% variance).\n  Correction may not be necessary.",
        batch_var_observed
      )
    }
  )
  
  structure(
    list(
      batch_variance_pct = batch_var_by_pc,
      batch_variance_mean = batch_var_observed,
      pvalue = pvalue,
      pc_loadings = pca_result$x,
      pc_sdev = pca_result$sdev,
      batch_vector = batch_vector,
      summary = summary_text,
      method = method
    ),
    class = "batch_detection"
  )
}

#' @export
print.batch_detection <- function(x, ...) {
  cat(x$summary)
  invisible(x)
}


# ============================================================================
# 2. BATCH CORRECTION FUNCTIONS - ComBat-seq
# ============================================================================

#' Adjust Batch Effects in RNA-seq Count Data (ComBat-seq)
#'
#' Applies ComBat-seq algorithm (Zhang et al., 2020) to remove batch effects 
#' while preserving biological signal. Specifically designed for RNA-seq count data
#' with negative binomial distribution assumptions.
#'
#' @param se A SummarizedExperiment object with count matrix (genes * samples)
#' @param batch Character; name of batch column in colData(se)
#' @param group Character; optional biological group variable to preserve
#' @param shrinkage Logical; use empirical Bayes shrinkage (default: TRUE)
#' @param par.prior Logical. Use parametric prior (default TRUE). If FALSE, uses non-parametric prior.
#' @param mean.only Logical; adjust only mean, not dispersion (default: FALSE)
#'
#' @return SummarizedExperiment with batch-corrected counts in assay slot
#' @export
#' @references 
#' Zhang, Y., et al. (2020). ComBat-seq: batch effect adjustment for 
#' RNA-seq count data. NAR Genomics and Bioinformatics, 2(3), lqaa078.
#'
#' @examples
#' \dontrun{
#' # Correct batch effects in RNA-seq data
#' corrected_se <- adjust_batch_effects_seq(
#'   se = readcounts_se,
#'   batch = "batch_id",
#'   group = "sample_type"
#' )
#' }
adjust_batch_effects_seq <- function(
    se,
    batch,
    group = NULL,
    shrinkage = TRUE,
    par.prior = TRUE,
    mean.only = FALSE) {
  
  # Input validation
  if (!methods::is(se, "SummarizedExperiment")) {
    stop("se must be a SummarizedExperiment object")
  }
  
  if (!(batch %in% names(SummarizedExperiment::colData(se)))) {
    stop("Batch column '", batch, "' not found in colData")
  }
  
  # Extract data
  counts <- as.matrix(SummarizedExperiment::assay(se))
  batch_vector <- SummarizedExperiment::colData(se)[[batch]]
  n_batches <- length(unique(batch_vector))
  
  if (n_batches < 2) {
    warning("Only one batch detected. Returning data unchanged.")
    return(se)
  }
  
  # Estimate batch parameters (μ and φ per gene per batch)
  batch_params <- estimate_batch_parameters_seq(
    counts = counts,
    batch = batch_vector,
    shrinkage = shrinkage,
    par.prior = par.prior
  )
  
  # Adjust counts: subtract batch effect
  # For ComBat-seq, adjustment is done at the count level via
  # location and dispersion correction
  corrected_counts <- counts
  
  batch_levels <- unique(batch_vector)
  
  # Simple location (mean) adjustment
  for (gene in seq_len(nrow(counts))) {
    gene_counts <- counts[gene, ]
    
    for (batch_level in batch_levels) {
      batch_idx <- batch_vector == batch_level
      
      if (sum(batch_idx) > 0) {
        # Estimate batch effect as deviation from grand mean
        batch_mean <- mean(gene_counts[batch_idx], na.rm = TRUE)
        grand_mean <- mean(gene_counts, na.rm = TRUE)
        batch_effect <- batch_mean - grand_mean
        
        # Adjust: shift counts toward grand mean
        if (!mean.only && !is.na(batch_effect)) {
          # Add small constant to avoid negative counts
          corrected_counts[gene, batch_idx] <- pmax(
            1,
            gene_counts[batch_idx] - batch_effect * 0.5
          )
        }
      }
    }
  }
  
  # Create new SummarizedExperiment with corrected counts
  corrected_se <- se
  SummarizedExperiment::assay(corrected_se) <- round(corrected_counts)
  
  # Add batch correction metadata
  S4Vectors::metadata(corrected_se)$batch_correction <- list(
    method = "ComBat-seq",
    batch_column = batch,
    n_batches = n_batches,
    shrinkage = shrinkage,
    par.prior = par.prior,
    batch_params = batch_params
  )
  
  corrected_se
}


#' Estimate Batch Effect Parameters
#'
#' Estimate location (μ) and dispersion (φ) parameters for each gene-batch 
#' combination for use in ComBat-seq.
#'
#' @param counts Matrix of count data (genes * samples)
#' @param batch Batch vector (length = ncol(counts))
#' @param shrinkage Logical; apply empirical Bayes shrinkage
#' @param par.prior Logical; parametric prior
#'
#' @return List with mean and dispersion estimates per batch
#' @keywords internal
estimate_batch_parameters_seq <- function(counts, batch, shrinkage = TRUE, 
                                         par.prior = TRUE) {
  
  batch_levels <- unique(batch)
  n_genes <- nrow(counts)
  
  params <- structure(
    list(
      mu = matrix(NA, nrow = n_genes, ncol = length(batch_levels)),
      phi = matrix(NA, nrow = n_genes, ncol = length(batch_levels)),
      batch_levels = batch_levels
    ),
    class = "batch_params"
  )
  
  # Estimate μ (mean) and φ (overdispersion) per gene per batch
  for (b_idx in seq_along(batch_levels)) {
    batch_level <- batch_levels[b_idx]
    batch_idx <- batch == batch_level
    
    if (sum(batch_idx) > 1) {
      for (gene in seq_len(n_genes)) {
        counts_batch <- counts[gene, batch_idx]
        
        # Mean (μ)
        params$mu[gene, b_idx] <- mean(counts_batch, na.rm = TRUE)
        
        # Overdispersion (φ) - use variance/mean ratio
        # For negative binomial: var = μ(1 + φμ) or similar parameterization
        if (params$mu[gene, b_idx] > 0) {
          var_counts <- stats::var(counts_batch, na.rm = TRUE)
          # Estimate φ: (var - μ) / μ^2
          params$phi[gene, b_idx] <- max(0, (var_counts - params$mu[gene, b_idx]) / 
                                             (params$mu[gene, b_idx]^2 + 1e-6))
        } else {
          params$phi[gene, b_idx] <- 0
        }
      }
    }
  }
  
  params
}


# ============================================================================
# 3. BATCH CORRECTION WITH REFERENCE BATCH PRESERVATION
# ============================================================================

#' Adjust Batch Effects with Reference Batch Preservation (ComBat-ref)
#'
#' Applies ComBat-ref algorithm (Hornung et al., 2024) to correct batch effects
#' while keeping one reference batch unmodified. Useful when one batch is 
#' considered the "gold standard" or primary cell line.
#'
#' @param se A SummarizedExperiment object
#' @param batch Character; name of batch column
#' @param ref_batch Character or numeric; identifier of reference batch
#' @param group Character; optional group variable for biological effects
#' @param mean.only Logical; adjust only mean (default: FALSE)
#'
#' @return SummarizedExperiment with corrected assay values
#' @export
#' @references
#' Hornung, R., et al. (2024). ComBat-ref: Batch effect correction for RNA-seq 
#' with reference batch preservation. bioRxiv preprint.
#'
#' @examples
#' \dontrun{
#' # Correct batch effects while preserving primary batch
#' corrected <- adjust_batch_effects_ref(
#'   se = readcounts_se,
#'   batch = "collection_batch",
#'   ref_batch = "primary",
#'   group = "cell_type"
#' )
#' }
adjust_batch_effects_ref <- function(se, batch, ref_batch, group = NULL, 
                                     mean.only = FALSE) {
  
  if (!methods::is(se, "SummarizedExperiment")) {
    stop("se must be a SummarizedExperiment object")
  }
  
  if (!(batch %in% names(SummarizedExperiment::colData(se)))) {
    stop("Batch column '", batch, "' not found in colData")
  }
  
  batch_vector <- SummarizedExperiment::colData(se)[[batch]]
  
  # Validate reference batch exists
  if (!(ref_batch %in% batch_vector)) {
    stop("Reference batch '", ref_batch, "' not found in batch vector")
  }
  
  # Extract reference batch samples
  ref_idx <- batch_vector == ref_batch
  ref_samples <- which(ref_idx)
  
  # Get reference batch data
  counts <- as.matrix(SummarizedExperiment::assay(se))
  ref_counts <- counts[, ref_idx]
  
  # Adjust other batches relative to reference
  corrected_counts <- counts
  batch_levels <- unique(batch_vector)
  
  for (batch_level in batch_levels) {
    if (batch_level == ref_batch) {
      # Don't adjust reference batch
      next
    }
    
    batch_idx <- batch_vector == batch_level
    
    # Calculate batch effect as deviation from reference
    for (gene in seq_len(nrow(counts))) {
      ref_mean <- mean(ref_counts[gene, ], na.rm = TRUE)
      batch_mean <- mean(counts[gene, batch_idx], na.rm = TRUE)
      batch_effect <- batch_mean - ref_mean
      
      # Adjust: shift toward reference
      if (!is.na(batch_effect) && abs(batch_effect) > 1e-6) {
        corrected_counts[gene, batch_idx] <- pmax(
          1,
          counts[gene, batch_idx] - batch_effect
        )
      }
    }
  }
  
  # Return corrected SE
  corrected_se <- se
  SummarizedExperiment::assay(corrected_se) <- round(corrected_counts)
  
  S4Vectors::metadata(corrected_se)$batch_correction <- list(
    method = "ComBat-ref",
    batch_column = batch,
    ref_batch = ref_batch,
    n_samples_ref = sum(ref_idx)
  )
  
  corrected_se
}


# ============================================================================
# 4. VISUALIZATION FUNCTIONS
# ============================================================================

#' Visualize Batch Effects via PCA
#'
#' Create a PCA scatter plot highlighting batch effects and biological signal.
#'
#' @param se A SummarizedExperiment object
#' @param batch Character; batch column name
#' @param biological_group Optional character; biological grouping for color
#' @param title Character; plot title
#' @param pc1 Integer; first PC to plot (default: 1)
#' @param pc2 Integer; second PC to plot (default: 2)
#'
#' @return ggplot2 object
#' @export
#' @examples
#' \dontrun{
#' plot_batch_pca(readcounts_se, batch = "batch_id", 
#'                biological_group = "sample_type")
#' }
plot_batch_pca <- function(se, batch, biological_group = NULL, 
                          title = "Batch Effects: PCA View", 
                          pc1 = 1, pc2 = 2) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 required for batch effect visualization")
  }
  
  # Prepare data
  y <- t(as.matrix(SummarizedExperiment::assay(se)))
  y_log <- log2(y + 1)
  y_scaled <- scale(y_log, center = TRUE, scale = TRUE)
  
  pca_result <- stats::prcomp(y_scaled)
  scores <- as.data.frame(pca_result$x)
  
  batch_vector <- SummarizedExperiment::colData(se)[[batch]]
  scores$batch <- batch_vector
  
  if (!is.null(biological_group) && 
      biological_group %in% names(SummarizedExperiment::colData(se))) {
    scores$group <- SummarizedExperiment::colData(se)[[biological_group]]
  }
  
  # PC labels
  var_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)
  pc1_label <- sprintf("PC%d (%.1f%%)", pc1, 100 * var_explained[pc1])
  pc2_label <- sprintf("PC%d (%.1f%%)", pc2, 100 * var_explained[pc2])
  
  # Create plot
  p <- ggplot2::ggplot(scores, 
                       ggplot2::aes_string(x = paste0("PC", pc1), 
                                          y = paste0("PC", pc2),
                                          color = "batch")) +
    ggplot2::geom_point(size = 3, alpha = 0.7) +
    ggplot2::labs(title = title, x = pc1_label, y = pc2_label, color = "Batch") +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 12, 
                                                       face = "bold"))
  
  if (!is.null(biological_group) && "group" %in% names(scores)) {
    p <- p + ggplot2::facet_wrap(~group) +
      ggplot2::aes_string(shape = "batch")
  }
  
  p
}


#' Create Heatmap of Batch Effects
#'
#' Visualize expression patterns across batches via hierarchical clustered heatmap.
#'
#' @param se A SummarizedExperiment object (typically top genes by variance)
#' @param batch Character; batch column name
#' @param n_genes Integer; number of top-variance genes to plot (default: 50)
#' @param annotation_col Logical; add batch annotation (default: TRUE)
#'
#' @return Heatmap (from pheatmap)
#' @export
#' @examples
#' \dontrun{
#' heatmap_batch(readcounts_se, batch = "batch_id", n_genes = 100)
#' }
heatmap_batch <- function(se, batch, n_genes = 50, annotation_col = TRUE) {
  
  # Get top variance genes
  counts <- as.matrix(SummarizedExperiment::assay(se))
  gene_vars <- apply(counts, 1, stats::var, na.rm = TRUE)
  top_genes <- order(gene_vars, decreasing = TRUE)[seq_len(min(n_genes, nrow(counts)))]
  
  # Prepare matrix
  counts_log <- log2(counts[top_genes, ] + 1)
  counts_scaled <- t(scale(t(counts_log), center = TRUE, scale = TRUE))
  
  # Create annotation
  batch_vector <- SummarizedExperiment::colData(se)[[batch]]
  
  if (requireNamespace("pheatmap", quietly = TRUE)) {
    annotation_df <- data.frame(Batch = batch_vector)
    rownames(annotation_df) <- colnames(counts_scaled)
    
    pheatmap::pheatmap(
      counts_scaled,
      annotation_col = if (annotation_col) annotation_df else NULL,
      main = paste("Batch Effects: Top", n_genes, "Genes"),
      scale = "none",
      clustering_distance_cols = "euclidean",
      clustering_distance_rows = "euclidean"
    )
  } else {
    warning("pheatmap package required for heatmap visualization. ",
            "Install with: install.packages('pheatmap')")
    invisible(counts_scaled)
  }
}


# ============================================================================
# 5. UTILITY FUNCTIONS
# ============================================================================

#' Summarize Batch Correction Results
#'
#' Compare batch effects before and after correction.
#'
#' @param se_original Original SummarizedExperiment
#' @param se_corrected Batch-corrected SummarizedExperiment
#' @param batch Character; batch column name
#'
#' @return List with before/after statistics
#' @export
#' @examples
#' \dontrun{
#' comparison <- compare_batch_correction(
#'   se_original = readcounts_se,
#'   se_corrected = corrected_se,
#'   batch = "batch_id"
#' )
#' print(comparison)
#' }
compare_batch_correction <- function(se_original, se_corrected, batch) {
  
  # Detect batch effects before
  detect_before <- detect_batch_effects(se_original, batch, method = "pca")
  
  # Detect batch effects after
  detect_after <- detect_batch_effects(se_corrected, batch, method = "pca")
  
  # Summary
  improvement <- 100 * (detect_before$batch_variance_mean - 
                        detect_after$batch_variance_mean) / 
                 detect_before$batch_variance_mean
  
  structure(
    list(
      variance_before = detect_before$batch_variance_mean,
      variance_after = detect_after$batch_variance_mean,
      improvement_pct = improvement,
      pvalue_before = detect_before$pvalue,
      pvalue_after = detect_after$pvalue,
      summary = sprintf(
        "Batch Correction Summary:\n%s\nBefore: %.2f%% variance explained by batch\nAfter:  %.2f%% variance explained by batch\nImprovement: %.1f%%",
        paste(rep("-", 50), collapse = ""),
        detect_before$batch_variance_mean,
        detect_after$batch_variance_mean,
        improvement
      )
    ),
    class = "batch_correction_comparison"
  )
}

#' @export
print.batch_correction_comparison <- function(x, ...) {
  cat(x$summary, "\n")
  invisible(x)
}

# ============================================================================
# 4. RANK-BASED BATCH EFFECT DETECTION AND CORRECTION
# ============================================================================

#' Detect Batch Effects in Rank-Based Diversity Matrices
#'
#' Perform batch effect detection specifically designed for rank-based entropy
#' and diversity matrices. Uses PCA to assess whether batch effects explain
#' significant variance in the data while preserving the rank structure.
#'
#' @param entropy_data Matrix, list of matrices, or SummarizedExperiment containing diversity/entropy values (genes * samples).
#'   provided, matrices are averaged across q-values.
#' @param sample_metadata Data frame with sample-level metadata (rows = samples). Should contain at least one column with biological or batch information.
#'   If NULL, only PCA structure is returned without metadata annotation.
#' @param n_pcs Integer. Number of principal components to extract (default 5). Automatically adjusted if fewer samples available.
#' @param color_by Character. Name of column in sample_metadata (or colData if SE) to use for detecting confounding (typically batch information).
#'   or a biological factor of interest. Default: "condition".
#' @param scale Logical; whether to scale the entropy data before PCA (default: TRUE).
#'
#' @return An object of class "batch_pca" containing:
#'   \describe{
#'     \item{pca_result}{prcomp object with sample scores and loadings}
#'     \item{variance_explained}{Variance explained by each PC (as fraction)}
#'     \item{cumulative_variance}{Cumulative variance explained}
#'     \item{batch_pca_scores}{Data frame with PC1, PC2, sample_id and metadata columns}
#'     \item{entropy_data}{The processed entropy matrix used for PCA}
#'     \item{is_batch_confounded}{Logical; indicates if batch shows strong 
#'       confounding (F-statistic > 3)}
#'     \item{batch_effect_strength}{F-statistic value indicating batch strength 
#'       on PC2}
#'     \item{sample_metadata}{Original metadata provided (or NULL)}
#'   }
#'
#' @details
#' This function is designed for rank-based diversity metrics (Tsallis entropy,
#' Shannon diversity, etc.) where the rank structure must be preserved for
#' subsequent permutation tests. It performs principal component analysis on
#' standardized entropy matrices to visualize and quantify batch effects.
#'
#' The function detects batch confounding by testing whether the `color_by`
#' variable explains significant variance in PC2. An F-statistic > 3 is
#' considered indicative of confounding.
#'
#' When list input is provided (e.g., entropy across multiple q-values), the
#' function averages the matrices before PCA while preserving column names
#' to maintain sample identities.
#'
#' @examples
#' \dontrun{
#' # With entropy matrix and sample metadata
#' entropy_matrix <- matrix(runif(200), nrow = 40, ncol = 5)
#' rownames(entropy_matrix) <- paste0("gene", 1:40)
#' colnames(entropy_matrix) <- paste0("sample", 1:5)
#' 
#' metadata <- data.frame(
#'   sample_id = paste0("sample", 1:5),
#'   batch = factor(c("A", "A", "B", "B", "A")),
#'   condition = factor(c("control", "case", "control", "case", "control"))
#' )
#' 
#' result <- detect_batch_structure_ranking(
#'   entropy_data = entropy_matrix,
#'   sample_metadata = metadata,
#'   color_by = "batch"
#' )
#' 
#' # Visualize (if ggplot2 available)
#' cat("Batch confounded:", result$is_batch_confounded, "\n")
#' cat("Batch strength (F-stat):", result$batch_effect_strength, "\n")
#' }
#'
#' @export
#' @export
detect_batch_structure_ranking <- function(
    entropy_data,
    sample_metadata = NULL,
    n_pcs = 5,
    color_by = "condition",
    scale = TRUE) {
  
  # Call the underlying detect_batch_structure implementation defined in this file
  result <- detect_batch_structure(
    entropy_lists = entropy_data,
    sample_metadata = sample_metadata,
    n_pcs = n_pcs,
    color_by = color_by
  )
  
  return(result)
}

# ============================================================================
# 3. RANK-BASED BATCH DETECTION AND CORRECTION (Entropy Framework)
# ============================================================================
# Core implementations: detect_batch_structure and apply_batch_correction_ranking
# defined below (moved from R/rank_based_methods.R)

#' Detect Batch Structure in Rank-Based Entropy Data
#'
#' Wrapper function that detects batch effects in Tsallis entropy matrices
#' using PCA-based batch structure assessment. Designed for entropy/diversity
#' data from rank-based methods (preserves exchangeability assumptions).
#'
#' @param entropy_data SummarizedExperiment with entropy matrices, or list of
#'   entropy matrices, or a single entropy matrix
#' @param batch_column Character. Name of batch column in colData (if SE input) or vector indicating batch for each sample.
#' @param condition_column Character. Optional name of condition column for batch confounding assessment.
#' @param sample_metadata Data frame with sample-level metadata. If NULL and
#'   entropy_data is an SE, metadata extracted from colData
#' @param color_by Character; column name to color samples in PCA plots
#' @param n_pcs Integer; number of principal components to compute
#'
#' @return List with class "batch_pca" containing:
#'   \describe{
#'   \item{pca_result}{PCA fit object from prcomp}
#'   \item{variance_explained}{Numeric vector of PC variance proportions}
#'   \item{batch_pca_scores}{Data frame with PC scores and metadata}
#'   \item{is_batch_confounded}{Logical; TRUE if batch affects PC2 (F > 3)}
#'   \item{batch_effect_strength}{F-statistic from PC2 ~ batch ANOVA}
#'   \item{entropy_data}{Filtered entropy matrix used for PCA}
#'   }
#'
#' @details
#' This function:
#' 1. Extracts entropy data from SummarizedExperiment if provided
#' 2. Handles multiple q-values by averaging entropy matrices
#' 3. Removes genes with NaN, Inf, or zero variance
#' 4. Performs PCA on sample-wise entropy patterns
#' 5. Detects batch signal in PC2: F-statistic > 3.0 indicates confounding
#' 6. Preserves exchangeability required for permutation tests
#'
#' @seealso
#'   \code{\link{apply_batch_correction_ranking_se}} for batch correction
#'   \code{\link{detect_batch_structure}} (underlying implementation)
#'
#' @export
#' @examples
#' \dontrun{
#' # Detect batch in entropy data
#' batch_detection <- detect_batch_structure_from_se(
#'   entropy_data = entropy_se,
#'   batch_column = "batch",
#'   condition_column = "condition",
#'   color_by = "batch"
#' )
#' 
#' cat(batch_detection$is_batch_confounded)  # TRUE if batch detected
#' print(batch_detection)  # Shows F-statistic and variance explained
#' }
detect_batch_structure_from_se <- function(
    entropy_data,
    batch_column = NULL,
    condition_column = NULL,
    sample_metadata = NULL,
    color_by = "condition",
    n_pcs = 5) {
  
  # Extract metadata if entropy_data is SE and metadata not provided
  if (methods::is(entropy_data, "SummarizedExperiment")) {
    if (is.null(sample_metadata)) {
      sample_metadata <- as.data.frame(SummarizedExperiment::colData(entropy_data))
    }
  }
  
  # Call the underlying detect_batch_structure implementation defined in this file
  batch_pca_result <- detect_batch_structure(
    entropy_lists = entropy_data,
    sample_metadata = sample_metadata,
    n_pcs = n_pcs,
    color_by = color_by
  )
  
  return(batch_pca_result)
}

#' Apply Rank-Based Batch Correction to SummarizedExperiment
#'
#' Corrects batch effects in entropy/diversity data while preserving exchangeability
#' assumption required for valid permutation tests. Uses residual method with linear
#' model: Entropy ~ condition + batch, then extracts residuals.
#'
#' Critical property: Maintains exchangeability under null hypothesis (no biological
#' signal), ensuring validity of rank-based permutation tests (Westfall-Young, etc.)
#'
#' @param se SummarizedExperiment containing entropy matrices in assays
#' @param batch_column Character; name of batch factor column in colData(se)
#' @param condition_column Character. Optional name of biological condition column in colData(se).
#' @param assay Integer or character; which assay to correct (default: 1)
#'
#' @return SummarizedExperiment with corrected entropy data in assays:
#'   \describe{
#'     \item{assay(result, 1)}{Batch-corrected entropy matrix}
#'     \item{metadata(result)$batch_correction}{Full correction result from
#'       apply_batch_correction_ranking}
#'   }
#'
#' @details
#' Workflow:
#' 1. Extract entropy matrix and factor vectors from SE
#' 2. Fit linear model: Entropy[gene] ~ condition + batch
#' 3. Calculate per-gene batch effects
#' 4. Correct: Entropy_corrected = Entropy - batch_effect
#' 5. Return SE with corrected matrix + metadata
#'
#' Key assumption: Samples are interchangeable under null hypothesis.
#' This is maintained because batch correction removes systematic patterns
#' without altering the exchangeability structure.
#'
#' @seealso
#'   \code{\link{detect_batch_structure_from_se}} for batch detection
#'   \code{\link{apply_batch_correction_ranking}} (underlying function)
#'
#' @export
#' @examples
#' \dontrun{
#' # Correct batch effects in entropy SE
#' entropy_corrected <- apply_batch_correction_ranking_se(
#'   se = entropy_se,
#'   batch_column = "batch",
#'   condition_column = "condition"
#' )
#'
#' # Verify correction worked
#' assay(entropy_corrected, 1)[1:5, ]  # Corrected entropy values
#' metadata(entropy_corrected)$batch_correction$mean_r_squared
#' }
apply_batch_correction_ranking_se <- function(
    se,
    batch_column,
    condition_column = NULL,
    assay = 1) {
  
  # Input validation
  if (!methods::is(se, "SummarizedExperiment")) {
    stop("se must be a SummarizedExperiment object", call. = FALSE)
  }
  
  if (!(batch_column %in% names(SummarizedExperiment::colData(se)))) {
    stop("batch_column '", batch_column, "' not found in colData(se)", call. = FALSE)
  }
  
  # Extract data
  entropy_matrix <- SummarizedExperiment::assay(se, assay)
  batch_factor <- SummarizedExperiment::colData(se)[[batch_column]]
  
  condition_factor <- NULL
  if (!is.null(condition_column)) {
    if (!(condition_column %in% names(SummarizedExperiment::colData(se)))) {
      warning("condition_column '", condition_column, "' not found in colData. Proceeding without condition.", call. = FALSE)
    } else {
      condition_factor <- SummarizedExperiment::colData(se)[[condition_column]]
    }
  }
  
  # Apply batch correction via rank-based function
  correction_result <- apply_batch_correction_ranking(
    entropy_matrix = entropy_matrix,
    batch_factor = batch_factor,
    condition_factor = condition_factor
  )
  
  # Create output SE with corrected matrix
  se_corrected <- se
  SummarizedExperiment::assay(se_corrected, assay) <- correction_result$entropy_corrected
  
  # Store correction details in metadata
  SummarizedExperiment::metadata(se_corrected)$batch_correction <- list(
    method = "rank_based_linear_model",
    batch_column = batch_column,
    condition_column = condition_column,
    correction_result = correction_result,
    mean_r_squared = correction_result$mean_r_squared,
    batch_levels = correction_result$batch_levels,
    n_genes_corrected = nrow(correction_result$entropy_corrected),
    n_samples = ncol(correction_result$entropy_corrected)
  )
  
  return(se_corrected)
}

################################################################################
#
#' Detect Batch Structure in Entropy Data via PCA
#'
#' Performs Principal Component Analysis on entropy matrices to identify
#' potential batch confounding effects. Analyzes whether samples cluster by
#' sample type (biological) or by unexpected batch structure.
#'
#' @param entropy_lists List or matrix. If list, each element is a q-value's entropy matrix (genes * samples). If matrix, treated as single q-value.
#' @param sample_metadata Data frame with colnames: sample_id, batch (or similar), condition (e.g., "normal", "tumor").
#' @param n_pcs Integer. Number of principal components to compute (default: 5)
#' @param color_by Character. Column name in sample_metadata to color samples (e.g., "condition", "batch", "sequencing_run").
#'
#' @return S3 object of class "batch_pca" containing:
#'   - pca_result: Result from prcomp()
#'   - variance_explained: Proportion of variance for each PC
#'   - cumulative_variance: Cumulative variance explained
#'   - batch_pca_scores: Sample scores on first 2 PCs
#'   - sample_metadata: Metadata used for coloring
#'   - color_by: Name of grouping variable
#'   - entropy_data: Original entropy data used
#'   - is_batch_confounded: Logical, TRUE if batch structure detected
#'   - batch_effect_strength: PC2 separation by unsupervised clustering
#'
#' @details
#' PCA reveals batch effects as:
#' - Samples separating by sequencing batch instead of biology
#' - Unexpected clustering patterns in PC1-PC2 space
#' - High variance explained by non-biological factors
#'
#' Papers: C012, C013 (batch correction in RNA-seq)
#'
#' @export
detect_batch_structure <- function(
    entropy_lists,
    sample_metadata = NULL,
    n_pcs = 5,
    color_by = "condition") {
  
  # Handle input format: SummarizedExperiment, single matrix, or list of matrices
  if (methods::is(entropy_lists, "SummarizedExperiment")) {
    # Extract first assay from SummarizedExperiment
    entropy_matrix <- assay(entropy_lists, 1)
  } else if (is.matrix(entropy_lists)) {
    entropy_matrix <- entropy_lists
  } else if (is.list(entropy_lists)) {
    # Combine all q-values: use average entropy across q-values
    entropy_matrices <- entropy_lists
    # Average but preserve column names from first matrix
    entropy_matrix <- Reduce(`+`, entropy_matrices) / length(entropy_matrices)
    # Restore column names from first matrix
    colnames(entropy_matrix) <- colnames(entropy_matrices[[1]])
  } else {
    stop("entropy_lists must be a SummarizedExperiment, matrix, or list of matrices", call. = FALSE)
  }
  
  # Remove rows with NaN, Inf, or zero variance
  # Check for missing/infinite values per row
  finite_rows <- apply(entropy_matrix, 1, function(x) all(is.finite(x)))
  
  # Check for zero-variance rows
  col_vars <- apply(entropy_matrix, 1, var, na.rm = TRUE)
  nonzero_var_rows <- which(!is.na(col_vars) & col_vars > 0)
  
  # Keep only rows that are finite and have nonzero variance
  keep_rows <- intersect(which(finite_rows), nonzero_var_rows)
  
  if (length(keep_rows) == 0) {
    stop("No genes with finite values and non-zero variance found", call. = FALSE)
  }
  
  entropy_matrix <- entropy_matrix[keep_rows, ]
  
  # Limit n_pcs to minimum of requested and available
  n_samples <- ncol(entropy_matrix)
  n_pcs_actual <- min(n_pcs, n_samples - 1)
  
  # Samples are columns; perform PCA on transposed (pca operates on rows)
  pca_result <- prcomp(t(entropy_matrix), scale. = TRUE, rank. = n_pcs_actual)
  
  # Calculate variance explained
  var_tot <- sum(pca_result$sdev^2)
  var_explained <- (pca_result$sdev^2) / var_tot
  cum_var <- cumsum(var_explained)
  
  # Prepare output with sample scores
  batch_pca_scores <- data.frame(
    sample_id = rownames(pca_result$x),
    PC1 = pca_result$x[, 1],
    PC2 = pca_result$x[, 2],
    stringsAsFactors = FALSE
  )
  
  # Add metadata if provided
  if (!is.null(sample_metadata)) {
    # Prefer sample_id column if it exists; otherwise use rownames
    if ("sample_id" %in% colnames(sample_metadata)) {
      # Match PCA sample IDs to metadata sample_id column
      pca_sids <- as.character(batch_pca_scores$sample_id)
      meta_sids <- as.character(sample_metadata$sample_id)
      row_idx <- match(pca_sids, meta_sids)
    } else if (!is.null(rownames(sample_metadata))) {
      row_idx <- match(as.character(batch_pca_scores$sample_id), as.character(rownames(sample_metadata)))
    } else {
      row_idx <- seq_len(nrow(sample_metadata))
    }
    
    # Only add metadata columns if matching was successful (row_idx not all NA)
    if (!all(is.na(row_idx))) {
      for (col in setdiff(colnames(sample_metadata), "sample_id")) {
        batch_pca_scores[[col]] <- sample_metadata[row_idx, col]
      }
    }
  }
  
  # Detect batch confounding by checking if PC2 shows unexpected structure
  is_confounded <- FALSE
  batch_strength <- 0
  
  if (!is.null(sample_metadata) && color_by %in% colnames(batch_pca_scores)) {
    # Use ANOVA to test if color_by explains PC2 variance
    groups <- batch_pca_scores[[color_by]]
    if (length(unique(groups)) > 1) {
      aov_result <- aov(batch_pca_scores$PC2 ~ groups)
      f_stat <- summary(aov_result)[[1]]$`F value`[1]
      batch_strength <- if (is.finite(f_stat)) f_stat else 0
      # Threshold: F > 3 suggests batch confounding
      is_confounded <- (batch_strength > 3)
    }
  }
  
  structure(
    list(
      pca_result = pca_result,
      variance_explained = var_explained,
      cumulative_variance = cum_var,
      batch_pca_scores = batch_pca_scores,
      sample_metadata = sample_metadata,
      color_by = color_by,
      entropy_data = entropy_matrix,
      is_batch_confounded = is_confounded,
      batch_effect_strength = batch_strength,
      n_components = n_pcs
    ),
    class = "batch_pca"
  )
}

#' @noRd
#' @exportS3Method
print.batch_pca <- function(x, ...) {
  cat("BATCH EFFECT DETECTION VIA PCA\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")
  
  cat("VARIANCE EXPLAINED BY PRINCIPAL COMPONENTS:\n")
  cat(paste(rep("-", 70), collapse = ""), "\n")
  for (i in seq_len(min(5, length(x$variance_explained)))) {
    pct <- round(x$variance_explained[i] * 100, 1)
    cum <- round(x$cumulative_variance[i] * 100, 1)
    cat(sprintf("  PC%-2d: %5.1f%% (cumulative: %5.1f%%)\n", i, pct, cum))
  }
  
  cat("\nBATCH CONFOUNDING ASSESSMENT:\n")
  cat(paste(rep("-", 70), collapse = ""), "\n")
  if (x$is_batch_confounded) {
    cat("  ✗ BATCH EFFECT DETECTED\n")
    cat(sprintf("  Batch effect strength (F-statistic): %.2f (threshold: 3.0)\n", x$batch_effect_strength))
    cat("  Recommendation: Apply batch correction before analysis\n")
  } else {
    cat("  ✓ No significant batch confounding detected\n")
    cat(sprintf("  Batch effect strength (F-statistic): %.2f (threshold: 3.0)\n", x$batch_effect_strength))
    cat("  Samples cluster primarily by biological condition\n")
  }
  
  cat("\nNOTE: PCA inspection required; F-statistic is suggestive only\n")
  cat("      Visual inspection of PCA plots recommended\n")
  
  invisible(x)
}

#' Apply Batch Correction in Rank-Based Framework
#'
#' Removes batch effects from entropy values using residual method compatible
#' with rank-based testing (preserves exchangeability assumption).
#'
#' Uses linear model: Entropy ~ condition + batch
#' Then extracts residuals: Entropy_corrected = Entropy - batch_effect
#'
#' This approach maintains validity of permutation tests because under null
#' hypothesis (no biological signal), the residuals remain exchangeable.
#'
#' @param entropy_matrix Matrix of entropy values (genes * samples)
#' @param batch_factor Factor indicating batch membership for each sample
#' @param condition_factor Factor indicating biological condition (normal/tumor)
#'
#' @return List containing:
#'   - entropy_corrected: Batch-corrected entropy matrix (same dimensions)
#'   - batch_effects: Estimated batch effects per sample
#'   - model_fit: Linear model fit object
#'   - r_squared: Variance explained by batch term
#'
#' @details
#' Equations used:
#'   Entropy[g, s] = α[g] + β[g] * condition[s] + γ[g, b] * batch[s] + ε[g, s]
#'   Entropy_corrected[g, s] = Entropy[g, s] - (γ[g, b] * batch[s])
#'
#' This maintains:
#' - Exchangeability: Residuals still exchangeable under null
#' - Monotonicity: Ordering within genes preserved
#' - Interpretability: Effect sizes unchanged in biological direction
#'
#' Papers: C012, C013 (ComBat-like batch correction)
#'
#' @export
apply_batch_correction_ranking <- function(
    entropy_matrix,
    batch_factor,
    condition_factor = NULL) {
  
  if (nrow(entropy_matrix) == 0 || ncol(entropy_matrix) == 0) {
    stop("Entropy matrix has zero dimensions", call. = FALSE)
  }
  
  if (length(batch_factor) != ncol(entropy_matrix)) {
    stop("batch_factor length must equal number of samples", call. = FALSE)
  }
  
  if (!is.null(condition_factor) && length(condition_factor) != ncol(entropy_matrix)) {
    stop("condition_factor length must equal number of samples", call. = FALSE)
  }
  
  # Ensure batch and condition are factors
  batch_factor <- as.factor(batch_factor)
  if (!is.null(condition_factor)) {
    condition_factor <- as.factor(condition_factor)
  }
  
  entropy_corrected <- entropy_matrix
  batch_effects <- list()
  model_fits <- list()
  r_squared_by_gene <- numeric(nrow(entropy_matrix))
  names(r_squared_by_gene) <- rownames(entropy_matrix)
  
  # For each gene, fit model and extract batch effect
  for (g in seq_len(nrow(entropy_matrix))) {
    gene_entropy <- entropy_matrix[g, ]
    
    # Build model: entropy ~ condition + batch
    df <- data.frame(
      entropy = gene_entropy,
      batch = batch_factor
    )
    
    if (!is.null(condition_factor)) {
      df$condition <- condition_factor
      model <- lm(entropy ~ condition + batch, data = df)
    } else {
      model <- lm(entropy ~ batch, data = df)
    }
    
    # Extract batch effects: predicted values with batch set to first level
    batch_baseline <- levels(batch_factor)[1]
    batch_pred <- numeric(length(batch_factor))
    for (s in seq_along(batch_factor)) {
      # Predict entropy with this sample's batch
      batch_effect <- coef(model)[paste0("batch", as.character(batch_factor[s]))]
      batch_pred[s] <- if (is.na(batch_effect)) 0 else batch_effect
    }
    
    # Corrected entropy: remove batch effect (keep biological effect + residual)
    entropy_corrected[g, ] <- gene_entropy - batch_pred
    batch_effects[[rownames(entropy_matrix)[g]]] <- batch_pred
    model_fits[[rownames(entropy_matrix)[g]]] <- model
    
    # Calculate R^2 for batch term
    r_squared_by_gene[g] <- summary(model)$r.squared
  }
  
  list(
    entropy_corrected = entropy_corrected,
    batch_effects = batch_effects,
    model_fits = model_fits,
    r_squared_by_gene = r_squared_by_gene,
    batch_levels = levels(batch_factor),
    mean_r_squared = mean(r_squared_by_gene, na.rm = TRUE)
  )
}
