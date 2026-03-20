#' Plot GAM q-curves from TSENATAnalysis object
#'
#' S4 wrapper that accepts a TSENATAnalysis object and generates GAM q-curve plots
#' for top genes identified by LM interaction analysis. Automatically extracts required
#' data from object slots.
#'
#' @param analysis \code{TSENATAnalysis} object containing:
#'   \itemize{
#'     \item \code{@se}: SummarizedExperiment with diversity values
#'     \item \code{@lm_results$lm_interaction}: Results from calculate_lm_interaction()
#'     \item \code{@config}: Configuration including condition_col if available
#'   }
#'
#' @param n_top \code{integer}. Number of top genes (by adjusted p-value) to plot 
#'   (default: 6). Only used if genes = NULL.
#'
#' @param genes \code{character} vector. Optional specific gene names to plot. 
#'   If provided, these genes are plotted directly regardless of significance.
#'   If NULL (default), top n_top significant genes are selected.
#'
#' @param condition_col \code{character}. Column name in colData(se) specifying 
#'   group assignments for samples. If NULL, attempts to auto-detect from @config
#'   (looks for "condition_col" or "condition"). If still NULL, defaults to "sample_type".
#'
#' @param sig_alpha \code{numeric}. Significance threshold for adjusted p-values 
#'   (default: 0.05). Only used if genes = NULL; filters lm_res to significant 
#'   genes before selecting top n.
#'
#' @param assay_name \code{character}. Name of the assay in se to extract 
#'   (default: "diversity").
#'
#' @param ... Additional arguments passed to \code{\link{plot_lm_interaction_gam}}.
#'
#' @return A single \code{ggplot} object with all selected genes arranged in a 
#'   grid layout. Can be saved with \code{ggplot2::ggsave()}.
#'
#' @details
#' This wrapper automatically:
#' 1. Extracts SummarizedExperiment from \code{@se} slot
#' 2. Extracts LM results from \code{@lm_results$lm_interaction} slot
#' 3. Detects condition_col from \code{@config} or uses default
#' 4. Calls \code{plot_lm_interaction_gam()} with extracted parameters
#'
#' **Parameter Resolution (condition_col):**
#' \enumerate{
#'   \item Explicit \code{condition_col} parameter (highest priority)
#'   \item \code{@config$condition_col} if available
#'   \item \code{@config$condition} if available  
#'   \item Default: "sample_type"
#' }
#'
#' @seealso \code{\link{plot_lm_interaction_gam}} for the underlying implementation,
#' \code{\link{calculate_lm_interaction_s4}} for running LM analysis on TSENATAnalysis.
#'
#' @examples
#' \dontrun{
#'   # Create TSENATAnalysis with diversity and LM results
#'   analysis <- TSENATAnalysis(se = se)
#'   analysis <- calculate_diversity_s4(analysis, q = seq(0.5, 2, by = 0.5))
#'   analysis <- calculate_lm_interaction_s4(analysis, condition_col = "condition")
#'   
#'   # Plot GAM curves for top genes
#'   plot <- plot_lm_interaction_gam_s4(
#'     analysis,
#'     n_top = 6,
#'     condition_col = "condition"
#'   )
#'   
#'   # Or plot specific genes
#'   plot <- plot_lm_interaction_gam_s4(
#'     analysis,
#'     genes = c("gene_1", "gene_3", "gene_5"),
#'     condition_col = "condition"
#'   )
#' }
#'
#' @export
plot_lm_interaction_gam_s4 <- function(
  analysis,
  n_top = 6,
  genes = NULL,
  condition_col = NULL,
  sig_alpha = 0.05,
  assay_name = "diversity",
  ...
) {
  # =========================================================================
  # INPUT VALIDATION
  # =========================================================================
  if (!is(analysis, "TSENATAnalysis")) {
    stop("'analysis' must be a TSENATAnalysis object", call. = FALSE)
  }
  
  # Check that LM results exist
  if (is.null(analysis@lm_results) || is.null(analysis@lm_results$lm_interaction)) {
    stop("[plot_lm_interaction_gam_s4] No LM interaction results found in @lm_results$lm_interaction. ",
         "Run calculate_lm_interaction_s4() first.",
         call. = FALSE)
  }
  
  lm_res <- analysis@lm_results$lm_interaction
  
  if (!is.data.frame(lm_res)) {
    stop("[plot_lm_interaction_gam_s4] @lm_results$lm_interaction must be a data.frame",
         call. = FALSE)
  }
  
  # Check that diversity results exist (needed for SE reconstruction)
  if (length(analysis@diversity_results) == 0) {
    stop("[plot_lm_interaction_gam_s4] No diversity results found in @diversity_results. ",
         "Run calculate_diversity_s4() first.",
         call. = FALSE)
  }
  
  # =========================================================================
  # AUTO-DETECT condition_col FROM @config IF NOT PROVIDED
  # =========================================================================
  if (is.null(condition_col)) {
    cd_cols <- colnames(colData(analysis@se))
    
    # Try Priority 1: @config$condition_col (validate it exists)
    if ("condition_col" %in% names(analysis@config)) {
      candidate <- analysis@config$condition_col
      if (candidate %in% cd_cols) {
        condition_col <- candidate
      }
    }
    
    # Try Priority 2: @config$condition (validate it exists)
    if (is.null(condition_col) && "condition" %in% names(analysis@config)) {
      candidate <- analysis@config$condition
      if (candidate %in% cd_cols) {
        condition_col <- candidate
      }
    }
    
    # Try Priority 3: Check for common column names in actual colData
    if (is.null(condition_col)) {
      if ("condition" %in% cd_cols) {
        condition_col <- "condition"
      }
      else if ("sample_type" %in% cd_cols) {
        condition_col <- "sample_type"
      }
      else {
        # Use first column as fallback
        if (length(cd_cols) > 0) {
          condition_col <- cd_cols[1]
        } else {
          stop("[plot_lm_interaction_gam_s4] No columns found in colData(se). ",
               "Cannot auto-detect condition_col.",
               call. = FALSE)
        }
      }
    }
  }
  
  # Validate that condition_col exists in colData
  if (!(condition_col %in% colnames(colData(analysis@se)))) {
    stop("[plot_lm_interaction_gam_s4] Specified condition_col='", condition_col, 
         "' not found in colData. Available columns: ",
         paste(colnames(colData(analysis@se)), collapse = ", "),
         call. = FALSE)
  }
  
  # =========================================================================
  # RECONSTRUCT COMBINED DIVERSITY SE FOR PLOTTING
  # (Same approach as in calculate_lm_interaction_s4)
  # =========================================================================
  # Extract q-values from diversity_results keys
  q_keys <- names(analysis@diversity_results)
  q_computed <- as.numeric(sub("^q_", "", q_keys))
  
  # Reconstruct combined diversity SE with all q-values
  diversity_combined <- tryCatch({
    verbose <- if ("verbose" %in% names(analysis@config)) {
      analysis@config$verbose
    } else {
      FALSE
    }
    
    calculate_diversity(
      x = analysis@se,
      q = sort(q_computed),
      norm = TRUE,
      verbose = verbose,
      bootstrap = FALSE
    )
  }, error = function(e) {
    stop(paste0("[plot_lm_interaction_gam_s4] Failed to reconstruct diversity SE:\n",
                conditionMessage(e)), call. = FALSE)
  })
  
  # =========================================================================
  # EXTRACT model_data FROM STORED RESULTS
  # =========================================================================
  model_data <- NULL
  if ("lm_interaction_model_data" %in% names(analysis@lm_results)) {
    model_data <- analysis@lm_results$lm_interaction_model_data
  }
  
  # =========================================================================
  # CALL plot_lm_interaction_gam WITH RECONSTRUCTED DIVERSITY SE
  # =========================================================================
  result <- tryCatch({
    plot_lm_interaction_gam(
      se = diversity_combined,
      lm_res = lm_res,
      condition_col = condition_col,
      n_top = n_top,
      genes = genes,
      sig_alpha = sig_alpha,
      assay_name = assay_name,
      model_data = model_data,
      ...
    )
  }, error = function(e) {
    stop(paste0("[plot_lm_interaction_gam_s4] Error in plot generation:\n", 
                conditionMessage(e)),
         call. = FALSE)
  })
  
  # =========================================================================
  # RETURN PLOT
  # =========================================================================
  # Track that plotting occurred
  if (is.list(analysis@metadata)) {
    analysis@metadata$function_calls <- c(
      analysis@metadata$function_calls,
      paste0("plot_lm_interaction_gam_s4[n_top=", n_top, ", condition_col=", 
             condition_col, "]")
    )
  }
  
  # Return the plot object directly (not the analysis object)
  result
}
