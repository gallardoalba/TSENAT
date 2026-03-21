#!/usr/bin/env R

# =========================================================================
# PRIVATE HELPER FUNCTIONS (must come before main roxygen block)
# =========================================================================

# -------------------------------------------------------------------------
# PRIVATE HELPER: Auto-detect group column and control group
# -------------------------------------------------------------------------

#' Auto-Detect Group Column and Control Group from SummarizedExperiment colData
#'
#' Automatically searches for a column containing group/condition assignments
#' and identifies the reference (control) group.
#'
#' **Detection Strategy for group_col:**
#' Searches colData in precedence order for common grouping column names:
#' 1. "group" (TSENAT default)
#' 2. "condition", "treatment", "sample_type" (common alternatives)
#'
#' **Detection Strategy for control_group:**
#' Once group column is found, identifies the reference group by:
#' 1. First looking for common control/reference names: "Normal", "Control", 
#'    "WT" (wild-type), "Reference", "Baseline", "wt"
#' 2. If no match, selects the unique group value with smallest sample count
#'    (typically the control/reference in case-control designs)
#' 3. If still no match, uses the first alphabetically sorted group name
#'
#' @param se SummarizedExperiment object with sample metadata in colData
#'
#' @return List with elements:
#'   \describe{
#'     \item{group_col}{Name of the colData column used for grouping, or NA_character_ if none detected}
#'     \item{control_group}{Name of the control/reference group, or NA_character_ if none detected}
#'     \item{groups}{Character vector of all unique groups found}
#'     \item{sample_counts}{Named integer vector of sample counts per group (names: group names)}
#'   }
#'
#' @keywords internal
#' @noRd
.auto_detect_groups <- function(se) {
  
  cd <- SummarizedExperiment::colData(se)
  cd_colnames <- colnames(cd)
  
  # Candidate column names (in priority order) for group/condition
  group_col_candidates <- c(
    "sample_type",      # TSENAT standard (created by map_metadata)
    "group",            # Alias for sample_type
    "condition",        # Common alternative
    "treatment",        # Experimental design
    "phenotype",        # Biological phenotype
    "batch",            # Last resort
    "category"          # Generic fallback
  )
  
  group_col <- NA_character_
  
  # Find first matching group column
  for (col_name in group_col_candidates) {
    if (col_name %in% cd_colnames) {
      group_col <- col_name
      break
    }
  }
  
  # If no group column found, return NAs
  if (is.na(group_col)) {
    return(list(
      group_col = NA_character_,
      control_group = NA_character_,
      groups = character(0),
      sample_counts = integer(0)
    ))
  }
  
  # Get group values and counts
  group_vec <- as.character(cd[[group_col]])
  unique_groups <- unique(group_vec)
  group_counts <- table(group_vec)
  
  # Candidate names for control/reference group (in priority order)
  control_candidates <- c(
    "Normal",           # TSENAT default
    "Control",          # Most common control label
    "WT",               # Wild-type (common in genetics)
    "wt",               # Lowercase variant
    "Reference",        # Explicit reference label
    "Baseline",         # Baseline condition
    "Wild-type",        # Full name variant
    "wild_type"         # Underscore variant
  )
  
  control_group <- NA_character_
  
  # Try to match control candidates
  for (control_name in control_candidates) {
    if (control_name %in% unique_groups) {
      control_group <- control_name
      break
    }
  }
  
  # If no match found, use heuristics:
  # 1. Select group with fewer samples (typical case-control design)
  # 2. Fallback to first alphabetically
  if (is.na(control_group)) {
    if (length(unique_groups) >= 2) {
      # Find group with minimum samples (typically control)
      min_samples_group <- names(group_counts)[which.min(group_counts)]
      control_group <- min_samples_group
    } else if (length(unique_groups) == 1) {
      control_group <- unique_groups[1]
    }
  }
  
  return(list(
    group_col = group_col,
    control_group = control_group,
    groups = unique_groups,
    sample_counts = as.vector(group_counts)
  ))
}


# -------------------------------------------------------------------------
# PRIVATE HELPER: Auto-detect paired samples from metadata
# -------------------------------------------------------------------------

#' Detect Paired Sample Structure from SummarizedExperiment colData
#'
#' Automatically searches for a column containing paired sample identifiers
#' (e.g., "pair_id", "paired_samples", "patient_id", "subject_id").
#' This enables pair-respecting bootstrap resampling in divergence calculations.
#'
#' **Detection Strategy:**
#' Searches colData in precedence order for common pairing column names:
#' 1. "paired_samples" (TSENAT default, matches readcounts metadata)
#' 2. "pair_id", "pair_samples", "subject_id", "patient_id" (common alternatives)
#'
#' Returns a mapping from sample names to pair identifiers, or NULL if no
#' pairing column is found. A valid pairing column has:
#' - Non-NA values for all samples
#' - At least 2 samples per pair
#' - Deterministic structure (e.g., all A's paired with another A sample nearby)
#'
#' **Database References (Papers validating auto-detection approach):**
#' - S102: "Experimental Control and Paired Design" - standardizes paired design annotation
#' - S107: "Related Sample Designs and Paired t-test" - validates paired structure detection
#'
#' @param se SummarizedExperiment object with sample metadata in colData
#'
#' @return List with elements:
#'   \describe{
#'     \item{pair_ids}{Character vector (names: sample names, values: pair identifiers)
#'       or NULL if no pairing detected}
#'     \item{column_name}{Name of the colData column used, or NA_character_ if none}
#'     \item{num_pairs}{Number of unique pairs (0 if none detected)}
#'     \item{samples_per_pair}{Vector of samples per pair (names: pair IDs, values: counts)}
#'   }
#'
#' @note Paired samples detected from any of: "paired_samples", "pair_id", "pair_samples",
#'   "subject_id", "patient_id". Returns NULL if none present or validation fails.
#'
#' @keywords internal
#' @noRd
.detect_pair_ids <- function(se) {
  
  cd <- SummarizedExperiment::colData(se)
  sample_names <- colnames(se)
  
  # Candidate column names (in priority order)
  candidate_cols <- c(
    "paired_samples",      # TSENAT default
    "pair_id",             # Common alternative
    "pair_samples",        # Variant
    "subject_id",          # Statistical defaults
    "patient_id"
  )
  
  for (col_name in candidate_cols) {
    if (col_name %in% colnames(cd)) {
      pair_col <- cd[[col_name]]
      
      # Validate: must be non-NA for all samples
      if (any(is.na(pair_col))) {
        next  # Skip if any NAs
      }
      
      # Valid pairing structure found
      pair_ids <- setNames(as.character(pair_col), sample_names)
      unique_pairs <- unique(pair_ids)
      samples_per_pair <- table(pair_ids)
      
      return(list(
        pair_ids = pair_ids,
        column_name = col_name,
        num_pairs = length(unique_pairs),
        samples_per_pair = samples_per_pair
      ))
    }
  }
  
  # No pairing detected
  return(list(
    pair_ids = NULL,
    column_name = NA_character_,
    num_pairs = 0,
    samples_per_pair = numeric(0)
  ))
}


#' Resample Data Respecting Paired Structure
#'
#' When resampling paired data, both members of a pair are selected or discarded
#' together. This preserves within-pair correlations critical for statistical validity
#' in matched designs (Efron & Tibshirani 1993).
#'
#' For k pairs:
#' - Draw k pair indices uniformly with replacement: pair_idx ~ U(1:k)
#' - For each drawn pair i, include both samples from pair i
#' - This maintains pairing structure across bootstrap replicates
#'
#' **Statistical Justification (Papers C016, S102-S109):**
#' - C016: Bootstrap for confidence intervals requires preserving data structure
#' - S102: "Paired Design" - paired resampling required for matched samples
#' - S107: "Related Sample Designs" - within-pair correlation invalidates independent resampling
#'
#' @param control_samples Vector of counts for control group
#' @param treatment_samples Vector of counts for treatment group
#' @param pair_ids Named character vector: names = sample names, values = pair IDs
#' @param group_col Character vector: group assignment (one per sample)
#' @param control_group Character: name of control group
#'
#' @return List with elements:
#'   \describe{
#'     \item{control_resampled}{Resampled control group counts}
#'     \item{treatment_resampled}{Resampled treatment group counts}
#'   }
#'
#' @keywords internal
#' @noRd
.resample_paired_data <- function(
    control_samples,
    treatment_samples,
    pair_ids,
    group_col,
    control_group) {
  
  # Map sample names to indices
  all_samples <- c(names(control_samples), names(treatment_samples))
  all_groups <- c(
    rep(control_group, length(control_samples)),
    rep(setdiff(unique(group_col), control_group), length(treatment_samples))
  )
  
  # Get unique pairs involved
  pairs_in_data <- unique(pair_ids[all_samples])
  num_pairs <- length(pairs_in_data)
  
  # Resample pairs with replacement
  resampled_pairs <- sample(seq_len(num_pairs), size = num_pairs, replace = TRUE)
  resampled_pair_ids <- pairs_in_data[resampled_pairs]
  
  # Collect samples for each resampled pair
  resampled_control <- numeric(0)
  resampled_treatment <- numeric(0)
  
  for (pair_id in resampled_pair_ids) {
    # Get both samples from this pair
    pair_samples <- names(pair_ids)[pair_ids == pair_id]
    
    for (sample_name in pair_samples) {
      if (sample_name %in% names(control_samples)) {
        resampled_control <- c(resampled_control, control_samples[sample_name])
      } else if (sample_name %in% names(treatment_samples)) {
        resampled_treatment <- c(resampled_treatment, treatment_samples[sample_name])
      }
    }
  }
  
  # BUGFIX #3: Validate balanced groups after paired resampling
  # Ensures control and treatment have equal sizes (required for divergence computation)
  if (length(resampled_control) != length(resampled_treatment)) {
    stop("Paired bootstrap produced unequal group sizes (",
         length(resampled_control), " control vs ", 
         length(resampled_treatment), " treatment). ",
         "Check for unbalanced or incomplete pairs in input data.")
  }
  
  return(list(
    control_resampled = resampled_control,
    treatment_resampled = resampled_treatment
  ))
}


# =========================================================================
# PRIVATE HELPER: Classify Q-Pattern
# =========================================================================

#' Classify a per-q divergence spectrum into biological pattern types
#'
#' When Tsallis divergence has been computed across multiple q values
#' for a gene, the resulting vector can be summarised by its trend across the
#' spectrum. This helper compares divergence in the rare-region (q < 1) vs
#' abundant-region (q >= 1).
#'
#' \describe{
#'   \item{RARE_DRIVEN}{Divergence higher at low q (q < 1);
#'     indicates changes driven by low-abundance isoforms.}
#'   \item{ABUNDANT_DRIVEN}{Divergence higher at high q (q >= 1);
#'     indicates shifts among the most abundant transcripts.}
#'   \item{BALANCED}{Similar divergence across rare and abundant regions.}
#' }
#'
#' If the input vector is too short, contains only NAs, or classification
#' cannot be performed, NA is returned.
#'
#' @param per_q_divs Named numeric vector of divergences. Names should be of form
#'   "q_0.01", "q_0.5", "q_1.0", etc.
#' @param ratio_threshold Numeric; ratio threshold for classification (default: 1.3).
#'   RARE_DRIVEN if rare_median / abundant_median > threshold.
#'
#' @return Character scalar: "RARE_DRIVEN", "ABUNDANT_DRIVEN", "BALANCED", or NA.
#'
#' @keywords internal
#' @noRd
.classify_q_pattern <- function(per_q_divs, ratio_threshold = 1.3) {
  # Input validation
  if (!is.numeric(per_q_divs) || length(per_q_divs) < 2) {
    return(NA_character_)
  }
  
  # Get names
  nm <- names(per_q_divs)
  if (is.null(nm) || any(is.na(nm))) {
    return(NA_character_)
  }
  
  # Extract q values from names: try "q_0.5", "q_0_5", etc.
  q_vals <- NA
  
  # Try format: "q_0.5" (standard with dot)
  if (all(grepl("^q_", nm))) {
    q_vals <- suppressWarnings(as.numeric(gsub("^q_", "", nm)))
  } else if (all(grepl("^q", nm))) {
    # Try other formats
    q_vals <- suppressWarnings(as.numeric(gsub("^q[_.]", "", nm)))
  }
  
  # If we still can't extract numeric q values, return NA
  if (any(is.na(q_vals))) {
    return(NA_character_)
  }
  
  # Check if all divergence values are NA
  if (all(is.na(per_q_divs))) {
    return(NA_character_)
  }
  
  # Need at least 2 non-NA pairs for comparison
  valid_pairs <- !is.na(per_q_divs)
  if (sum(valid_pairs) < 2) {
    return(NA_character_)
  }
  
  # Split into rare-region (q < 1) and abundant-region (q >= 1)
  rare_mask <- q_vals < 1
  abund_mask <- q_vals >= 1
  
  # Calculate median divergence in each region
  rare_div_median <- NA
  abund_div_median <- NA
  
  if (sum(rare_mask & valid_pairs) > 0) {
    rare_div_median <- median(per_q_divs[rare_mask & valid_pairs], na.rm = TRUE)
  }
  
  if (sum(abund_mask & valid_pairs) > 0) {
    abund_div_median <- median(per_q_divs[abund_mask & valid_pairs], na.rm = TRUE)
  }
  
  # If we have both regions, compare them
  if (!is.na(rare_div_median) && !is.na(abund_div_median) && abund_div_median > 0) {
    ratio <- rare_div_median / abund_div_median
    
    if (ratio > ratio_threshold) {
      return("RARE_DRIVEN")
    } else if (ratio < 1 / ratio_threshold) {
      return("ABUNDANT_DRIVEN")
    } else {
      return("BALANCED")
    }
  }
  
  # Fallback: if only one region available, can't classify
  return(NA_character_)
}


# =========================================================================
# MAIN FUNCTION: Calculate Bootstrap Divergence
# =========================================================================

#' Calculate Bootstrap Divergence Confidence Intervals Across Genes
#'
#' **NEW ARCHITECTURE: Transcript-level counts -> Gene-level aggregation -> Tsallis divergence**
#' 
#' Computes bootstrap confidence intervals for Tsallis divergence comparing
#' two groups across multiple genes. Automatically aggregates transcript-level counts
#' to gene-level (per Paper I033: gene-level analysis for information-theoretic diversity).
#' Supports both sequential and parallel computation, with optional support for paired sample designs.
#' 
#' Returns a SummarizedExperiment object containing:
#' - **assay**: genes * q matrix of divergence estimates (one per q value)
#' - **rowData**: gene metadata including per-q divergence estimates, CIs, pattern classification
#' - **colData**: one row per q value with q-specific metadata
#' - **metadata**: processing parameters and summary statistics
#'
#' **INPUT & OUTPUT ARCHITECTURE:**
#' ```
#' calculate_divergence(se, res=NULL, ...)  
#'   Input:  SummarizedExperiment (raw TRANSCRIPT-level counts)
#'           Each row is a transcript; rowData must have gene_names/gene_name column
#'   Step 1: Auto-aggregates transcripts -> genes via colSums
#'   Step 2: Computes divergence for each gene across q values
#'   Output: SummarizedExperiment with:
#'           - assay: genes * q_values matrix (divergence estimates)
#'           - rowData: gene_name, per_q_pattern, estimate_q*, lower_ci_q*, etc.
#'           - colData: one row per q value
#'           - metadata: parameters, timing, sample sizes
#' ```
#' Matches `calculate_diversity()` input/output pattern: transcript counts SE -> gene-level derivative SE
#'
#' **DESIGN PRINCIPLE - Transcript-to-Gene Aggregation:**
#' Following Paper I033 ("Application of information theoretical approaches to assess diversity 
#' in single-cell transcriptomics"), divergence analysis operates on GENE-LEVEL expression profiles.
#' When input is transcript-level data (typical RNA-seq output), this function automatically:
#'   1. Identifies all transcripts for each gene (via rowData gene_names column)
#'   2. Sums counts across transcripts for each gene
#'   3. Computes divergence on aggregated gene-level counts
#' This ensures statistical validity (one observation per gene per sample) and biological relevance.
#'
#' @param se SummarizedExperiment object with transcript-level counts
#'           (assay called "counts", rowData with gene identifier columns)
#' @param group_col Character; colData column for group membership (optional).
#'            If NULL, auto-detects in this order: "group", "condition", "treatment", 
#'            "sample_type". If no match found, an error is raised.
#'            (default: NULL, auto-detect)
#' @param control_group Character; reference group name (optional).
#'            If NULL, auto-detects by: (1) looking for "Normal", "Control", "WT", etc.,
#'            or (2) selecting the group with fewer samples (typical case-control),
#'            or (3) first alphabetically.
#'            (default: NULL, auto-detect)
#' @param q Tsallis parameter (scalar or vector) (default: 1)
#' @param paired Logical; if TRUE or if paired_samples column detected, uses paired sample design.
#'               With bootstrap=TRUE, automatically detects paired samples from metadata
#'               column names (searched in order: "paired_samples", "pair_id", "pair_samples",
#'               "subject_id", "patient_id") and applies pair-respecting bootstrap resampling
#'               to preserve within-pair correlations (Papers C016, S102-S109).
#'               (default: FALSE)
#' @param bootstrap Logical; if TRUE, computes bootstrap confidence intervals (~2-3 sec/gene).
#'                  If FALSE, computes point estimates only (~0.02-0.05 sec/gene).
#'                  (default: FALSE)
#' @param nboot Number of bootstrap replicates (default: 1000)
#'               Note: ignored if bootstrap=FALSE
#' @param ci Confidence level (default: 0.95)
#' @param method Bootstrap method: "percentile" or "bca" (default: "percentile")
#' @param log_base Logarithm base (default: exp(1), natural log)
#' @param norm Logical or character; normalization/standardization mode (default: TRUE).
#'        Backward compatible: TRUE = "range", FALSE = "none".
#'        Options:
#'        - "none": Raw divergence values, no standardization
#'        - "range": Range standardization [0,1] per q (classic approach)
#'        - "zscore": Z-score standardization per q: (D_q - mean) / sd
#'          Useful for cross-study comparison; results in mean=0, sd=1
#'        - "log_odds_ratio": Log ratio relative to theoretical maximum
#'          D_norm = log(D_q / D_max) where D_max depends on q-value
#'          Interpretation: 0 = theoretical max, <0 = below max
#'        - "relative_reference": Ratio to reference group (requires group_col)
#'          Interpretation: Reference = 1, >1 higher than reference
#' @param pseudocount Pseudocount for stability (default: 0.5)
#' @param nthreads Number of CPU threads for parallel processing (default: 1).
#'                 Set to > 1 to parallelize gene-level bootstrap computations.
#'                 nthreads=NULL auto-detects available cores minus 1.
#' @param progress Logical; show progress bar and timing (default: TRUE)
#' @param seed Random seed (optional; NULL for non-reproducible)
#'
#' @return SummarizedExperiment object with:
#'   **assays** (genes * q matrices):
#'     - divergence: Divergence point estimates for each gene
#'   
#'   **rowData** (data frame with one row per gene):
#'     - gene_name: Gene identifier
#'     - estimate: Point estimate of divergence
#'     - lower_ci: Lower CI bound (NA if bootstrap=FALSE)
#'     - upper_ci: Upper CI bound (NA if bootstrap=FALSE)
#'     - ci_width: Width of confidence interval (NA if bootstrap=FALSE)
#'     - q: Tsallis parameter(s) used
#'     - nboot: Number of bootstrap replicates
#'     - method: Bootstrap method ("percentile", "bca", or NA)
#'     - computation_time_sec: Wall-clock time per gene (seconds)
#'     - error: Error message if computation failed, NA_character_ otherwise
#'   
#'   **colData**:
#'     - Inherited from input `se` (sample grouping, pairing info, etc.)
#'   
#'   **metadata** (list):
#'     - summary_stats: list with counts (total, successful, failed)
#'     - elapsed_time_sec: Total computation time
#'     - avg_time_per_gene: Average time per gene
#'     - genes_per_minute: Processing rate
#'     - bootstrap_config: list with bootstrap parameters (nboot, ci, method)
#'     - computation_mode: "sequential" or "parallel"
#'
#' @details
#' **Paired Sample Auto-Detection (NEW FEATURE):**
#' When bootstrap=TRUE, the function automatically detects paired sample metadata from colData:
#' - Searches for columns: "paired_samples", "pair_id", "pair_samples", "subject_id", "patient_id"
#' - If found, uses **pair-respecting bootstrap resampling**:
#'   * Resamples pair indices (not individual samples) with replacement
#'   * Preserves within-pair correlations critical for matched designs
#'   * Maintains statistical validity in paired experimental designs
#' - If paired=TRUE but no pairing detected, falls back to independent bootstrap with warning
#' 
#' **Scientific Justification (Papers validating auto-detection):**
#' - Papers C016: Bootstrap for confidence intervals requires preserving data structure
#' - Papers S102-S109: Paired design standards and statistical methods
#' - Papers I002-I004: Tsallis divergence mathematical foundation
#' 
#' **Computation Mode Selection:**
#' The function automatically selects between sequential and parallel processing:
#' - nthreads=1 (default): Direct sequential loop, minimal overhead
#' - nthreads > 1 & num_genes >= 5: Parallel PSOCK cluster
#' - nthreads > 1 & num_genes < 5: Falls back to sequential (overhead not warranted)
#'
#' **Performance Characteristics:**
#' - With bootstrap=TRUE (default):
#'   - Sequential: ~2-3 seconds per gene (nboot=1000, percentile method)
#'   - Parallel overhead: ~1-2 seconds initial cluster setup
#'   - Break-even point: ~10-20 genes
#' - With bootstrap=FALSE (point estimates only):
#'   - Sequential: ~0.02-0.05 seconds per gene (50-100* faster)
#'
#' **Gene Filtering:**
#' - Always process all genes in se
#'
#' **Database Verification (tsenat_papers.db):**
#' - Tsallis divergence mathematical foundation: Papers I001-I004 validate
#'   divergence formula and q-parameter effects
#' - Bootstrap methodology: Papers C016, S018, S030 validate percentile and BCa
#'   bootstrap for entropy/divergence estimates with confidence level >= 0.95
#' - Transcript aggregation: Paper C105 validates gene-level aggregation
#' - Divergence normalization: Papers C112, S196, S201 validate normalization
#'   approaches for effect size comparability (S197 - DESeq2 independent filtering)
#' @keywords internal
#' @noRd

calculate_divergence <- function(
    se,
    group_col = NULL,
    control_group = NULL,
    q = 1,
    paired = FALSE,
    bootstrap = FALSE,
    nboot = "auto",
    ci = 0.95,
    method = "percentile",
    norm = TRUE,
    log_base = exp(1),
    pseudocount = 0.5,
    nthreads = 1,
    progress = FALSE,
    verbose = TRUE,
    seed = NULL) {

  # =========================================================================
  # INPUT VALIDATION
  # =========================================================================

  # Normalize norm parameter: coerce logical to character for backward compatibility
  if (is.logical(norm)) {
    norm <- if (norm) "range" else "none"
  }
  
  # Validate norm parameter
  norm <- match.arg(norm, choices = c("none", "range", "zscore", 
                                      "log_odds_ratio", "relative_reference"))

  # BUGFIX #4: Auto-sort q parameter for consistent output and q-spectrum analysis
  # Sorts q values in ascending order (rare -> abundant: 0.5 -> 2)
  q <- sort(as.numeric(q))
  if (any(q <= 0)) {
    stop("q parameter must be positive. ",
         "Note: q should be in range (0, 3] for typical use. ",
         "Got: ", paste(q, collapse = ", "))
  }

  if (!methods::is(se, "SummarizedExperiment")) {
    stop("se must be a SummarizedExperiment object")
  }

  # =========================================================================
  # AUTO-DETECT GROUP COLUMN AND CONTROL GROUP
  # =========================================================================
  
  if (is.null(group_col) || is.null(control_group)) {
    auto_groups <- .auto_detect_groups(se)
    
    if (is.null(group_col)) {
      if (is.na(auto_groups$group_col)) {
        stop("Could not auto-detect group column in colData. ",
             "Available columns: ", 
             paste(colnames(SummarizedExperiment::colData(se)), collapse = ", "),
             ". Please specify 'group_col' explicitly.",
             call. = FALSE)
      }
      group_col <- auto_groups$group_col
      if (progress) {
        message("[calculate_divergence] Auto-detected group_col='", group_col, "'")
      }
    }
    
    if (is.null(control_group)) {
      if (is.na(auto_groups$control_group)) {
        stop("Could not auto-detect control_group. Found groups: ",
             paste(auto_groups$groups, collapse = ", "),
             ". Please specify 'control_group' explicitly.",
             call. = FALSE)
      }
      control_group <- auto_groups$control_group
      if (progress) {
        message("[calculate_divergence] Auto-detected control_group='", control_group, "'")
      }
    }
  }

  # Get list of UNIQUE genes to process (not transcripts)
  rd <- SummarizedExperiment::rowData(se)
  
  # Identify which column contains gene symbols/names (prioritize human-readable gene_names over gene_ids)
  # Preference order: gene_name (human-readable) > gene_names (legacy) > gene_id (ensembl/ids) > genes (legacy ids)
  gene_col_candidates <- c("gene_name","gene_id")
  gene_col <- NA_character_
  if (!is.null(rd)) {
    for (col in gene_col_candidates) {
      if (col %in% colnames(rd)) {
        gene_col <- col
        break
      }
    }
  }
  
  # Extract unique gene names/symbols
  if (!is.na(gene_col) && !is.null(rd)) {
    all_gene_names <- unique(as.character(rd[[gene_col]]))
  } else {
    all_gene_names <- rownames(se)
  }
  
  if (is.null(all_gene_names) || length(all_gene_names) == 0) {
    stop("se must have gene identifiers in rowData or rownames")
  }

  # Process unique genes only
  genes_to_process <- all_gene_names
  gene_indices <- seq_along(all_gene_names)

  num_genes <- length(gene_indices)
  
  if (num_genes == 0) {
    stop("No genes to process. ",
         "se gene names (first 3): ", paste(head(all_gene_names, 3), collapse=", "))
  }

  # =========================================================================
  # BOOTSTRAP CONFIGURATION
  # =========================================================================

  if (!is.logical(bootstrap) || length(bootstrap) != 1) {
    stop("bootstrap must be a logical (TRUE/FALSE)")
  }

  if (!bootstrap) {
    nboot <- 0
  }

  # AUTO-SELECT NBOOT WHEN "auto"
  if (bootstrap && identical(nboot, "auto")) {
    num_genes <- nrow(se)
    use_bca <- method == "bca"
    nboot <- suggest_nboot(num_genes, use_bca = use_bca, nthreads = nthreads)
    if (progress) {
      cat("Auto-selected nboot =", nboot, "for", num_genes, "genes\n")
    }
  }

  # THREAD CONFIGURATION
  # =========================================================================

  if (is.null(nthreads)) {
    nthreads <- parallel::detectCores() - 1
    nthreads <- max(1, nthreads)
  }

  if (!is.numeric(nthreads) || nthreads < 1) {
    stop("nthreads must be a positive integer")
  }
  nthreads <- as.integer(nthreads)

  use_parallel <- (num_genes >= 5 && nthreads > 1)

  # =========================================================================
  # PAIRED SAMPLE DETECTION (NEW FEATURE)
  # =========================================================================
  
  pair_ids <- NULL
  pairing_info <- ""
  
  if (bootstrap) {
    # Auto-detect paired samples when bootstrap=TRUE (or when paired=TRUE)
    pair_detected <- .detect_pair_ids(se)
    
    if (pair_detected$num_pairs > 0) {
      pair_ids <- pair_detected$pair_ids
      pairing_info <- sprintf(" [paired: %d unique pairs from '%s' column]", 
                              pair_detected$num_pairs, 
                              pair_detected$column_name)
      
      if (paired == FALSE && progress) {
        cat("NOTE: Paired sample structure detected in '", 
            pair_detected$column_name, "' column.\n",
            "      Using pair-respecting bootstrap resampling.\n\n", sep = "")
      }
    } else {
      if (paired == TRUE && progress) {
        cat("WARNING: paired=TRUE but no pair ID column detected in colData.\n",
            "         Using independent bootstrap resampling instead.\n\n")
      }
    }
  }

  if (progress) {
    mode_desc <- if (bootstrap) {
      paste0("bootstrap with ", nboot, " replicates (", method, ")", pairing_info)
    } else {
      "point estimates only"
    }

    mode_str <- if (use_parallel) "Parallel" else "Sequential"
    cat(mode_str, "mode: ", num_genes, " genes", 
        if (use_parallel) paste0(" on ", nthreads, " threads") else "",
        " [", mode_desc, "]\n\n", sep = "")
  }

  # =========================================================================
  # COMPUTATION (SEQUENTIAL OR PARALLEL)
  # =========================================================================

  start_time <- Sys.time()
  results_list <- list()

  if (!use_parallel) {
    # ====== SEQUENTIAL PROCESSING ======
    for (i in seq_along(gene_indices)) {
      gene_idx <- gene_indices[i]
      target_gene <- all_gene_names[gene_idx]
      gene_name <- target_gene
      gene_start <- Sys.time()

      tryCatch({
        # Get gene-level counts via transcript aggregation
        # (Following Paper I033: gene-level analysis from transcript-level data)
        
        # Find ALL transcripts for this gene
        if (!is.na(gene_col) && !is.null(rd)) {
          gene_transcript_indices <- which(as.character(rd[[gene_col]]) == target_gene)
        } else {
          # Fallback: match by rownames
          gene_transcript_indices <- which(rownames(se) == target_gene)
        }
        
        if (length(gene_transcript_indices) == 0) {
          results_list[[i]] <- list(
            gene_name = gene_name,
            results_per_q = rep(list(list(estimate = NA_real_, lower_ci = NA_real_, upper_ci = NA_real_, method = NA_character_)), length(q)),
            computation_time_sec = as.numeric(Sys.time() - gene_start, units = "secs"),
            error = "No transcripts found for gene"
          )
          next
        }
        
        # Aggregate counts across all transcripts for this gene
        counts_matrix <- as.matrix(SummarizedExperiment::assay(se, "counts")[gene_transcript_indices, , drop = FALSE])
        counts_gene <- colSums(counts_matrix)
        
        groups <- se[[group_col]]
        
        x <- counts_gene[groups == control_group]
        y <- counts_gene[groups != control_group]

        if (length(x) == 0 || length(y) == 0) {
          results_list[[i]] <- list(
            gene_name = gene_name,
            estimate = NA_real_,
            lower_ci = NA_real_,
            upper_ci = NA_real_,
            ci_width = NA_real_,
            q = paste(q, collapse = ","),
            nboot = nboot,
            method = if (bootstrap) method else NA_character_,
            computation_time_sec = as.numeric(Sys.time() - gene_start, units = "secs"),
            error = "Insufficient group samples"
          )
          if (progress && (i %% 10 == 0)) {
            cat("[", i, "/", num_genes, "] Processed\n")
          }
          next
        }

        # Compute divergence for each q value separately and store all results
        gene_results <- list()
        for (j in seq_along(q)) {
          q_val <- q[j]
          
          # Build arguments conditionally to avoid passing NULL pair_ids if not detected
          bootstrap_args <- list(
            x = x, y = y,
            q = q_val, nboot = nboot, ci = ci, method = method,
            log_base = log_base, pseudocount = pseudocount,
            gene_name =  gene_name, print_results = FALSE, seed = seed,
            paired = !is.null(pair_ids)
          )
          
          # Only add pair_ids if it exists (non-NULL)
          if (!is.null(pair_ids)) {
            bootstrap_args$pair_ids <- pair_ids
          }
          
          result <- do.call(calculate_divergence_bootstrap, bootstrap_args)
          gene_results[[j]] <- result
        }

        gene_elapsed <- as.numeric(Sys.time() - gene_start, units = "secs")

        results_list[[i]] <- list(
          gene_name = gene_name,
          results_per_q = gene_results,
          computation_time_sec = gene_elapsed,
          error = NA_character_
        )

        if (progress && (i %% 10 == 0)) {
          elapsed <- as.numeric(Sys.time() - start_time, units = "secs")
          rate <- (i / elapsed) * 60
          cat("[", i, "/", num_genes, "] (", sprintf("%.1f genes/min", rate), ")\n")
        }

      }, error = function(e) {
        results_list[[i]] <<- list(
          gene_name = gene_name,
          results_per_q = rep(list(list(estimate = NA_real_, lower_ci = NA_real_, upper_ci = NA_real_, method = NA_character_)), length(q)),
          computation_time_sec = as.numeric(Sys.time() - gene_start, units = "secs"),
          error = as.character(e$message)
        )
      })
    }

  } else {
    # ====== PARALLEL PROCESSING ======

    cl <- parallel::makeCluster(nthreads, type = "PSOCK")
    on.exit(parallel::stopCluster(cl), add = TRUE)

    # Export all necessary variables to cluster nodes
    parallel::clusterExport(cl, 
      c("calculate_divergence_bootstrap", ".tsallis_divergence_scalar",
        "rd", "gene_col", "all_gene_names", "se", "q", "nboot", "ci", "method",
        "log_base", "pseudocount", "seed", "group_col", "control_group", "pair_ids"),
      envir = environment())

    parallel::clusterCall(cl, function() {
      requireNamespace("SummarizedExperiment", quietly = TRUE)
    })

    results_list <- parallel::parLapply(cl, seq_along(gene_indices), function(i) {
      gene_idx <- gene_indices[i]
      target_gene <- all_gene_names[gene_idx]
      gene_name <- target_gene
      gene_start <- Sys.time()

      tryCatch({
        # Get gene-level counts via transcript aggregation
        # (Following Paper I033: gene-level analysis from transcript-level data)
        
        # Find ALL transcripts for this gene
        if (!is.na(gene_col) && !is.null(rd)) {
          gene_transcript_indices <- which(as.character(rd[[gene_col]]) == target_gene)
        } else {
          # Fallback: match by rownames
          gene_transcript_indices <- which(rownames(se) == target_gene)
        }
        
        if (length(gene_transcript_indices) == 0) {
          return(list(
            gene_name = gene_name,
            results_per_q = rep(list(list(estimate = NA_real_, lower_ci = NA_real_, upper_ci = NA_real_, method = NA_character_)), length(q)),
            computation_time_sec = as.numeric(Sys.time() - gene_start, units = "secs"),
            error = "No transcripts found for gene"
          ))
        }
        
        # Aggregate counts across all transcripts for this gene
        counts_matrix <- as.matrix(SummarizedExperiment::assay(se, "counts")[gene_transcript_indices, , drop = FALSE])
        counts_gene <- colSums(counts_matrix)
        
        groups <- se[[group_col]]

        x <- counts_gene[groups == control_group]
        y <- counts_gene[groups != control_group]

        if (length(x) == 0 || length(y) == 0) {
          return(list(
            gene_name = gene_name,
            results_per_q = rep(list(list(estimate = NA_real_, lower_ci = NA_real_, upper_ci = NA_real_, method = NA_character_)), length(q)),
            computation_time_sec = as.numeric(Sys.time() - gene_start, units = "secs"),
            error = "Insufficient group samples"
          ))
        }

        # Compute divergence for each q value separately and store all results
        gene_results <- list()
        for (j in seq_along(q)) {
          q_val <- q[j]
          
          #Build arguments conditionally to avoid passing NULL pair_ids if not detected
          bootstrap_args <- list(
            x = x, y = y,
            q = q_val, nboot = nboot, ci = ci, method = method,
            log_base = log_base, pseudocount = pseudocount,
            gene_name = gene_name, print_results = FALSE, seed = seed,
            paired = !is.null(pair_ids)
          )
          
          # Only add pair_ids if it exists (non-NULL)
          if (!is.null(pair_ids)) {
            bootstrap_args$pair_ids <- pair_ids
          }
          
          result <- do.call(calculate_divergence_bootstrap, bootstrap_args)
          gene_results[[j]] <- result
        }

        gene_elapsed <- as.numeric(Sys.time() - gene_start, units = "secs")

        list(
          gene_name = gene_name,
          results_per_q = gene_results,
          computation_time_sec = gene_elapsed,
          error = NA_character_
        )
      }, error = function(e) {
        list(
          gene_name = gene_name,
          results_per_q = rep(list(list(estimate = NA_real_, lower_ci = NA_real_, upper_ci = NA_real_, method = NA_character_)), length(q)),
          computation_time_sec = as.numeric(Sys.time() - gene_start, units = "secs"),
          error = as.character(e$message)
        )
      })
    })
  }

  elapsed <- as.numeric(Sys.time() - start_time, units = "secs")

  # =========================================================================
  # RESULTS COMPILATION - Build genes * q matrix
  # =========================================================================

  # Extract results into proper per-q structure
  num_q_vals <- length(q)
  max_genes <- num_genes
  
  # Initialize matrices: genes * q
  assay_matrix <- matrix(NA_real_, nrow = max_genes, ncol = num_q_vals,
                         dimnames = list(NULL, paste0("q_", q)))
  
  row_data_df <- data.frame(
    gene_name = character(max_genes),
    error = character(max_genes),
    computation_time_sec = numeric(max_genes),
    stringsAsFactors = FALSE
  )
  
  # Add columns for each q value's metadata
  for (j in seq_len(num_q_vals)) {
    row_data_df[[paste0("estimate_q", q[j])]] <- NA_real_
    row_data_df[[paste0("lower_ci_q", q[j])]] <- NA_real_
    row_data_df[[paste0("upper_ci_q", q[j])]] <- NA_real_
    row_data_df[[paste0("ci_width_q", q[j])]] <- NA_real_
    row_data_df[[paste0("method_q", q[j])]] <- NA_character_
    row_data_df[[paste0("nboot_q", q[j])]] <- NA_integer_
  }
  
  for (i in seq_len(max_genes)) {
    res <- results_list[[i]]
    row_data_df$gene_name[i] <- res$gene_name
    row_data_df$error[i] <- res$error
    row_data_df$computation_time_sec[i] <- res$computation_time_sec
    
    # Extract results for each q value
    for (j in seq_len(num_q_vals)) {
      if (j <= length(res$results_per_q)) {
        result_j <- res$results_per_q[[j]]
        # Check if result is valid (estimate is not NA)
        # BUGFIX: Ensure we handle legitimate zero divergences (they can occur)
        # Only skip if estimate is truly NA or NULL
        if (!is.null(result_j) && !is.na(result_j$estimate)) {
          assay_matrix[i, j] <- result_j$estimate
          q_col_name <- paste0("estimate_q", q[j])
          row_data_df[[q_col_name]][i] <- result_j$estimate
          row_data_df[[paste0("lower_ci_q", q[j])]][i] <- result_j$lower_ci
          row_data_df[[paste0("upper_ci_q", q[j])]][i] <- result_j$upper_ci
          row_data_df[[paste0("ci_width_q", q[j])]][i] <- result_j$upper_ci - result_j$lower_ci
          row_data_df[[paste0("method_q", q[j])]][i] <- result_j$method
          row_data_df[[paste0("nboot_q", q[j])]][i] <- nboot
        } else if (progress && i <= 3) {
          # Log first few genes with NA estimates for debugging rare q-values
        }
      }
    }
  }
  
  rownames(row_data_df) <- row_data_df$gene_name
  rownames(assay_matrix) <- row_data_df$gene_name
  
  # Classify per-q patterns: RARE_DRIVEN, ABUNDANT_DRIVEN, or BALANCED
  row_data_df$per_q_pattern <- NA_character_
  if (length(q) > 1) {
    for (i in seq_len(nrow(row_data_df))) {
      # Get divergence values across q for this gene
      per_q_divs <- assay_matrix[i, ]
      names(per_q_divs) <- paste0("q_", q)
      
      # Only classify if we have valid values
      if (sum(!is.na(per_q_divs)) >= 2) {
        pattern <- classify_q_pattern(per_q_divs)
        row_data_df$per_q_pattern[i] <- if (is.na(pattern)) "UNCLASSIFIED" else pattern
      }
    }
  }
  
  # Populate generic estimate/lower_ci/upper_ci columns using reference q value (q=1)
  # This is needed for downstream functions like effect_sizes_divergence()
  q_ref <- 1.0
  q_idx <- which.min(abs(q - q_ref))  # Find closest q to 1.0
  if (length(q_idx) > 0 && q_idx <= length(q)) {
    ref_q <- q[q_idx]
    estimate_col <- paste0("estimate_q", ref_q)
    lower_ci_col <- paste0("lower_ci_q", ref_q)
    upper_ci_col <- paste0("upper_ci_q", ref_q)
    ci_width_col <- paste0("ci_width_q", ref_q)
    
    # Copy per-q columns to generic columns for downstream compatibility
    if (estimate_col %in% colnames(row_data_df)) {
      row_data_df$estimate <- row_data_df[[estimate_col]]
      row_data_df$lower_ci <- row_data_df[[lower_ci_col]]
      row_data_df$upper_ci <- row_data_df[[upper_ci_col]]
      row_data_df$ci_width <- row_data_df[[ci_width_col]]
    }
  }
  
  # Apply normalization if requested (BEFORE creating SE)
  if (norm != "none") {
    if (progress) {
      cat(sprintf("Applying '%s' normalization to divergence estimates...\n", norm))
    }
    
    if (norm == "range") {
      # Range normalization [0,1] per q value (normalize within each column)
      for (j in seq_len(ncol(assay_matrix))) {
        col_vals <- assay_matrix[, j]
        valid_vals <- col_vals[!is.na(col_vals)]
        if (length(valid_vals) > 1) {
          min_val <- min(valid_vals, na.rm = TRUE)
          max_val <- max(valid_vals, na.rm = TRUE)
          if (max_val > min_val) {
            assay_matrix[, j] <- (col_vals - min_val) / (max_val - min_val)
            
            # BUGFIX #5 & Hidden Bug: Apply same normalization to estimate AND CI bounds (not just CIs)
            estimate_col <- paste0("estimate_q", q[j])
            lower_col <- paste0("lower_ci_q", q[j])
            upper_col <- paste0("upper_ci_q", q[j])
            
            if (estimate_col %in% colnames(row_data_df)) {
              row_data_df[[estimate_col]] <- (row_data_df[[estimate_col]] - min_val) / (max_val - min_val)
              row_data_df[[lower_col]] <- (row_data_df[[lower_col]] - min_val) / (max_val - min_val)
              row_data_df[[upper_col]] <- (row_data_df[[upper_col]] - min_val) / (max_val - min_val)
            }
          }
        }
      }
    } else if (norm == "zscore") {
      # Z-score normalization per q value (normalize within each column)
      for (j in seq_len(ncol(assay_matrix))) {
        col_vals <- assay_matrix[, j]
        valid_vals <- col_vals[!is.na(col_vals)]
        if (length(valid_vals) > 1) {
          mean_val <- mean(valid_vals, na.rm = TRUE)
          sd_val <- sd(valid_vals, na.rm = TRUE)
          if (sd_val > 0) {
            assay_matrix[, j] <- (col_vals - mean_val) / sd_val
            
            # BUGFIX #5 & Hidden Bug: Apply same normalization to estimate AND CI bounds (not just CIs)
            estimate_col <- paste0("estimate_q", q[j])
            lower_col <- paste0("lower_ci_q", q[j])
            upper_col <- paste0("upper_ci_q", q[j])
            
            if (estimate_col %in% colnames(row_data_df)) {
              row_data_df[[estimate_col]] <- (row_data_df[[estimate_col]] - mean_val) / sd_val
              row_data_df[[lower_col]] <- (row_data_df[[lower_col]] - mean_val) / sd_val
              row_data_df[[upper_col]] <- (row_data_df[[upper_col]] - mean_val) / sd_val
            }
          }
        }
      }
    } else if (norm == "log_odds_ratio") {
      # Log ratio to theoretical maximum
      max_div <- 1.0
      valid_mask <- !is.na(assay_matrix) & assay_matrix > 0
      if (any(valid_mask)) {
        assay_matrix[valid_mask] <- log(assay_matrix[valid_mask] / max_div)
        
        # BUGFIX #5 & Hidden Bug: Apply same normalization to estimate AND CI bounds (not just CIs)
        for (j in seq_len(ncol(assay_matrix))) {
          estimate_col <- paste0("estimate_q", q[j])
          lower_col <- paste0("lower_ci_q", q[j])
          upper_col <- paste0("upper_ci_q", q[j])
          
          if (estimate_col %in% colnames(row_data_df)) {
            mask <- !is.na(row_data_df[[estimate_col]]) & row_data_df[[estimate_col]] > 0
            row_data_df[[estimate_col]][mask] <- log(row_data_df[[estimate_col]][mask] / max_div)
            
            mask <- !is.na(row_data_df[[lower_col]]) & row_data_df[[lower_col]] > 0
            row_data_df[[lower_col]][mask] <- log(row_data_df[[lower_col]][mask] / max_div)
            
            mask <- !is.na(row_data_df[[upper_col]]) & row_data_df[[upper_col]] > 0
            row_data_df[[upper_col]][mask] <- log(row_data_df[[upper_col]][mask] / max_div)
          }
        }
      }
    } else if (norm == "relative_reference") {
      # Ratio to reference group mean per q value
      for (j in seq_len(ncol(assay_matrix))) {
        col_vals <- assay_matrix[, j]
        valid_vals <- col_vals[!is.na(col_vals)]
        if (length(valid_vals) > 0) {
          reference_mean <- mean(valid_vals, na.rm = TRUE)
          if (reference_mean > 0) {
            assay_matrix[, j] <- col_vals / reference_mean
            
            # BUGFIX #5 & Hidden Bug: Apply same normalization to estimate AND CI bounds (not just CIs)
            estimate_col <- paste0("estimate_q", q[j])
            lower_col <- paste0("lower_ci_q", q[j])
            upper_col <- paste0("upper_ci_q", q[j])
            
            if (estimate_col %in% colnames(row_data_df)) {
              row_data_df[[estimate_col]] <- row_data_df[[estimate_col]] / reference_mean
              row_data_df[[lower_col]] <- row_data_df[[lower_col]] / reference_mean
              row_data_df[[upper_col]] <- row_data_df[[upper_col]] / reference_mean
            }
          }
        }
      }
    }
  }

  # Summary statistics
  num_success <- sum(is.na(row_data_df$error))
  num_errors <- num_genes - num_success

  if (progress) {
    cat("\n", paste(rep("=", 70), collapse = ""), "\n", sep = "")
    cat("DIVERGENCE COMPUTATION COMPLETE\n")
    cat(paste(rep("=", 70), collapse = ""), "\n\n", sep = "")

    cat("Summary:\n")
    cat("  Genes processed:        ", num_genes, "\n")
    cat("  Successful:             ", num_success, "\n")
    cat("  Failed:                 ", num_errors, "\n")
    cat("  Total elapsed time:     ", sprintf("%.1f seconds", elapsed), "\n")
    cat("  Average per gene:       ", sprintf("%.2f seconds", elapsed / num_genes), "\n")
    cat("  Genes per minute:       ", sprintf("%.1f", (num_genes / elapsed) * 60), "\n\n")

    if (num_errors > 0) {
      cat("Failed genes:\n")
      failed <- row_data_df[!is.na(row_data_df$error), ]
      for (i in seq_len(min(10, nrow(failed)))) {
        cat(sprintf("  [%d] %s: %s\n", i, failed$gene_name[i], failed$error[i]))
      }
      if (num_errors > 10) {
        cat("  ... and", num_errors - 10, "more\n")
      }
      cat("\n")
    }
  }

  # =========================================================================
  # CREATE SUMMARIZED EXPERIMENT
  # =========================================================================

  # Create SummarizedExperiment with:
  # - assay: genes * 1 matrix (one divergence value per gene)
  # - rowData: gene metadata (estimates, CIs, computation details)
  # - colData: "divergence" as the single pseudo-sample (represents the divergence assay)
  # - metadata: processing parameters

  # Initialize assays list with divergence
  assays_list <- list(divergence = assay_matrix)

  # Create colData for each q value column in the assay
  col_data_output <- data.frame(
    q_value = q,
    sample_type = rep("divergence_estimate", length(q)),
    row.names = paste0("q_", q)
  )

  result_se <- SummarizedExperiment::SummarizedExperiment(
    assays = assays_list,
    rowData = row_data_df,
    colData = col_data_output,
    metadata = list(
      summary_stats = list(
        total_genes = num_genes,
        successful = num_success,
        failed = num_errors
      ),
      elapsed_time_sec = elapsed,
      avg_time_per_gene = elapsed / num_genes,
      genes_per_minute = (num_genes / elapsed) * 60,
      bootstrap_config = list(
        nboot = nboot,
        ci = ci,
        method = method
      ),
      normalization = norm,
      computation_mode = if (use_parallel) "parallel" else "sequential"
    )
  )

  return(result_se)
}


# -------------------------------------------------------------------------
# HELPER FUNCTION: Bootstrap CI Computation
# -------------------------------------------------------------------------

#' Calculate Tsallis Divergence Bootstrap Confidence Intervals
#'
#' Internal helper function for compute_divergence_bootstrap.
#' Computes bootstrap confidence intervals for Tsallis divergence.
#'
#' @keywords internal
#' @noRd
calculate_divergence_bootstrap <- function(
    x, y,
    q = 1,
    nboot = 1000,
    ci = 0.95,
    method = "percentile",
    log_base = exp(1),
    pseudocount = 0.5,
    gene_name = NA_character_,
    print_results = FALSE,
    seed = NULL,
    paired = FALSE,
    pair_ids = NULL) {

  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Compute point estimate
  point_est <- .tsallis_divergence_scalar(x, y, q, pseudocount, log_base)

  # Bootstrap confidence interval
  if (nboot > 0) {
    bootstrap_dist <- numeric(nboot)
    
    # Determine if using paired resampling
    use_paired_bootstrap <- paired && !is.null(pair_ids)

    for (b in seq_len(nboot)) {
      if (use_paired_bootstrap) {
        # Paired bootstrap: resample pair indices to preserve within-pair correlation
        # while maintaining equal group sizes
        x_names <- names(x)
        y_names <- names(y)
        all_names <- c(x_names, y_names)
        pair_ids_subset <- pair_ids[all_names]
        
        # Get unique pairs and resample indices with replacement
        unique_pairs <- unique(pair_ids_subset)
        num_pairs <- length(unique_pairs)
        resampled_pair_indices <- sample(seq_len(num_pairs), size = num_pairs, replace = TRUE)
        
        # Collect samples from resampled pairs, maintaining group structure
        x_boot <- numeric(0)
        y_boot <- numeric(0)
        
        for (idx in resampled_pair_indices) {
          pair_id <- unique_pairs[idx]
          pair_mask <- pair_ids_subset == pair_id
          pair_samples <- names(pair_ids_subset)[pair_mask]
          
          for (sample in pair_samples) {
            if (sample %in% x_names) {
              x_boot <- c(x_boot, x[sample])
            } else if (sample %in% y_names) {
              y_boot <- c(y_boot, y[sample])
            }
          }
        }
      } else {
        # Independent bootstrap: standard resampling (each group independently)
        x_boot <- sample(x, size = length(x), replace = TRUE)
        y_boot <- sample(y, size = length(y), replace = TRUE)
      }

      bootstrap_dist[b] <- .tsallis_divergence_scalar(x_boot, y_boot, q, pseudocount, log_base)
    }

    alpha <- (1 - ci) / 2

    if (method == "percentile") {
      lower_ci <- stats::quantile(bootstrap_dist, probs = alpha, names = FALSE)
      upper_ci <- stats::quantile(bootstrap_dist, probs = 1 - alpha, names = FALSE)
    } else if (method == "bca") {
      # BCA not appropriate for divergence (requires two-sample jackknife)
      # Fall back to percentile method which is valid for any divergence
      lower_ci <- stats::quantile(bootstrap_dist, probs = alpha, names = FALSE)
      upper_ci <- stats::quantile(bootstrap_dist, probs = 1 - alpha, names = FALSE)
    }
  } else {
    lower_ci <- NA_real_
    upper_ci <- NA_real_
  }

  return(list(
    estimate = point_est,
    lower_ci = lower_ci,
    upper_ci = upper_ci,
    q = q,
    nboot = nboot,
    method = if (nboot > 0) method else NA_character_
  ))
}


#' Compute Tsallis Divergence Between Two Count Vectors
#'
#' Computes scalar Tsallis divergence D_q(p || q) using the Furuichi formula.
#'
#' @keywords internal
#' @noRd
.tsallis_divergence_scalar <- function(x, y, q_val, pseudocount = 0.5, log_base = exp(1)) {
  # Validate input vectors
  if (length(x) == 0 || length(y) == 0) {
    return(NA_real_)
  }
  
  # Normalize to probabilities
  p <- (x + pseudocount) / (sum(x) + length(x) * pseudocount)
  r <- (y + pseudocount) / (sum(y) + length(y) * pseudocount)

  if (any(is.na(p)) || any(is.na(r))) {
    return(NA_real_)
  }
  
  # BUGFIX #2: Add explicit safeguard for near-zero probabilities
  # Prevents NaN/Inf from log(0) or extremely small values in power operations
  min_prob <- 1e-10
  p[p < min_prob] <- min_prob
  r[r < min_prob] <- min_prob
  
  # Re-normalize to maintain probability constraint (sum = 1)
  p <- p / sum(p)
  r <- r / sum(r)

  # Compute Tsallis divergence using correct formula from Paper I004
  # D_q(p||r) with D_q >= 0 and equality iff p = r
  # BUGFIX: Ensure formula is applied correctly for all q values
  
  if (abs(q_val - 1) < 0.01) {
    # KL divergence (special case q -> 1): lim_{q->1} D_q = sum(p*log(p/r))
    div <- sum(p * log(p / r), na.rm = TRUE)
  } else if (q_val > 0 && q_val != 1) {
    # Standard Tsallis divergence formula: D_q(p||r) = (1/(q-1)) * (1 - sum(p^q * r^(1-q)))
    # This ensures D_q >= 0 and is asymmetric in p, r
    # CRITICAL: Ensure p and r vectors are properly aligned
    p_power <- p^q_val
    r_power <- r^(1 - q_val)
    
    # Check for numerical issues (inf, nan, underflow)
    if (any(is.nan(p_power)) || any(is.infinite(p_power)) ||
        any(is.nan(r_power)) || any(is.infinite(r_power))) {
      # Log-space computation for numerical stability when q is far from 1
      log_p_power <- q_val * log(p + 1e-10)
      log_r_power <- (1 - q_val) * log(r + 1e-10)
      sum_term <- sum(exp(log_p_power + log_r_power), na.rm = TRUE)
    } else {
      sum_term <- sum(p_power * r_power, na.rm = TRUE)
    }
    
    div <- (1 - sum_term) / (q_val - 1)
  } else {
    # Invalid q value
    return(NA_real_)
  }

  if (is.nan(div) || !is.finite(div)) {
    return(NA_real_)
  }

  # Apply log_base normalization CONSISTENTLY for all q values
  # This ensures consistent scaling across multi-q spectrum analysis
  if (log_base != exp(1)) {
    div <- div / log(log_base)
  }

  # BUG FIX: Handle sign correctly for q < 1
  # When q < 1, (q_val - 1) is negative, so the formula naturally produces
  # a positive divergence. We must take absolute value and ensure non-negativity.
  # Divergence should always be >= 0.
  return(abs(div))
}


#' Classify a per-q divergence spectrum into biological pattern types
#'
#' When Tsallis divergence has been computed across multiple q values
#' for a gene, the resulting vector can be summarised by its trend across the
#' spectrum. This helper computes the Pearson correlation between the numeric
#' q values and their corresponding divergences and interprets the
#' slope as one of three categories:
#' \describe{
#'   \item{RARE_DRIVEN}{Divergence decreases with q (high at q=0.5);
#'     indicates changes driven by low-abundance isoforms.}
#'   \item{ABUNDANT_DRIVEN}{Divergence increases with q (high at q=2);
#'     indicates shifts among the most abundant transcripts.}
#'   \item{BALANCED}{Little or no trend across q; effects are proportional.}
#' }
#'
#' If the input vector is too short, contains only NAs, or the
#' correlation cannot be calculated, NA is returned.
#'
#' @param per_q_divs Named numeric vector of divergences. Names should be of form
#'   "q_0.5", "q_1.0", etc. (or similar with separators "_", "=", or ".").
#' @param threshold Numeric; correlation threshold for classification. Absolute
#'   correlation values below this threshold are classified as BALANCED (default: 0.5).
#' @return Character scalar giving the pattern type, or NA if classification
#'   cannot be performed.
#' @keywords internal
#' @noRd
#' @examples
#' per_q <- c(q_0.5=0.5, q_1=0.3, q_2=0.1)
#' classify_q_pattern(per_q)
#'
#' # handling missing values
#' classify_q_pattern(c(q_0.5=NA, q_1=0.2))
classify_q_pattern <- function(per_q_divs, threshold = 0.5) {
  # Input validation
  if (!is.numeric(per_q_divs) || length(per_q_divs) < 2) {
    return(NA_character_)
  }
  
  # Validate names format
  nm <- names(per_q_divs)
  if (is.null(nm) || any(is.na(nm))) {
    return(NA_character_)
  }
  
  # Handle flexible name formats: try multiple patterns
  q_vals <- NA
  
  # Try format: "q_0.5" (standard)
  if (all(grepl("^q_", nm))) {
    q_vals <- as.numeric(gsub("^q_", "", nm))
  } else if (all(grepl("^q[_=.]", nm))) {
    # Try other separators
    q_vals <- as.numeric(gsub("^q[_=.]", "", gsub("_", ".", nm)))
  }
  
  # Fallback: try to extract numeric directly after "q"
  if (all(is.na(q_vals))) {
    q_vals <- suppressWarnings(as.numeric(gsub("^q", "", nm)))
  }
  
  # If we still can't extract numeric q values, return NA
  if (any(is.na(q_vals))) {
    return(NA_character_)
  }
  
  # Check if all divergence values are NA
  if (all(is.na(per_q_divs))) {
    return(NA_character_)
  }
  
  # Need at least 2 non-NA pairs for correlation
  valid_pairs <- !is.na(per_q_divs)
  if (sum(valid_pairs) < 2) {
    return(NA_character_)
  }
  
  # IMPROVED CLASSIFICATION: Compare q-regions instead of just endpoints
  # Split into rare-region (q < 1) and abundant-region (q >= 1)
  rare_mask <- q_vals < 1
  abund_mask <- q_vals >= 1
  
  # Calculate median divergence in each region
  if (sum(rare_mask & valid_pairs) > 0) {
    rare_div_median <- median(per_q_divs[rare_mask & valid_pairs], na.rm = TRUE)
  } else {
    rare_div_median <- NA
  }
  
  if (sum(abund_mask & valid_pairs) > 0) {
    abund_div_median <- median(per_q_divs[abund_mask & valid_pairs], na.rm = TRUE)
  } else {
    abund_div_median <- NA
  }
  
  # If we have both regions, compare them
  if (!is.na(rare_div_median) && !is.na(abund_div_median)) {
    ratio <- rare_div_median / abund_div_median
    
    # Use a proper ratio threshold (not the correlation threshold)
    # Ratio threshold should be > 1 to distinguish RARE from ABUNDANT
    # Use 1.3 as the ratio threshold (30% difference = sensitive but not overly permissive)
    ratio_threshold <- 1.3
    
    # RARE_DRIVEN: rare region has notably higher divergence (ratio > threshold)
    # ABUNDANT_DRIVEN: abundant region has notably higher divergence (ratio < 1/threshold)
    # BALANCED: similar divergence across regions (ratio near 1)
    
    if (ratio > ratio_threshold) {
      return("RARE_DRIVEN")
    } else if (ratio < 1 / ratio_threshold) {
      return("ABUNDANT_DRIVEN")
    } else {
      return("BALANCED")
    }
  }
  
  # Fallback: Use original correlation-based approach
  slope <- suppressWarnings(cor(q_vals, per_q_divs, use = "complete.obs"))
  
  # If correlation is NA (e.g., constant divergence), treat as balanced
  if (is.na(slope)) {
    return("BALANCED")
  }
  
  # Classify based on slope and threshold
  abs_slope <- abs(slope)
  
  if (abs_slope < threshold) {
    "BALANCED"
  } else if (slope < -threshold) {
    "RARE_DRIVEN"
  } else if (slope > threshold) {
    "ABUNDANT_DRIVEN"
  } else {
    "BALANCED"
  }
}


#' Merge LMM Interaction Results with Tsallis Divergence Effect Sizes
#'
#' Combines LMM interaction test p-values with pre-computed Tsallis divergence
#' effect sizes and bootstrap confidence intervals across ALL q values. This is a 
#' **data merger**, not a model fitter--all statistical computation happens upstream in:
#' - `calculate_lm_interaction()` -> LMM p-values
#' - `calculate_divergence()` -> Divergence estimates and CIs for multiple q
#'
#' This function merges the results into a single data frame for downstream
#' interpretation. When multiple q values are present, effect sizes are computed
#' for each q to capture the full biological spectrum (rare->abundant isoforms).
#'
#' **Architecture:**
#' ```
#' Input 1: LMM results (from calculate_lm_interaction)
#'   - gene names
#'   - adj_p_interaction values
#'
#' Input 2: Divergence SE (from calculate_divergence)
#'   - gene names
#'   - divergence estimates and bootstrap CIs for each q value
#'
#' Output: Merged data frame
#'   - gene name
#'   - statistical significance (p-value)
#'   - effect magnitude for EACH q: D_q, lower_ci_q, upper_ci_q
#' ```
#'
#' @param lm_res A data frame of LMM interaction test results from `calculate_lm_interaction()`,
#'   with columns: `gene` (character, gene name), `adj_p_interaction` (numeric, multiple-test adjusted p-value).
#'   Genes with adj_p_interaction below `significance_threshold` are included.
#'
#' @param divergence_results_se A SummarizedExperiment from `calculate_divergence()`,
#'   containing rowData with columns: `gene_name` and either:
#'   - Generic: `estimate`, `lower_ci`, `upper_ci` (single q-value results), OR
#'   - Per-q: `estimate_q*`, `lower_ci_q*`, `upper_ci_q*` (multiple q-values)
#'   All divergence-related data is self-contained in this object.
#'
#' @param significance_threshold Numeric; p-value threshold for filtering significant
#'   genes (default: 0.05). Only genes with adj_p_interaction < threshold are included.
#'
#' @param enrich_per_q_pattern Logical; if TRUE (default), adds a 'per_q_pattern' column
#'   to the output data frame containing comma-separated divergence values across the
#'   q spectrum for each gene. This column enables visualization and classification
#'   of whether treatment effects are driven by rare (low-q) or abundant (high-q)
#'   isoforms. Set to FALSE to reduce output size if this annotation is not needed.
#'
#' @param verbose Logical; if TRUE, print detailed validation and merge statistics
#'   to console (default: TRUE). Shows counts of passed, skipped, and failed genes.
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{\code{interaction_results}}{Data frame (genes * columns) with merged results:
#'       - `gene`: Gene name (character)
#'       - `p_value_interaction`: LMM adjusted p-value for q:group interaction
#'       - `slope_diff`: q:group interaction slope coefficient (if present in lm_res)
#'       - For EACH q value found: 
#'         - `effect_size_D_q*`: Absolute Tsallis divergence at q
#'         - `D_q*_lower_ci`: Bootstrap lower confidence bound
#'         - `D_q*_upper_ci`: Bootstrap upper confidence bound
#'     }
#'     \item{\code{validation_stats}}{List with merge quality metrics:
#'       - `total_genes`: Total significant genes from LMM
#'       - `passed_lmm`: Successfully merged with divergence data
#'       - `failed_missing_divergence`: Missing or NA divergence estimate
#'       - `other_errors`: Other processing failures
#'       - `q_values`: Numeric vector of q-values processed
#'     }
#'   }
#'
#' @details
#' **Interpretation of Effect Sizes:**
#'
#' Tsallis divergence D_q quantifies the information-theoretic distance between
#' control and treatment isoform distributions at each q-value:
#' - D_q > 0.05: Small effect size
#' - D_q > 0.10: Medium effect size (meaningful biological significance)
#' - D_q > 0.20: Large effect size
#'
#' **Different q-values capture different biological scales:**
#' - q=0.5: Rare (low-abundance) isoforms dominate
#' - q=1.0: Shannon entropy (balanced across abundances)
#' - q=2.0: Common (high-abundance) isoforms dominate
#'
#' Examining the divergence spectrum across q reveals whether treatment effects
#' are driven by rare transcripts (high D at low q) or abundant transcripts (high D at high q).
#'
#' For paired designs, divergence is computed separately within each pair,
#' then averaged to account for pairing structure.
#'
#' **Database References:**
#' - Papers I002-I004: Tsallis divergence mathematical foundation
#' - Papers C016: Bootstrap CI computation respecting data structure
#' - Papers S197: Quality filtering and effect size thresholds
#'
#' @keywords internal
#' @noRd
effect_sizes_divergence <- function(
    lm_res,
    divergence_results_se,
    significance_threshold = 0.05,
    enrich_per_q_pattern = TRUE,
    verbose = TRUE) {

  # =========================================================================
  # INPUT VALIDATION
  # =========================================================================

  if (!is.data.frame(lm_res)) {
    stop("lm_res must be a data frame")
  }

  if (!("gene" %in% colnames(lm_res)) || !("adj_p_interaction" %in% colnames(lm_res))) {
    stop("lm_res must have columns 'gene' and 'adj_p_interaction'")
  }

  if (!methods::is(divergence_results_se, "SummarizedExperiment")) {
    stop("divergence_results_se must be a SummarizedExperiment from calculate_divergence()")
  }

  # Validate divergence SE has required columns in rowData
  rd <- SummarizedExperiment::rowData(divergence_results_se)
  
  if (!("gene_name" %in% colnames(rd))) {
    stop("divergence_results_se rowData must have 'gene_name' column")
  }

  # =========================================================================
  # ALIGN DATASETS
  # =========================================================================
  # Filter lm_res to include only genes present in divergence_results_se
  # This prevents row mismatches and ensures robust merging
  
  # Extract gene identifiers from divergence_results_se
  if ("gene_name" %in% colnames(rd)) {
    divergence_genes <- as.character(rd$gene_name)
  } else {
    divergence_genes <- as.character(rownames(divergence_results_se))
  }
  
  # Get lm_res gene identifiers (try rownames first, then gene column)
  lm_res_genes <- rownames(lm_res)
  if (length(lm_res_genes) == 0 || all(is.na(lm_res_genes))) {
    lm_res_genes <- as.character(lm_res$gene)
  }
  
  # Find matching genes and filter lm_res
  matching_idx <- lm_res_genes %in% divergence_genes
  n_before_filter <- nrow(lm_res)
  n_after_filter <- sum(matching_idx)
  
  if (n_after_filter > 0) {
    lm_res <- lm_res[matching_idx, , drop = FALSE]
    if (verbose) {
      cat("[effect_sizes_divergence] Gene alignment:\n")
      cat("  - lm_res before filtering: ", n_before_filter, " genes\n", sep = "")
      cat("  - lm_res after filtering: ", n_after_filter, " genes\n", sep = "")
      cat("  - Genes filtered out: ", n_before_filter - n_after_filter, "\n\n", sep = "")
    }
  } else {
    if (verbose) {
      cat("[effect_sizes_divergence] WARNING: No matching genes found between lm_res and divergence_results_se\n\n")
    }
  }
  
  # Auto-detect available q values from per-q columns
  estimate_cols <- grep("^estimate_q", colnames(rd), value = TRUE)
  
  if (length(estimate_cols) == 0) {
    # Check for generic columns (single q result)
    if (!all(c("estimate", "lower_ci", "upper_ci") %in% colnames(rd))) {
      stop("divergence_results_se rowData must have either:\n",
           "  - Generic columns: 'estimate', 'lower_ci', 'upper_ci', OR\n",
           "  - Per-q columns: 'estimate_q*', 'lower_ci_q*', 'upper_ci_q*'")
    }
    q_values <- NA_real_
    use_generic <- TRUE
    if (verbose) {
      cat("[effect_sizes_divergence] Using generic divergence columns (single q-value results)\n")
    }
  } else {
    # Extract q values from column names
    q_values <- as.numeric(sub("estimate_q", "", estimate_cols))
    q_values <- sort(q_values)  # Sort for consistent output
    use_generic <- FALSE
    
    if (verbose) {
      cat("[effect_sizes_divergence] Detected per-q columns for q values: ",
          paste(q_values, collapse = ", "), "\n")
    }
  }

  # =========================================================================
  # FILTERING
  # =========================================================================

  significant_genes <- lm_res$gene[lm_res$adj_p_interaction < significance_threshold]

  if (verbose) {
    cat("\n**Filtering effect size analysis to significant genes:**\n")
    cat("- Genes with adj_p_interaction <", significance_threshold, ":", 
        length(significant_genes), "\n\n")
  }

  # Initialize empty results data frame with columns for each q value
  interaction_results <- data.frame(
    gene = character(0),
    p_value_interaction = numeric(0),
    slope_diff = numeric(0),
    stringsAsFactors = FALSE
  )
  
  # Add per-q effect size columns
  if (use_generic) {
    interaction_results$effect_size_D <- numeric(0)
    interaction_results$D_lower_ci <- numeric(0)
    interaction_results$D_upper_ci <- numeric(0)
  } else {
    for (q_val in q_values) {
      q_label <- gsub("\\.", "_", as.character(q_val))  # Replace . with _ for column names
      interaction_results[[paste0("effect_size_D_q", q_label)]] <- numeric(0)
      interaction_results[[paste0("D_q", q_label, "_lower_ci")]] <- numeric(0)
      interaction_results[[paste0("D_q", q_label, "_upper_ci")]] <- numeric(0)
    }
  }

  if (length(significant_genes) == 0) {
    if (verbose) {
      cat("No genes with significant q*group interaction detected.\n")
    }
    return(list(
      interaction_results = interaction_results,
      validation_stats = list(
        total_genes = 0,
        passed_lmm = 0,
        failed_missing_divergence = 0,
        other_errors = 0,
        q_values = q_values
      )
    ))
  }

  # =========================================================================
  # MERGE RESULTS
  # =========================================================================

  # For each significant gene, combine LMM p-value with divergence effect sizes
  validation_stats <- list(
    total_genes = length(significant_genes),
    passed_lmm = 0,
    failed_missing_divergence = 0,
    other_errors = 0,
    q_values = q_values
  )

  if (verbose) {
    cat("\nMerging LMM results with divergence effect sizes...\n\n")
  }

  # Determine which column in lm_res to use for matching gene names
  # Prefer gene_name if available (set by calculate_lm_interaction), fall back to gene column
  use_gene_name_col <- "gene_name" %in% colnames(lm_res)
  
  if (verbose) {
    cat("[effect_sizes_divergence] Gene name matching strategy:\n")
    cat("  - gene_name column in lm_res:", use_gene_name_col, "\n")
    if (use_gene_name_col) {
      cat("  - lm_res$gene (first 5):", paste(head(lm_res$gene, 5), collapse=", "), "\n")
      cat("  - lm_res$gene_name (first 5):", paste(head(lm_res$gene_name, 5), collapse=", "), "\n")
    } else {
      cat("  - lm_res$gene (first 5):", paste(head(lm_res$gene, 5), collapse=", "), "\n")
    }
    cat("  - divergence gene_name (first 5):", paste(head(rd$gene_name, 5), collapse=", "), "\n")
    cat("\n")
  }

  if (verbose) {
    cat("\n[effect_sizes_divergence] MERGE STARTING\n")
    cat("  - significant_genes count:", length(significant_genes), "\n")
    cat("  - lm_res rows:", nrow(lm_res), "\n")
    cat("  - divergence rowData rows:", nrow(rd), "\n")
    cat("  - use_gene_name_col:", use_gene_name_col, "\n")
  }
  
  for (i in seq_along(significant_genes)) {
    gene_id <- significant_genes[i]
    
    # Get LMM info
    lmm_row <- lm_res[lm_res$gene == gene_id, ]
    if (nrow(lmm_row) == 0) {
      validation_stats$other_errors <- validation_stats$other_errors + 1
      next
    }

    # Determine the gene name to use for matching against divergence_results_se
    match_name <- if (use_gene_name_col && !is.na(lmm_row$gene_name[1])) {
      lmm_row$gene_name[1]
    } else {
      gene_id
    }

    p_interaction <- lmm_row$adj_p_interaction[1]
    slope_diff <- if ("slope_diff" %in% colnames(lmm_row)) {
      lmm_row$slope_diff[1]
    } else {
      NA_real_
    }

    if (verbose && i <= min(3, length(significant_genes))) {
      cat("  [Gene ", i, "] gene_id='", gene_id, "' match_name='", match_name, "'", sep="")
    }

    # Get divergence info from SE
    # Try to match by gene_name first (most reliable), then by rownames
    div_row <- rd[rd$gene_name == match_name, ]
    
    if (nrow(div_row) == 0 && !use_gene_name_col) {
      # If match_name is an ID and we have rownames, try matching rownames
      div_row <- rd[rownames(rd) == match_name, ]
    }

    if (verbose && i <= min(3, length(significant_genes))) {
      cat(" -> found ", nrow(div_row), " row(s)\n", sep="")
    }

    if (nrow(div_row) == 0) {
      validation_stats$failed_missing_divergence <- validation_stats$failed_missing_divergence + 1
      if (verbose && i > min(3, length(significant_genes))) {
        # Only show skipped messages for genes after the debug ones
        cat("  [Skipped] ", match_name, " - divergence data not found\n", sep = "")
      }
      next
    }

    # Check if any divergence estimates are available
    if (use_generic) {
      # Single q result
      if (is.na(div_row$estimate[1])) {
        validation_stats$failed_missing_divergence <- validation_stats$failed_missing_divergence + 1
        if (verbose) {
          cat("  [Skipped] ", gene_name, " - divergence estimate is NA\n", sep = "")
        }
        next
      }
      
      new_row <- data.frame(
        gene = match_name,
        p_value_interaction = p_interaction,
        slope_diff = slope_diff,
        effect_size_D = abs(div_row$estimate[1]),
        D_lower_ci = div_row$lower_ci[1],
        D_upper_ci = div_row$upper_ci[1],
        stringsAsFactors = FALSE
      )
    } else {
      # Per-q results - check if any are available
      any_valid <- FALSE
      new_row <- data.frame(
        gene = match_name,
        p_value_interaction = p_interaction,
        slope_diff = slope_diff,
        stringsAsFactors = FALSE
      )
      
      for (q_val in q_values) {
        estimate_col <- paste0("estimate_q", q_val)
        lower_col <- paste0("lower_ci_q", q_val)
        upper_col <- paste0("upper_ci_q", q_val)
        q_label <- gsub("\\.", "_", as.character(q_val))
        
        if (!is.na(div_row[[estimate_col]][1])) {
          any_valid <- TRUE
          new_row[[paste0("effect_size_D_q", q_label)]] <- abs(div_row[[estimate_col]][1])
          new_row[[paste0("D_q", q_label, "_lower_ci")]] <- div_row[[lower_col]][1]
          new_row[[paste0("D_q", q_label, "_upper_ci")]] <- div_row[[upper_col]][1]
        } else {
          # Set to NA for this q
          new_row[[paste0("effect_size_D_q", q_label)]] <- NA_real_
          new_row[[paste0("D_q", q_label, "_lower_ci")]] <- NA_real_
          new_row[[paste0("D_q", q_label, "_upper_ci")]] <- NA_real_
        }
      }
      
      if (!any_valid) {
        validation_stats$failed_missing_divergence <- validation_stats$failed_missing_divergence + 1
        if (verbose) {
          cat("  [Skipped] ", match_name, " - all divergence estimates are NA\n", sep = "")
        }
        next
      }
    }

    # Add to results
    interaction_results <- rbind(interaction_results, new_row)
    validation_stats$passed_lmm <- validation_stats$passed_lmm + 1

    if (verbose) {
      if (use_generic) {
        div_val <- abs(div_row$estimate[1])
        ci_text <- if (!is.na(div_row$lower_ci[1])) {
          sprintf(" CI=[%.4f, %.4f]", div_row$lower_ci[1], div_row$upper_ci[1])
        } else {
          ""
        }
      } else {
        # Show summary across q values
        div_vals <- sapply(q_values, function(q) {
          estimate_col <- paste0("estimate_q", q)
          div_row[[estimate_col]][1]
        })
        div_summary <- paste(sprintf("%.3f", abs(div_vals)), collapse = ", ")
        div_val <- max(abs(div_vals), na.rm = TRUE)
        ci_text <- ""
      }
      
      cat("  [SUCCESS] ", match_name, " - p=", 
          format(p_interaction, digits = 3), ", D_spectrum=[", 
          if (use_generic) format(div_val, scientific = TRUE, digits = 3)
          else div_summary, "]", ci_text, "\n", sep = "")
    }
  }

  # =========================================================================
  # SUMMARY
  # =========================================================================

  # SUMMARY: Print merge results (only if verbose)
  if (verbose) {
    cat("\n[effect_sizes_divergence] MERGE COMPLETED\n")
    cat("  - Total significant genes:", validation_stats$total_genes, "\n")
    cat("  - Passed merge:", validation_stats$passed_lmm, "\n")
    cat("  - Failed (missing divergence):", validation_stats$failed_missing_divergence, "\n")
    cat("  - Other errors:", validation_stats$other_errors, "\n")
    cat("  - interaction_results rows:", nrow(interaction_results), "\n\n")
    cat("\n**Effect Size Merge Summary:**\n\n")
    cat("- Total significant genes:", validation_stats$total_genes, "\n")
    cat("- Passed merge:", validation_stats$passed_lmm, "\n")
    cat("- Failed (missing divergence):", validation_stats$failed_missing_divergence, "\n")
    cat("- Other errors:", validation_stats$other_errors, "\n\n")

    if (nrow(interaction_results) > 0) {
      cat("**Effect Size Distribution Across q Values:**\n")
      
      if (use_generic) {
        div_col <- "effect_size_D"
        cat("- Mean D:", round(mean(interaction_results[[div_col]], na.rm = TRUE), 4), "\n")
        cat("- Median D:", round(median(interaction_results[[div_col]], na.rm = TRUE), 4), "\n")
        cat("- Range: [", 
            round(min(interaction_results[[div_col]], na.rm = TRUE), 4), ", ",
            round(max(interaction_results[[div_col]], na.rm = TRUE), 4), "]\n")
      } else {
        for (q_val in q_values) {
          q_label <- gsub("\\.", "_", as.character(q_val))
          div_col <- paste0("effect_size_D_q", q_label)
          if (div_col %in% colnames(interaction_results)) {
            valid_vals <- interaction_results[[div_col]][!is.na(interaction_results[[div_col]])]
            if (length(valid_vals) > 0) {
              cat("- q=", q_val, ": mean=", round(mean(valid_vals, na.rm = TRUE), 4),
                  ", median=", round(median(valid_vals, na.rm = TRUE), 4), "\n")
            }
          }
        }
      }
      
      cat("- Interpretation: D > 0.05 = small, D > 0.1 = medium, D > 0.2 = large\n\n")
    }
  }

  # =========================================================================
  # ENRICH RESULTS: Add per_q_pattern column
  # =========================================================================

  if (enrich_per_q_pattern && nrow(interaction_results) > 0) {
    # Extract gene names from divergence_results_se rowData and assay matrix
    div_assay <- SummarizedExperiment::assay(divergence_results_se)
    div_rd <- as.data.frame(SummarizedExperiment::rowData(divergence_results_se))

    if (nrow(div_assay) > 0 && ncol(div_assay) > 0) {
      # Map genes from divergence_results_se
      div_gene_names <- if ("gene_name" %in% colnames(div_rd)) {
        div_rd$gene_name
      } else {
        rownames(div_assay)
      }

      # Create per_q_pattern column: classify divergence patterns
      # RARE_DRIVEN = divergence higher at low q (rare isoforms drive changes)
      # ABUNDANT_DRIVEN = divergence higher at high q (abundant isoforms drive changes)  
      # BALANCED = similar divergence across diversity scales
      per_q_patterns <- character(nrow(interaction_results))
      for (i in seq_len(nrow(interaction_results))) {
        gene_name <- interaction_results$gene[i]
        gene_idx <- which(div_gene_names == gene_name)

        if (length(gene_idx) > 0) {
          # Get divergence values for this gene across q values
          divs <- div_assay[gene_idx[1], ]
          
          # Create named vector for classify_q_pattern
          # Column names in divs should be like "q_0.01", "q_0.5", "q_1.0", etc.
          per_q_patterns[i] <- .classify_q_pattern(divs)
          
          # If classification failed, return "UNCLASSIFIED"
          if (is.na(per_q_patterns[i])) {
            per_q_patterns[i] <- "UNCLASSIFIED"
          }
        }
      }
      interaction_results$per_q_pattern <- per_q_patterns
    }
  }

  return(list(
    interaction_results = interaction_results,
    validation_stats = validation_stats
  ))
}
