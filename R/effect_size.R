#' Calculate Effect Sizes for Wilcoxon Tests
#'
#' Compute effect size measures for Wilcoxon rank sum and signed rank tests.
#' Supports multiple effect size metrics including Cliff's delta, r-value,
#' and rank-biserial correlation.
#'
#' @param x Matrix with splicing diversity values or other numeric data. Ignored if \code{ts_se} is provided.
#' @param samples Character vector with length equal to number of columns in input dataset, specifying category of each sample (two groups expected).
#'   Ignored if \code{ts_se} is provided; extracted from colData automatically.
#' @param ts_se Optional SummarizedExperiment containing diversity values and sample metadata. When provided, x and samples are extracted automatically from the diversity assay.
#'   assay and colData. If \code{NULL}, \code{x} and \code{samples} must be provided.
#' @param paired Logical. If TRUE, compute paired effect sizes (signed rank test, default FALSE).
#' @param pairs Optional character vector with length equal to number of columns in input dataset, specifying pairing identifier for each sample. 
#'   When provided with \code{paired = TRUE}, samples are matched based on this 
#'   pairing information rather than column order.
#' @param res Optional \code{data.frame} with statistical test results containing columns
#'   'pvalue' and 'padj'. When provided, combines effect sizes with
#'   p-values for integrated interpretation (see Details).
#' @param nthreads Number of threads for parallel processing (default: 1).
#' @param return_summary Logical. If TRUE and res is provided, print summary combining p-values and effect sizes for top genes (default TRUE).
#'
#' @return A \code{data.frame} with columns:
#'   \describe{
#'   \item{cliffs_delta}{Cliff's delta effect size (range: [-1, 1]).}
#'   \item{r_value}{Standardized effect size Z / √N (range: [-1, 1]).}
#'   \item{rank_biserial}{Rank-biserial correlation (for unpaired tests).}
#'   \item{effect_magnitude}{Interpretation: 'negligible', 'small', 'medium', 'large'.}
#'   }
#'
#' @details
#' **Cliff's Delta:** Non-parametric effect size measuring the dominance probability (Cliff, 1993).
#'
#' \deqn{\delta = \frac{\#(x_i > y_j) - \#(x_i < y_j)}{n_x \cdot n_y}}{delta = [#(x > y) - #(x < y)] / (n_x * n_y)}
#'
#' Interpretation thresholds based on Robinson & Leach (1998) recommendations:
#' - \eqn{|\delta| < 0.147}{|delta| < 0.147}: negligible effect
#' - \eqn{0.147 \leq |\delta| < 0.330}{0.147 <= |delta| < 0.33}: small effect
#' - \eqn{0.330 \leq |\delta| < 0.474}{0.33 <= |delta| < 0.474}: medium effect
#' - \eqn{|\delta| \geq 0.474}{|delta| >= 0.474}: large effect
#'
#' **r-value:** Standardized effect size from Wilcoxon test statistic.
#'
#' \deqn{r = \frac{Z}{\sqrt{N}}}{r = Z / sqrt(N)}
#'
#' where \eqn{Z}{Z} is the standardized rank statistic and \eqn{N = n_1 + n_2}{N = n1 + n2}.
#'
#' Interpretation thresholds:
#' - \eqn{|r| < 0.1}{|r| < 0.1}: negligible effect
#' - \eqn{0.1 \leq |r| < 0.3}{0.1 <= |r| < 0.3}: small effect
#' - \eqn{0.3 \leq |r| < 0.5}{0.3 <= |r| < 0.5}: medium effect
#' - \eqn{|r| \geq 0.5}{|r| >= 0.5}: large effect
#'
#' **Rank-Biserial Correlation:** Effect size for paired tests (signed rank test).
#'
#' \deqn{r_{rb} = 1 - \frac{2R}{n(n+1)}}{r_rb = 1 - (2R) / (n(n+1))}
#'
#' where \eqn{R}{R} is the sum of ranks of positive differences and \eqn{n}{n} is the number of paired observations.
#'
#' **Integration with results:** When \code{res} is provided, the function automatically
#' combines effect sizes with p-values and displays an interpretation guide showing:
#' - Genes with both low p-values AND large effect sizes (confident discoveries)
#' - Genes with low p-values but small effect sizes (marginal biological relevance)
#' - Genes with high p-values but detectable effect sizes (underpowered detection)
#'
#' **Automatic data extraction from SummarizedExperiment:** When \code{ts_se} is provided:
#' - Diversity values are extracted from the "diversity" assay
#' - Sample group information is extracted from the "sample_type" column of colData
#'   (created by internal metadata mapping), or falls back to "condition" if available
#' - Pairing information is automatically extracted from the "sample_base" column if available
#'
#' @references
#' Cliff, N. (1993). Dominance statistics: Ordinal analyses for ordinal data.
#' Journal of Modelling in Management, 8, 53-70.
#' Reference: ES001
#'
#' Kerby, D. S. (2014). The simple difference formula: An approach to teaching
#' nonparametric correlation. Comprehensive Psychology, 3, 11.IT.3.1.
#' Reference: ES002
#'
#' Hollander, M., Wolfe, D. A., & Chicken, E. (2014).
#' Nonparametric Statistics: Theory and Methods (3rd ed.).
#' World Scientific Publishing Company.
#' Reference: NP001
#'
#' @seealso
#' See data preparation documentation for methods to add sample type information to colData.
#'
#' @export
#' @examples
#' # Create a matrix of values (3 genes x 6 samples)
#' mat <- matrix(c(
#'   0.5, 0.6, 0.55, 0.8, 0.85, 0.9,  # gene1: treatment higher
#'   0.7, 0.75, 0.72, 0.6, 0.55, 0.5  # gene2: control higher
#' ), nrow = 2, byrow = TRUE)
#' samples <- c('Control', 'Control', 'Control', 'Treatment', 'Treatment', 'Treatment')
#'
#' # Calculate effect sizes
#' es <- calculate_effect_sizes(mat, samples, paired = FALSE)
#' head(es)
calculate_effect_sizes <- function(x = NULL, samples = NULL, ts_se = NULL, paired = FALSE, 
                                    pairs = NULL, res = NULL, nthreads = 1, return_summary = TRUE) {
    # If ts_se is provided, extract data automatically
    if (!is.null(ts_se)) {
        if (!(inherits(ts_se, "SummarizedExperiment") || 
              inherits(ts_se, "RangedSummarizedExperiment"))) {
            stop("ts_se must be a SummarizedExperiment object")
        }
        
        # Extract diversity assay
        if ("diversity" %in% SummarizedExperiment::assayNames(ts_se)) {
            x <- as.matrix(SummarizedExperiment::assay(ts_se, "diversity"))
        } else {
            stop("'diversity' assay not found in ts_se")
        }
        
        # Extract samples from colData
        coldata <- SummarizedExperiment::colData(ts_se)
        if ("sample_type" %in% colnames(coldata)) {
            samples <- coldata$sample_type
        } else if ("condition" %in% colnames(coldata)) {
            samples <- coldata$condition
        } else {
            stop("'sample_type' or 'condition' column not found in colData")
        }
        
        # Extract pairs if available and not already provided
        if (is.null(pairs) && "sample_base" %in% colnames(coldata)) {
            pairs <- coldata$sample_base
        }
    }
    
    # Validate that x and samples are provided (either directly or via ts_se)
    if (is.null(x) || is.null(samples)) {
        stop("Either 'x' and 'samples' must be provided, or 'ts_se' must be a SummarizedExperiment")
    }
    
    # Determine groups
    groups <- unique(sort(samples))
    if (length(groups) != 2) {
        stop("`samples` must contain exactly two groups for effect size calculations.")
    }

    g1_idx <- which(samples == groups[1])
    g2_idx <- which(samples == groups[2])
    n1 <- length(g1_idx)
    n2 <- length(g2_idx)

    # Function to compute effect sizes for a single feature
    .effect_size_one <- function(i) {
        tryCatch({
            x1 <- as.numeric(x[i, g1_idx])
            x2 <- as.numeric(x[i, g2_idx])

            if (paired && !is.null(pairs)) {
                # Paired based on explicit pairing
                unique_pairs <- unique(pairs)
                diffs <- numeric(0)
                for (p in unique_pairs) {
                    g1_samples <- which(pairs == p & samples == groups[1])
                    g2_samples <- which(pairs == p & samples == groups[2])
                    if (length(g1_samples) == 1 && length(g2_samples) == 1) {
                        diffs <- c(diffs, x[i, g1_samples] - x[i, g2_samples])
                    }
                }
                return(.tsenat_effect_size_paired(diffs, length(diffs)))
            } else if (paired) {
                # Paired based on position
                if (n1 != n2) {
                    return(data.frame(cliffs_delta = NA_real_, r_value = NA_real_,
                        rank_biserial = NA_real_, effect_magnitude = NA_character_))
                }
                diffs <- x1 - x2
                return(.tsenat_effect_size_paired(diffs, n1))
            } else {
                # Unpaired
                return(.tsenat_effect_size_unpaired(x1, x2))
            }
        }, error = function(e) {
            data.frame(cliffs_delta = NA_real_, r_value = NA_real_, rank_biserial = NA_real_,
                effect_magnitude = NA_character_)
        })
    }

    # Apply across all features
    results <- .tsenat_bplapply(seq_len(nrow(x)), .effect_size_one, nthreads = nthreads)

    # Combine results
    combined <- do.call(rbind, results)
    rownames(combined) <- rownames(x)

    # If res provided and return_summary is TRUE, combine with p-values and print summary
    if (!is.null(res) && return_summary && nrow(combined) > 0) {
        top_n <- min(10, nrow(res))
        
        # Match gene indices
        common_genes <- intersect(rownames(combined), head(rownames(res), top_n))
        
        if (length(common_genes) > 0) {
            # Build combined results data frame
            combined_results <- data.frame(
                Gene = common_genes,
                Raw_Pvalue = res$pvalue[match(common_genes, rownames(res))],
                Adjusted_Pvalue = res$padj[match(common_genes, rownames(res))],
                Cliffs_Delta = combined[common_genes, "cliffs_delta"],
                r_value = combined[common_genes, "r_value"],
                Effect_Size = combined[common_genes, "effect_magnitude"],
                stringsAsFactors = FALSE
            )
            
            cat("\n")
            cat("Summary of Top Genes: Statistical Significance + Effect Sizes\n")
            cat("=============================================================\n\n")
            
            # Format and print combined_results with NA handling
            # Replace NA values with "NA" for printing to avoid logical operation issues
            combined_results_print <- combined_results
            for (i in seq_len(ncol(combined_results_print))) {
              combined_results_print[[i]][is.na(combined_results_print[[i]])] <- "NA"
            }
            print(combined_results_print, na.print = "NA")
            
            cat("\n\nInterpretation Guide:\n")
            cat("- Genes with LOW p-values AND large effect sizes: Confident discoveries\n")
            cat("- Genes with LOW p-values but SMALL effect sizes: Marginal biological relevance\n")
            cat("- Genes with HIGH p-values but large effect sizes: Underpowered detection\n\n")
            
            # Display effect size interpretation thresholds
            wilcoxon_guidelines <- wilcoxon_effect_size_guidelines()
            
            # Filter with NA-safe approach
            cliffs_rows <- which(!is.na(wilcoxon_guidelines$Metric) & wilcoxon_guidelines$Metric == "Cliff's Delta")
            if (length(cliffs_rows) > 0) {
              cat("Cliff's Delta Interpretation Thresholds:\n")
              cat("========================================\n")
              print(wilcoxon_guidelines[cliffs_rows, ], na.print = "NA")
            }
        }
    }

    return(combined)
}

.tsenat_cliffs_delta <- function(x1, x2) {
    # Remove NAs
    x1 <- x1[!is.na(x1)]
    x2 <- x2[!is.na(x2)]

    if (length(x1) == 0 || length(x2) == 0) {
        return(NA_real_)
    }

    # Count dominance: for each x1[i], count how many x2[j] values it exceeds
    greater <- sum(outer(x1, x2, ">"))
    less <- sum(outer(x1, x2, "<"))

    delta <- (greater - less) / (length(x1) * length(x2))
    return(delta)
}

.tsenat_r_value <- function(x1, x2) {
    # Remove NAs
    x1 <- x1[!is.na(x1)]
    x2 <- x2[!is.na(x2)]

    n1 <- length(x1)
    n2 <- length(x2)

    if (n1 == 0 || n2 == 0) {
        return(NA_real_)
    }

    # Perform Wilcoxon test to get statistic
    test_result <- tryCatch({
        wilcox.test(x1, x2, paired = FALSE, exact = FALSE)
    }, error = function(e) NULL)

    if (is.null(test_result)) {
        return(NA_real_)
    }

    # Extract U statistic and compute Z
    U <- test_result$statistic
    N <- n1 + n2

    # Expected value and standard deviation of U
    expected_U <- (n1 * n2) / 2
    var_U <- (n1 * n2 * (N + 1)) / 12
    sd_U <- sqrt(var_U)

    # Z-score
    Z <- (U - expected_U) / sd_U

    # r = Z / sqrt(N)
    r <- Z / sqrt(N)
    
    # Clamp to valid correlation range [-1, 1]
    r <- pmax(-1, pmin(1, r))

    return(r)
}

.tsenat_rank_biserial <- function(x1, x2) {
    # Remove NAs
    x1 <- x1[!is.na(x1)]
    x2 <- x2[!is.na(x2)]

    n1 <- length(x1)
    n2 <- length(x2)
    N <- n1 + n2

    if (n1 == 0 || n2 == 0) {
        return(NA_real_)
    }

    # Rank all values combined
    all_vals <- c(x1, x2)
    ranks <- rank(all_vals)

    # Sum of ranks for first group
    R1 <- sum(ranks[1:n1])

    # Rank-biserial formula: r_rb = (R1 - n1*(n1+1)/2) / (n1*n2)
    # This is the standard Mann-Whitney rank-biserial correlation
    r_rb <- (R1 - n1 * (n1 + 1) / 2) / (n1 * n2)

    return(r_rb)
}

.tsenat_interpret_effect_magnitude <- function(delta, type = "delta") {
    abs_delta <- abs(delta)

    if (is.na(abs_delta)) {
        return(NA_character_)
    }

    if (type == "delta") {
        # Cliff's delta interpretation thresholds
        if (abs_delta < 0.147) return("negligible")
        if (abs_delta < 0.33) return("small")
        if (abs_delta < 0.474) return("medium")
        return("large")
    } else if (type == "r") {
        # r-value interpretation thresholds
        if (abs_delta < 0.1) return("negligible")
        if (abs_delta < 0.3) return("small")
        if (abs_delta < 0.5) return("medium")
        return("large")
    } else {
        return(NA_character_)
    }
}

.tsenat_effect_size_unpaired <- function(x1, x2) {
    cliffs_delta <- .tsenat_cliffs_delta(x1, x2)
    r_value <- .tsenat_r_value(x1, x2)
    rank_biserial <- .tsenat_rank_biserial(x1, x2)

    # Interpretation based on Cliff's delta
    magnitude <- .tsenat_interpret_effect_magnitude(cliffs_delta, type = "delta")

    data.frame(
        cliffs_delta = cliffs_delta,
        r_value = r_value,
        rank_biserial = rank_biserial,
        effect_magnitude = magnitude,
        stringsAsFactors = FALSE
    )
}

.tsenat_effect_size_paired <- function(diffs, n) {
    diffs <- diffs[!is.na(diffs)]

    if (length(diffs) == 0) {
        return(data.frame(
            cliffs_delta = NA_real_,
            r_value = NA_real_,
            rank_biserial = NA_real_,
            effect_magnitude = NA_character_
        ))
    }

    # For paired data, use signed-rank test for r-value
    test_result <- tryCatch({
        wilcox.test(diffs, mu = 0, paired = FALSE, exact = FALSE)
    }, error = function(e) NULL)

    r_value <- NA_real_
    if (!is.null(test_result)) {
        U <- test_result$statistic
        N <- length(diffs)
        expected_U <- N * (N + 1) / 4
        var_U <- (N * (N + 1) * (2 * N + 1)) / 24
        sd_U <- sqrt(var_U)
        Z <- (U - expected_U) / sd_U
        r_value <- Z / sqrt(N)
        # Clamp to valid correlation range [-1, 1]
        r_value <- pmax(-1, pmin(1, r_value))
    }

    # Cliff's delta for paired data: use 0 as reference (median of differences)
    cliffs_delta <- NA_real_
    if (length(diffs) > 0) {
        pos <- sum(diffs > 0)
        neg <- sum(diffs < 0)
        cliffs_delta <- (pos - neg) / length(diffs)
    }

    # Rank-biserial for paired: simplified version
    rank_biserial <- NA_real_
    if (length(diffs) > 0) {
        abs_diffs <- abs(diffs)
        ranks <- rank(abs_diffs)
        rank_pos_diffs <- sum(ranks[diffs > 0])
        R <- rank_pos_diffs
        n_valid <- length(diffs)
        rank_biserial <- (2 * R) / (n_valid * (n_valid + 1)) - 1
    }

    # Interpretation
    magnitude <- .tsenat_interpret_effect_magnitude(cliffs_delta, type = "delta")

    data.frame(
        cliffs_delta = cliffs_delta,
        r_value = r_value,
        rank_biserial = rank_biserial,
        effect_magnitude = magnitude,
        stringsAsFactors = FALSE
    )
}

#' Effect Size Interpretation Guidelines
#'
#' @description
#' Provides standardized interpretation thresholds for effect size metrics
#' (Cliff's delta and r-value) based on established statistical guidelines.
#' This is an internal utility function.
#'
#' @return Data frame with effect size interpretation thresholds for Wilcoxon test metrics.
#'   Columns include: Metric (effect size measure name), Magnitude (interpretation level),
#'   Range (threshold range), and Description (practical meaning).
#'
#' @references
#' Cliff, N. (1993). Dominance statistics: Ordinal analyses for ordinal data.
#' Journal of Modelling in Management, 8, 53-70.
#' Reference: ES001
#'
#' Kerby, D. S. (2014). The simple difference formula: An approach to teaching
#' nonparametric correlation. Comprehensive Psychology, 3, 11.IT.3.1.
#' Reference: ES002
#'
#' @noRd
wilcoxon_effect_size_guidelines <- function() {
    cliffs_thresholds <- data.frame(
        Metric = rep("Cliff's Delta", 4),
        Magnitude = c("Negligible", "Small", "Medium", "Large"),
        Range = c("|δ| < 0.147", "0.147 <= |δ| < 0.33", "0.33 <= |δ| < 0.474", "|δ| >= 0.474"),
        Description = c(
            "Effect not meaningful",
            "Small practical effect",
            "Medium practical effect",
            "Large practical effect"
        ),
        stringsAsFactors = FALSE
    )

    r_thresholds <- data.frame(
        Metric = rep("r-value", 4),
        Magnitude = c("Negligible", "Small", "Medium", "Large"),
        Range = c("|r| < 0.1", "0.1 <= |r| < 0.3", "0.3 <= |r| < 0.5", "|r| >= 0.5"),
        Description = c(
            "Effect not meaningful",
            "Small practical effect",
            "Medium practical effect",
            "Large practical effect"
        ),
        stringsAsFactors = FALSE
    )

    rbind(cliffs_thresholds, r_thresholds)
}


