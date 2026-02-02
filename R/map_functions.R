#' Map external coldata into a SummarizedExperiment   This helper maps an
#' external `coldata` table (with sample IDs and a  condition/label column) into
#' a `SummarizedExperiment` produced by  `calculate_diversity()`. It assigns a
#' `sample_type` column to  `colData(ts_se)` and falls back to
#' `infer_sample_group()` for  unmapped samples.
#' @param ts_se A `SummarizedExperiment` object with assay columns named
#' possibly including `_q=` suffixes.
#' @param coldata A data.frame with sample metadata or `NULL`.
#' @param coldata_sample_col Name of the column in `coldata` with sample IDs.
#' @param coldata_condition_col Name of the column in `coldata` with condition/labels.
#' @return The input `ts_se` with `colData(ts_se)$sample_type` set when possible.
#' @export
#' @examples
#' data("tcga_brca_luma_dataset", package = "TSENAT")  rc <-
#' as.matrix(tcga_brca_luma_dataset[1:20, -1, drop = FALSE])  gs <-
#' tcga_brca_luma_dataset$genes[1:20]  se <- calculate_diversity(rc, gs, q =
#' 0.1, norm = TRUE)  sample_names <- sub("_q=.*", "",
#' colnames(SummarizedExperiment::assay(se)))  coldata_df <- data.frame(  Sample
#' = sample_names,  Condition = rep(c("A", "B"), length.out = ncol(se))  )
#' map_coldata_to_se(se, coldata_df)
map_coldata_to_se <- function(ts_se, coldata, coldata_sample_col = "Sample", coldata_condition_col = "Condition") {
  if (is.null(coldata)) {
    return(ts_se)
  }
  if (!is.data.frame(coldata) || !all(c(coldata_sample_col, coldata_condition_col) %in% colnames(coldata))) {
    return(ts_se)
  }
  sample_base_names <- sub("_q=.*", "", colnames(SummarizedExperiment::assay(ts_se)))
  st_map <- setNames(as.character(coldata[[coldata_condition_col]]), as.character(coldata[[coldata_sample_col]]))
  sample_types <- unname(st_map[sample_base_names])
  missing_idx <- which(is.na(sample_types))
  if (length(missing_idx) > 0) {
    sample_types[missing_idx] <- infer_sample_group(sample_base_names[missing_idx], coldata = coldata)
  }
  SummarizedExperiment::colData(ts_se)$sample_type <- sample_types
  return(ts_se)
}

# Helper to map sample group labels from sample names.
# Supports underscore-suffixed labels (e.g. Sample_N, Sample_T) and
# TCGA barcodes. Mapping from tokens to group labels is configurable
# via `suffix_map` and `tcga_map`; when no mapping is provided the
# raw suffix token or TCGA two-digit code is returned.
#' Infer sample group from sample names
#' @param sample_names Character vector of sample names.
#' @param suffix_sep Character separator to detect suffix groups (default "_").
#' @param suffix_map Named character vector mapping suffix tokens (case-
#' insensitive)  to group labels, e.g. c(N = "Normal", T = "Tumor").
#' @param tcga_map Named character vector mapping TCGA two-digit codes to group
#' labels, e.g. c("01" = "Tumor", "11" = "Normal").
#' @param coldata Optional data.frame or named vector providing mapping from
#' sample names to conditions. If a data.frame, specify `coldata_sample_col`
#' and `coldata_condition_col` for the relevant columns.
#' @param coldata_sample_col Column name in `coldata` indicating sample IDs
#' (default "Sample").
#' @param coldata_condition_col Column name in `coldata` for condition/label
#' (default "Condition").
#' @param prefer_suffix Logical; if TRUE, prefer suffix-based inference when
#' both patterns match.
#' @param default Character scalar returned when no mapping applies (default
#' NA_character_).
#' @return Character vector of group labels (or the `default` value) with same
#' length as `sample_names`.
#' @examples
#' # Basic usage: returns raw suffix tokens or TCGA two-digit codes when  # no
#' mapping is supplied  infer_sample_group(c("S1_N", "TCGA-XX-01A"))   # Provide
#' a suffix_map to translate tokens  #' infer_sample_group(c("S1_N", "S2_T"),
#' suffix_map = c(N = "Normal", T =  #' "Tumor"))   # Provide a TCGA mapping and
#' prefer TCGA codes over suffixes  tcga_map <- c("01" = "Tumor")
#' infer_sample_group(c("TCGA-XX-01A", "Sample_N"),  tcga_map = tcga_map,
#' prefer_suffix = FALSE  )
#' @export
infer_sample_group <- function(sample_names,
                               suffix_sep = "_",
                               suffix_map = NULL,
                               tcga_map = NULL,
                               coldata = NULL,
                               coldata_sample_col = "Sample",
                               coldata_condition_col = "Condition",
                               prefer_suffix = TRUE,
                               default = NA_character_) {
  if (is.null(sample_names)) {
    return(character(0))
  }
  n <- length(sample_names)
  res <- rep(default, n)

  # Helper to map tokens using a named map (case-insensitive on names)
  map_token <- function(token, map) {
    if (is.na(token) || token == "") {
      return(NA_character_)
    }
    if (is.null(map) || length(map) == 0) {
      return(NA_character_)
    }
    nms <- toupper(names(map))
    toku <- toupper(token)
    ix <- match(toku, nms)
    if (!is.na(ix)) {
      return(unname(map[ix]))
    }
    return(NA_character_)
  }

  # If coldata mapping is supplied and is usable, prefer its mapping
  if (!is.null(coldata)) {
    if (is.data.frame(coldata) && all(c(coldata_sample_col, coldata_condition_col) %in% colnames(coldata))) {
      map <- setNames(as.character(coldata[[coldata_condition_col]]), as.character(coldata[[coldata_sample_col]]))
      mapped <- unname(map[sample_names])
      # For entries that map successfully, fill results
      ok <- !is.na(mapped) & mapped != ""
      if (any(ok)) res[ok] <- mapped[ok]
    } else if (is.vector(coldata) && !is.null(names(coldata))) {
      mapped <- unname(coldata[sample_names])
      ok <- !is.na(mapped) & mapped != ""
      if (any(ok)) res[ok] <- mapped[ok]
    }
  }

  # Suffix-based detection
  if (!is.null(suffix_sep) && nzchar(suffix_sep) && any(grepl(suffix_sep, sample_names))) {
    groups <- vapply(sample_names, function(s) {
      if (grepl(suffix_sep, s, fixed = TRUE)) {
        tail <- sub(paste0(".*", suffix_sep), "", s)
        mapped <- map_token(tail, suffix_map)
        if (!is.na(mapped)) {
          return(mapped)
        }
        # if no mapping provided, return the raw suffix token
        return(tail)
      }
      NA_character_
    }, character(1), USE.NAMES = FALSE)
    if (prefer_suffix) {
      res[!is.na(groups)] <- groups[!is.na(groups)]
    } else {
      # Only fill res where not already assigned
      to_fill <- which(is.na(res) & !is.na(groups))
      res[to_fill] <- groups[to_fill]
    }
  }

  # TCGA barcode detection: extract two-digit sample type code like '01' or '11'
  has_tcga <- grepl("^TCGA-", sample_names)
  if (any(has_tcga)) {
    codes <- rep(NA_character_, n)
    codes[has_tcga] <- sub(".*-(\\d{2})[A-Z]?.*", "\\1", sample_names[has_tcga])
    codes[!grepl("^\\d{2}$", codes)] <- NA_character_
    groups_tcga <- vapply(codes, function(c) {
      if (is.na(c)) {
        return(NA_character_)
      }
      mapped <- map_token(c, tcga_map)
      if (!is.na(mapped)) {
        return(mapped)
      }
      # if no mapping, return raw TCGA code
      return(c)
    }, character(1), USE.NAMES = FALSE)
    # Fill in results where still default (or based on preference)
    if (prefer_suffix) {
      to_fill <- which(is.na(res) | (!is.na(default) & res == default))
    } else {
      to_fill <- which(is.na(res) | (!is.na(default) & res == default))
    }
    if (length(to_fill) > 0) res[to_fill] <- groups_tcga[to_fill]
  }

  # Ensure length and replace any remaining NA with default
  res[is.na(res)] <- default
  return(unname(res))
}
