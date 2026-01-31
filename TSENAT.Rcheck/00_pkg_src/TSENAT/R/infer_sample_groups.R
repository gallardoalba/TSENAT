## Helper to infer sample group labels from sample names
## Supports underscore-suffixed labels (e.g. Sample_N, Sample_T)
## and TCGA barcodes (maps '01' -> 'Tumor', '11' -> 'Normal').
#' Infer sample group from sample names
#'
#' @param sample_names Character vector of sample names.
#' @return Character vector of group labels (e.g. 'Normal'/'Tumor') or NA.
#' @export
infer_sample_group <- function(sample_names) {
  if (is.null(sample_names)) return(rep(NA_character_, 0))

  # If names contain underscores, assume suffix after last underscore is group
  if (any(grepl("_", sample_names))) {
    groups <- vapply(sample_names, function(s) {
      if (grepl("_", s)) sub(".*_", "", s) else NA_character_
    }, character(1))
    # map common short labels to Normal/Tumor
    groups <- ifelse(toupper(groups) %in% c("N", "NORMAL"), "Normal",
                     ifelse(toupper(groups) %in% c("T", "TUMOR"), "Tumor", groups))
    return(unname(groups))
  }

  # TCGA barcode detection: extract two-digit sample type code like '01' or '11'
  if (any(grepl("^TCGA-", sample_names))) {
    codes <- sub(".*-(\\d{2})[A-Z]?.*", "\\1", sample_names)
    # if extraction failed for some entries, set NA
    codes[!grepl("^\\d{2}$", codes)] <- NA
    groups <- ifelse(codes == "11", "Normal",
                     ifelse(codes == "01", "Tumor", NA_character_))
    return(unname(groups))
  }

  # default: unable to infer
  return(rep(NA_character_, length(sample_names)))
}
