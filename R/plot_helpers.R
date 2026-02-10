# Internal plot helpers

.tsenat_format_label <- function(lbl) {
    if (is.null(lbl)) return(NULL)
    s <- gsub("_", " ", lbl)
    s <- gsub("\\s+", " ", s)
    s <- trimws(s)
    s <- tolower(s)
    if (nchar(s) == 0) return(s)
    if (nchar(s) == 1) return(toupper(s))
    paste0(toupper(substr(s, 1, 1)), substr(s, 2, nchar(s)))
}

.tsenat_prepare_ma_plot_df <- function(df, fold_col, mean_cols, x_label, y_label) {
    # Detect x-axis values
    if (length(mean_cols) >= 2) {
        xvals <- rowMeans(df[, mean_cols[seq_len(2)], drop = FALSE], na.rm = TRUE)
        x_label <- x_label %||% paste0(mean_cols[1], " vs ", mean_cols[2])
    } else if (length(mean_cols) == 1) {
        xvals <- as.numeric(df[[mean_cols[1]]])
        x_label <- x_label %||% mean_cols[1]
    } else if ("mean" %in% colnames(df)) {
        xvals <- as.numeric(df$mean)
        x_label <- x_label %||% "Mean"
    } else {
        xvals <- seq_len(nrow(df))
        x_label <- x_label %||% "Index"
    }

    yvals <- as.numeric(df[[fold_col]])

    padj_candidates <- c("adjusted_p_values", "adj_p_value", "adj_p", "padj", "p.adjust")
    padj_col <- intersect(padj_candidates, colnames(df))
    padj_col <- if (length(padj_col)) padj_col[1] else NULL

    padj <- if (!is.null(padj_col)) as.numeric(df[[padj_col]]) else rep(1, length(yvals))
    padj[is.na(padj)] <- 1

    sig_flag <- ifelse(abs(yvals) > 0 & padj < 0.05, "significant", "non-significant")

    plot_df <- data.frame(genes = df$genes, x = xvals, y = yvals, padj = padj, significant = sig_flag, stringsAsFactors = FALSE)

    list(plot_df = plot_df, x_label = x_label, y_label = y_label)
}

.tsenat_prepare_volcano_df <- function(diff_df, x_col = NULL, padj_col = "adjusted_p_values", label_thresh = 0.1, padj_thresh = 0.05, title = NULL) {
    df <- as.data.frame(diff_df)
    cn <- colnames(df)

    # Auto-detect x-axis column if not specified
    if (is.null(x_col)) {
        diff_cols <- grep("_difference$", cn, value = TRUE, ignore.case = TRUE)
        if (length(diff_cols) > 0) {
            x_col <- diff_cols[1]
        } else {
            numeric_cols <- vapply(df, is.numeric, logical(1))
            p_cols <- grep("p_value|p.value", cn, ignore.case = TRUE)
            numeric_cols[p_cols] <- FALSE
            if (any(numeric_cols)) {
                x_col <- cn[which(numeric_cols)[1]]
            } else {
                stop("Could not find suitable column for x-axis. Specify 'x_col' explicitly.")
            }
        }
    }

    # Verify columns
    if (!(x_col %in% cn)) stop(sprintf("Column '%s' not found in diff_df", x_col))
    if (!(padj_col %in% cn)) stop(sprintf("Column '%s' not found in diff_df", padj_col))

    df$xval <- as.numeric(df[[x_col]])
    df$padj <- as.numeric(df[[padj_col]])
    df$padj[is.na(df$padj)] <- 1
    df$padj[df$padj <= 0] <- .Machine$double.xmin

    df$significant <- ifelse(abs(df$xval) >= label_thresh & df$padj < padj_thresh, "significant", "non-significant")
    df <- df[is.finite(df$xval) & is.finite(df$padj), ]

    if (nrow(df) == 0) stop("No valid points to plot")

    metric_label <- if (grepl("median", x_col, ignore.case = TRUE)) "Median" else if (grepl("mean", x_col, ignore.case = TRUE)) "Mean" else "Value"

    title_use <- title %||% "Volcano plot: fold-change vs significance"

    x_label_formatted <- .tsenat_format_label(x_col)
    padj_label_formatted <- .tsenat_format_label(padj_col)

    list(df = df, x_col = x_col, padj_col = padj_col, x_label_formatted = x_label_formatted, padj_label_formatted = padj_label_formatted, title_use = title_use)
}
