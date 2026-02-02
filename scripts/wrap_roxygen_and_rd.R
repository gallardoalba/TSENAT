# Wrap roxygen comment blocks in R/ files and wrap prose in man/*.Rd to 80 chars
if (!requireNamespace("stringr", quietly = TRUE)) install.packages("stringr", repos = "https://cloud.r-project.org")

wrap_roxygen_block <- function(block_lines, width = 80) {
    # block_lines include leading "#'"
    # remove leading "#'" and possible leading space
    txt <- gsub("^#'\\s?", "", block_lines)
    txt <- paste(txt, collapse = " ")
    txt <- gsub("\\s+", " ", trimws(txt))
    wrapped <- stringr::str_wrap(txt, width = width - 3)
    out <- sapply(wrapped, function(ln) paste0("#' ", ln))
    out
}

wrap_rd_prose <- function(lines, width = 80) {
    out <- character(0)
    buf <- character(0)
    flush_buf <- function() {
        if (length(buf) == 0) {
            return()
        }
        txt <- paste(gsub("\\s+", " ", trimws(buf)), collapse = " ")
        wrapped <- stringr::str_wrap(txt, width = width)
        out <<- c(out, wrapped)
        buf <<- character(0)
    }
    for (ln in lines) {
        # keep directive lines starting with \\ or % untouched
        if (grepl("^\\\\|^%|^\\\\begin|^\\\\end", ln) || grepl("^\\\\[A-Za-z]", ln)) {
            flush_buf()
            out <- c(out, ln)
        } else if (nzchar(trimws(ln))) {
            buf <- c(buf, ln)
        } else {
            flush_buf()
            out <- c(out, "")
        }
    }
    flush_buf()
    out
}

# Process R/ roxygen comments
rfiles <- list.files("R", pattern = "\\.R$", full.names = TRUE)
for (f in rfiles) {
    lines <- readLines(f, warn = FALSE)
    out <- character(0)
    i <- 1
    changed <- FALSE
    while (i <= length(lines)) {
        if (grepl("^#'", lines[i])) {
            # collect roxygen block
            j <- i
            block <- character(0)
            while (j <= length(lines) && grepl("^#'", lines[j])) {
                block <- c(block, lines[j])
                j <- j + 1
            }
            wrapped <- wrap_roxygen_block(block, width = 80)
            if (!identical(block, wrapped)) changed <- TRUE
            out <- c(out, wrapped)
            i <- j
        } else {
            out <- c(out, lines[i])
            i <- i + 1
        }
    }
    if (changed) writeLines(out, f)
}

# Process man/*.Rd files
rdfiles <- list.files("man", pattern = "\\.Rd$", full.names = TRUE)
for (f in rdfiles) {
    lines <- readLines(f, warn = FALSE)
    out <- wrap_rd_prose(lines, width = 80)
    writeLines(out, f)
}

cat("Wrapping of roxygen and Rd files completed.\n")
