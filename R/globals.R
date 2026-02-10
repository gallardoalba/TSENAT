# # Declare package global variables to satisfy R CMD check NOTES
if (getRversion() >= "2.15.1") {
    utils::globalVariables(c(
        "genes", "xval", "x", "y", "padj", "padj_num", "padj_clean",
        "label_flag", "sample_q", "qnum", "significant"
    ))
}
