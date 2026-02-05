# # Declare package global variables to satisfy R CMD check NOTES
if (getRversion() >= "2.15.1") {
    utils::globalVariables(c(
        "genes",
        "xval",
        "padj_num",
        "label_flag",
        "sample_q",
        "qnum"
    ))
}
