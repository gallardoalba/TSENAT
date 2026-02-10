# Prevent tests from creating `Rplots.pdf` in the package directory during
# non-interactive test runs by opening a temporary PDF device for the duration
# of the test suite. This device will be closed and the temp file removed when
# the R session exits.
if (!interactive()) {
    tmp_tests_pdf <- tempfile("TSENAT_tests_plot_", fileext = ".pdf")
    tryCatch({
        grDevices::pdf(tmp_tests_pdf)
    }, error = function(e) {
        # If we cannot open a pdf device, fall back to quietly doing nothing.
        message("helper-graphics: could not open temporary pdf device: ", conditionMessage(e))
    })

    # Ensure the device is closed and tempfile removed when the session ends
    cleanup_env <- new.env(parent = emptyenv())
    cleanup_env$tmp_tests_pdf <- tmp_tests_pdf
    reg.finalizer(cleanup_env, function(e) {
        tryCatch({
            while (grDevices::dev.cur() > 1) grDevices::dev.off()
        }, error = function(err) NULL)
        if (file.exists(e$tmp_tests_pdf)) unlink(e$tmp_tests_pdf)
    }, onexit = TRUE)
    # Keep the environment alive during the session
    assign(".TSENAT_test_graphics_env", cleanup_env, envir = .GlobalEnv)
}