choose_input_licor_files <- function()
{
    if (!interactive() | .Platform$OS.type != "windows") {
        stop(
            paste(
                "The `choose_input_licor_files` function is only available in",
                "interactive R sessions running in MS Windows"
            )
        )
    }

    utils::choose.files(
        default = "",
        caption = "Select Licor Excel input files",
        multi = TRUE,
        filters = matrix(
            c(
                "Excel files (*.xlsx)", "All files (*.*)",
                "*.xlsx", "*.*"
            ),
            ncol = 2
        ),
        index = 1
    )
}
