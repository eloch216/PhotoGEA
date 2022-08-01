choose_input_tdl_files <- function()
{
    if (!interactive() | .Platform$OS.type != "windows") {
        stop(
            paste(
                "The `choose_input_tdl_files` function is only available in",
                "interactive R sessions running in MS Windows"
            )
        )
    }

    utils::choose.files(
        default = "",
        caption = "Select TDL input files",
        multi = TRUE,
        filters = matrix(
            c(
                "TDL files (*.dat)", "All files (*.*)",
                "*.dat", "*.*"
            ),
            ncol = 2
        ),
        index = 1
    )
}
