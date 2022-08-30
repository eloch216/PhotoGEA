choose_input_files <- function() {
    if (!interactive()) {
        stop("The `choose_input_files` function is only available in interactive R sessions")
    }

    single_file_path <- file.choose()
    extension <- gsub('^.+\\.', '', single_file_path)
    directory <- dirname(single_file_path)

    list.files(
        path = directory,
        pattern = paste0('\\.', extension, '$'),
        full.names = TRUE
    )
}
