PhotoGEA_example_file_path <- function(example_file_name) {
    # Check inputs
    if (length(example_file_name) > 1 || !is.character(example_file_name)) {
        stop('example_file_name must be a single file name string')
    }

    # Try to find the file
    file_path <-
        system.file('extdata', example_file_name, package = 'PhotoGEA')

    # Check if it exists
    if (file_path == '') {
        msg <- paste0(
            'The PhotoGEA package does not include an example file called `',
            example_file_name,
            '`. If you are loading your own data files, do not use the ',
            '`PhotoGEA_example_file_path` function. Type ',
            '`?PhotoGEA_example_file_path` for more information.'
        )
        stop(msg)
    }

    # Return the path
    return(file_path)
}
