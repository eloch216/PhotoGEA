read_gasex_file <- function(
    file_name,
    timestamp_colname = NA,
    posix_options = list(),
    file_type = 'AUTO',
    instrument_type = 'AUTO',
    ...
)
{
    # Some basic input checking
    if (file_name == '') {
        stop(
            'The `file_name` input argument is an empty string. If ',
            '`file_name` was generated using `system.file`, this means ',
            'that the desired file could not be found.'
        )
    }

    if (!file.exists(file_name)) {
        stop('`', file_name, '` does not exist')
    }

    if (!is.list(posix_options)) {
        stop('posix_options must be a list')
    }

    if (length(posix_options) > 1 &&
        (is.null(names(posix_options)) || '' %in% names(posix_options))) {
        stop('all elements of posix_options must have names')
    }

    # Try to determine the file type from its name, if necessary
    file_type <- if (file_type == 'AUTO') {
        extension <- tools::file_ext(file_name)
        if (extension == '') {
            'plaintext'
        } else if (extension == 'xlsx') {
            'Excel'
        } else if (extension == 'dat') {
            'data'
        } else {
            stop(paste('Cannot automatically determine file_type for file:', file_name))
        }
    } else {
        file_type
    }

    # Try to determine the instrument type from the file type, if necessary
    instrument_type <- if (instrument_type == 'AUTO') {
        if (file_type == 'plaintext' || file_type == 'Excel') {
            'Licor LI-6800'
        } else if (file_type == 'data') {
            'CR3000'
        } else {
            stop(paste('Cannot automatically determine instrument_type for file:', file_name))
        }
    } else {
        instrument_type
    }

    # Try to load the file using the appropriate method
    gasex_exdf <- if (instrument_type == 'Licor LI-6800' && file_type == 'plaintext') {
        read_licor_6800_plaintext(file_name, ...)
    } else if (instrument_type == 'Licor LI-6800' && file_type == 'Excel') {
        read_licor_6800_Excel(file_name, ...)
    } else if (instrument_type == 'CR3000' && file_type == 'data') {
        read_cr3000(file_name, ...)
    } else {
        stop(paste('Unsupported (instrument_type file_type) option:', instrument_type, file_type))
    }

    # Make sure the timestamp column is properly interpreted
    if (!is.na(timestamp_colname)) {
        # Set the default POSIX options
        full_posix_options <- list(origin = '1970-01-01', tz = '')

        # Override any default options with user-supplied options
        for (i in seq_along(posix_options)) {
            full_posix_options[[names(posix_options)[i]]] <- posix_options[[i]]
        }

        # Add the `x` value
        full_posix_options[['x']] <- gasex_exdf$main_data[[timestamp_colname]]

        # Convert to POSIXlt
        gasex_exdf$main_data[[timestamp_colname]] <-
            do.call(as.POSIXlt, full_posix_options)
    }

    # Add a filename column to the exdf
    gasex_exdf <- set_variable(
        gasex_exdf,
        'file_name',
        category = 'read_gasex_file',
        value = file_name
    )

    # Add "extras" to the exdf
    gasex_exdf$file_name <- file_name
    gasex_exdf$file_type <- file_type
    gasex_exdf$instrument_type <- instrument_type
    gasex_exdf$timestamp_colname <- timestamp_colname

    return(gasex_exdf)
}
