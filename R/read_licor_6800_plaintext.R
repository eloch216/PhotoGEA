read_licor_6800_plaintext <- function(file_name, ...) {
    # First read the file as a set of lines. This will allow us to find the rows
    # where the [Header] and [Data] sections begin.
    fconn <- file(file_name)
    file_lines <- readLines(fconn)
    close(fconn)

    header_indx <- which(file_lines == '[Header]')
    data_indx <- which(file_lines == '[Data]')

    # Define a helping function for reading a row as a data frame and then
    # replacing any Unicode characters
    get_processed_row_data_frame <- function(skip) {
        row_data <- utils::read.delim(
            file_name,
            skip = skip,
            nrows = 1,
            header = FALSE
        )

        row_data[1, ] <- replace_unicode(row_data[1, ])

        row_data
    }

    # Now we can read the column names, categories, units, and values from the
    # [Data] section of the file. Make sure to replace any Unicode characters in
    # the column names, units, and categories.
    licor_variable_names <- get_processed_row_data_frame(data_indx + 1)
    licor_variable_units <- get_processed_row_data_frame(data_indx + 2)
    licor_variable_categories <- get_processed_row_data_frame(data_indx)

    licor_data <- utils::read.delim(
        file_name,
        skip = data_indx + 3,
        header = FALSE
    )

    # Convert the data to numeric values whenever possible
    licor_data <- as.data.frame(
        lapply(licor_data, try_as_numeric),
        stringsAsFactors = FALSE
    )

    # Apply column names
    colnames(licor_variable_units) <- licor_variable_names[1, ]
    colnames(licor_variable_categories) <- licor_variable_names[1, ]
    colnames(licor_data) <- licor_variable_names[1, ]

    # We can also read the header information
    preamble_raw <- utils::read.delim(
        file_name,
        skip = 1,
        nrows = data_indx - header_indx,
        header = FALSE
    )

    licor_preamble <- stats::setNames(
        as.data.frame(t(preamble_raw[, 2])),
        preamble_raw[, 1]
    )

    colnames(licor_preamble) <- replace_unicode(colnames(licor_preamble))
    licor_preamble[1, ] <- replace_unicode(licor_preamble[1, ])

    # Return an exdf object with information from the file
    return(
        exdf(
            licor_data,
            licor_variable_units,
            licor_variable_categories,
            file_name = file_name,
            instrument_type = 'Licor LI-6800',
            file_type = 'plaintext',
            preamble = licor_preamble,
            header_indx = header_indx,
            data_indx = data_indx
        )
    )
}
