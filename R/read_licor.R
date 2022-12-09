read_licor_6800_plaintext <- function(file_name) {
    # First read the file as a set of lines. This will allow us to find the rows
    # where the [Header] and [Data] sections begin.
    fconn <- file(file_name)
    file_lines <- readLines(fconn)
    close(fconn)

    header_indx <- which(file_lines == '[Header]')
    data_indx <- which(file_lines == '[Data]')

    # Define a helping function for reading a row as a data frame and then
    # replacing any Unicode characters
    get_processed_row_data_frame <- function(skip)
    {
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

read_licor_6800_Excel <- function(
    file_name,
    preamble_data_rows,
    variable_category_row,
    variable_name_row,
    variable_unit_row,
    data_start_row
)
{
    # Read the entire workbook into a single data frame
    rawdata <- openxlsx::readWorkbook(
        file_name,
        colNames = FALSE,
        skipEmptyRows = FALSE,
        skipEmptyCols = FALSE
    )

    # Search for the column and row indices that match 'obs'
    column_name <- 'obs'
    data_search <- sapply(rawdata, function(x) {match(column_name, x)})
    data_col <- which(!is.na(data_search))
    data_row <- data_search[data_col]

    # Get variable names, units, and categories
    licor_variable_names <- replace_unicode(rawdata[data_row, ])

    licor_variable_units <- as.data.frame(matrix(nrow = 1, ncol = ncol(rawdata)))
    licor_variable_units[1, ] <- replace_unicode(rawdata[data_row + 1, ])
    colnames(licor_variable_units) <- licor_variable_names

    licor_variable_categories <- licor_variable_units
    licor_variable_categories[1, ] <- replace_unicode(rawdata[data_row - 1, ])

    # Get the main data
    licor_data <- rawdata[seq(data_row + 2, nrow(rawdata)), ]

    # Convert the data to numeric values whenever possible
    licor_data <- as.data.frame(
        lapply(licor_data, try_as_numeric),
        stringsAsFactors = FALSE
    )

    # Apply column names
    colnames(licor_data) <- licor_variable_names
    colnames(licor_variable_categories) <- licor_variable_names
    colnames(licor_variable_units) <- licor_variable_names

    # Get the raw preamble data
    raw_preamble <- rawdata[seq_len(data_row - 2), ]

    # Define a helping function for extracting one row of the preamble data
    extract_preamble_row <- function(start_indx) {
        pr <- as.data.frame(matrix(nrow = 1, ncol = ncol(raw_preamble)))
        pr[1, ] <- replace_unicode(raw_preamble[start_indx + 1, ])
        colnames(pr) <- replace_unicode(raw_preamble[start_indx, ])
        pr[!is.na(colnames(pr))]
    }

    # Get all the rows and combine them into one data frame
    row_df_list <- lapply(
        seq(1, nrow(raw_preamble), by = 2),
        extract_preamble_row
    )
    licor_preamble <- do.call(cbind, row_df_list)

    return(
        exdf(
            licor_data,
            licor_variable_units,
            licor_variable_categories,
            file_name = file_name,
            instrument_type = 'Licor LI-6800',
            file_type = 'plaintext',
            preamble = licor_preamble,
            data_row = data_row
        )
    )
}

read_licor_file <- function(
    file_name,
    preamble_data_rows,
    variable_category_row,
    variable_name_row,
    variable_unit_row,
    data_start_row,
    timestamp_colname = NA,
    instrument_type = 'Licor LI-6800',
    file_type = 'AUTO'
)
{
    # Try to determine the file type from its name, if necessary
    file_type <- if (file_type == 'AUTO') {
        extension <- tools::file_ext(file_name)
        if (extension == '') {
            'plaintext'
        } else if (extension == 'xlsx') {
            'Excel'
        } else {
            stop(paste('Cannot automatically determine file_type for file:', file_name))
        }
    } else {
        file_type
    }

    # Try to load the file using the appropriate method
    licor_exdf <- if (instrument_type == 'Licor LI-6800' && file_type == 'plaintext') {
        read_licor_6800_plaintext(file_name)
    } else if (instrument_type == 'Licor LI-6800' && file_type == 'Excel') {
        read_licor_6800_Excel(
            file_name,
            preamble_data_rows,
            variable_category_row,
            variable_name_row,
            variable_unit_row,
            data_start_row
        )
    } else {
        stop(paste('Unsupported (instrument_type file_type) option:', instrument_type, file_type))
    }

    # Make sure the timestamp column is properly interpreted
    if (!is.na(timestamp_colname)) {
        licor_exdf$main_data[[timestamp_colname]] <- as.POSIXlt(
            licor_exdf$main_data[[timestamp_colname]],
            origin = "1970-01-01"
        )
    }

    return(licor_exdf)
}
