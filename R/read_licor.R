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
    # Define a helping function for reading one row (with column names) from an
    # Excel file, where the data row is specified and the column names are
    # assumed to be in the preceding row
    read_row_wc <- function(row_num)
    {
        row_data <- openxlsx::readWorkbook(
            file_name,
            colNames = TRUE,
            skipEmptyCols = FALSE,
            rows = c(row_num - 1, row_num)
        )
    }

    # Get the preamble data rows
    licor_preamble <- lapply(preamble_data_rows, read_row_wc)

    # Define a helping function for reading one row (without column names) from
    # an Excel file
    read_row_nc <- function(row_num)
    {
        row_data <- openxlsx::readWorkbook(
            file_name,
            colNames = FALSE,
            skipEmptyCols = FALSE,
            rows = row_num
        )
    }

    # Define a helping function for converting a row to a data frame
    row_to_data_frame <- function(rowdata)
    {
        data_frame <- data.frame(
            matrix(
                ncol = length(rowdata),
                nrow = 1
            ),
            stringsAsFactors = FALSE
        )
        data_frame[1,] <- rowdata
        return(data_frame)
    }

    # Define a helping function for reading a row, replacing any Unicode
    # characters, and converting the result to a data frame
    get_processed_row_data_frame <- function(row_num)
    {
        row_data <- read_row_nc(row_num)
        row_data <- replace_unicode(row_data)
        row_data <- row_to_data_frame(row_data)
        return(row_data)
    }

    # Get the column types, names and units, replacing any problematic Unicode
    # characters and converting the results to data frames when required
    licor_variable_categories <- get_processed_row_data_frame(variable_category_row)
    licor_variable_names <- get_processed_row_data_frame(variable_name_row)
    licor_variable_units <- get_processed_row_data_frame(variable_unit_row)

    # Get the main data
    licor_data <- openxlsx::readWorkbook(
        file_name,
        startRow = data_start_row,
        colNames = FALSE,
        skipEmptyCols = FALSE
    )

    # Convert the data to numeric values whenever possible
    licor_data <- as.data.frame(
        lapply(licor_data, try_as_numeric),
        stringsAsFactors = FALSE
    )

    # Apply column names
    colnames(licor_data) <- licor_variable_names[1,]
    colnames(licor_variable_categories) <- licor_variable_names[1,]
    colnames(licor_variable_units) <- licor_variable_names[1,]

    return(
        exdf(
            licor_data,
            licor_variable_units,
            licor_variable_categories,
            file_name = file_name,
            instrument_type = 'Licor LI-6800',
            file_type = 'plaintext',
            preamble = licor_preamble,
            preamble_data_rows = preamble_data_rows,
            variable_category_row = variable_category_row,
            variable_name_row = variable_name_row,
            variable_unit_row = variable_unit_row,
            data_start_row = data_start_row
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
