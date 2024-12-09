read_licor_6800_Excel <- function(
    file_name,
    column_name = 'obs',
    get_oxygen = TRUE,
    check_for_zero = c('A', 'gsw'),
    ...
)
{
    # Read the entire workbook into a single data frame
    rawdata <- openxlsx::readWorkbook(
        file_name,
        colNames = FALSE,
        skipEmptyRows = FALSE,
        skipEmptyCols = FALSE
    )

    # Search for the column and row indices that match the column name
    data_search <- sapply(rawdata, function(x) {match(column_name, x)})
    data_col <- which(!is.na(data_search))
    data_row <- data_search[data_col]

    if (length(data_row) < 1) {
        stop(paste0(
            'A column named `',
            column_name,
            '` could not be found in file `',
            file_name,
            '`'
        ))
    }

    # Get variable names, units, and categories
    licor_variable_names <- replace_unicode(rawdata[data_row, ])

    licor_variable_units <- as.data.frame(matrix(nrow = 1, ncol = ncol(rawdata)), stringsAsFactors = FALSE)
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

    # Create the exdf object
    exdf_obj <- exdf(
        licor_data,
        licor_variable_units,
        licor_variable_categories,
        preamble = licor_preamble,
        data_row = data_row
    )

    # Check for columns that are all zero
    all_zeros <- sapply(check_for_zero, function(cn) {
        all(exdf_obj[, cn] == 0)
    })

    if (any(all_zeros)) {
        all_zero_cols <- check_for_zero[all_zeros]
        msg <- paste0(
            'The following columns in Licor 6800 Excel file `', file_name,
            '` are all zero: ', paste(all_zero_cols, collapse = ', '),
            '.\nYou may need to open the file in Excel to "calculate" its ',
            'values; type `?read_licor_6800_Excel` for more information.'
        )
        stop(msg)
    }

    # Return the object, including oxygen information if necessary
    if (get_oxygen) {
        get_oxygen_from_preamble(exdf_obj)
    } else {
        exdf_obj
    }
}
