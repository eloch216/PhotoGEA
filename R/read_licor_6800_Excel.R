read_licor_6800_Excel <- function(
    file_name,
    column_name = 'obs',
    get_oxygen = TRUE,
    check_for_zero = c('A', 'gsw'),
    include_user_remark_column = TRUE,
    ...
)
{
    # Read the entire first sheet of the workbook into a single data frame
    rawdata <- openxlsx::readWorkbook(
        file_name,
        sheet = 1,
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

    # Read the entire second sheet of the workbook to get the remarks
    rawdata_remarks <- openxlsx::readWorkbook(
        file_name,
        sheet = 2,
        colNames = FALSE,
        skipEmptyRows = FALSE,
        skipEmptyCols = FALSE
    )

    # Replace any unicode
    rawdata_remarks[, 1] <- replace_unicode(rawdata_remarks[, 1])
    rawdata_remarks[, 2] <- replace_unicode(rawdata_remarks[, 2])

    # Find the user remark rows, whose first column values are formatted like
    # HH:MM:SS
    row_is_remark <- grepl(
      '^[[:digit:]]{2}:[[:digit:]]{2}:[[:digit:]]{2}$',
      rawdata_remarks[, 1]
    )

    # Extract the user remarks
    user_remarks <- rawdata_remarks[row_is_remark, ]

    colnames(user_remarks) <- c('remark_time', 'remark_value')
    rownames(user_remarks) <- NULL

    # Get the other entries in the remarks
    other_remarks <- t(rawdata_remarks[!row_is_remark, ])
    other_remarks <- data.frame(other_remarks)

    remarks           <- other_remarks[2, ]
    colnames(remarks) <- as.character(other_remarks[1, ])
    rownames(remarks) <- NULL

    # Create the exdf object
    exdf_obj <- exdf(
        licor_data,
        licor_variable_units,
        licor_variable_categories,
        preamble = cbind(remarks, licor_preamble),
        data_row = data_row,
        user_remarks = user_remarks
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

    # Add user remarks if necessary
    if (include_user_remark_column) {
        exdf_obj <- add_latest_remark(exdf_obj)
    }

    # Return the object, including oxygen information if necessary
    if (get_oxygen) {
        get_oxygen_from_preamble(exdf_obj)
    } else {
        exdf_obj
    }
}
