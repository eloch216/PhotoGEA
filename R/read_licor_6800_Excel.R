# Helping function for processing user constants that were logged as rows
process_user_constant_rows_6800_excel <- function(licor_file) {
    # User constants are included as two sequential rows. The first row has a
    # character value in its first column, and the second row has NA in its
    # first column. This can be used to identify the first row in each pair.
    uc_rows <- sapply(seq_len(nrow(licor_file) - 1), function(i) {
        !is.na(licor_file[i, 1]) && is.character(licor_file[i, 1]) && is.na(licor_file[i + 1, 1])
    })

    # Get the indices of user constant rows
    uc_row_indx <- which(uc_rows)

    # Process the user constant rows if there are any
    if (length(uc_row_indx) > 0) {
        # Update the main table following each user constant row
        for (i in uc_row_indx) {
            # The first column of the first row contains the category
            uc_category <- licor_file[i, 1]

            # The second column of the first row contains the name
            uc_name <- licor_file[i, 2]

            # The second column of the second row contains the value
            uc_val <- licor_file[i + 1, 2]

            # If this column is present in the main data table, update its value
            # in this row and all following rows
            uc_exists <- uc_name %in% colnames(licor_file) &&
                licor_file[['categories']][[uc_name]] == uc_category

            if (uc_exists) {
                new_vals                           <- licor_file[, uc_name]
                new_vals[seq(i, nrow(licor_file))] <- uc_val
                licor_file[, uc_name]              <- new_vals
            }
        }

        # Remove the user constant rows from the data table
        all_uc_rows <- rep_len(FALSE, nrow(licor_file))

        all_uc_rows[uc_row_indx]     <- TRUE
        all_uc_rows[uc_row_indx + 1] <- TRUE

        licor_file <- licor_file[!all_uc_rows, , TRUE]

        # Make sure there are no row names
        rownames(licor_file$main_data) <- NULL

        # The user-constant rows may have prevented some columns from being
        # properly identified as having numeric values, so try to convert them
        # again
        for (i in seq_len(ncol(licor_file))) {
            licor_file[['main_data']][, i] <-
                try_as_numeric(licor_file[['main_data']][, i])
        }
    }

    # Return the (potentially) modified version of the exdf object
    licor_file
}

# Helping function for extracting one row of the preamble data and returning it
# as an exdf object
extract_6800_excel_preamble_row <- function(preamble_df, start_indx) {
    # Initialize a placeholder data frame
    pr <- as.data.frame(matrix(nrow = 1, ncol = ncol(preamble_df)))

    # Read the column values
    pr[1, ] <- replace_unicode(preamble_df[start_indx + 1, ])

    # Read the column names
    colnames(pr) <- replace_unicode(preamble_df[start_indx, ])

    # Remove any empty columns
    pr <- pr[!is.na(colnames(pr))]

    # Convert the data to numeric values whenever possible
    pr <- as.data.frame(
        lapply(pr, try_as_numeric),
        stringsAsFactors = FALSE
    )

    # The first column indicates the category
    pr_category <- colnames(pr)[1]
    pr[, 1] <- NULL

    # Form a data frame with the categories
    pr_category_df <- pr
    pr_category_df[1, ] <- pr_category

    # Create and return an exdf; the preamble rows do not specify units
    exdf(pr, categories = pr_category_df)
}

read_licor_6800_Excel <- function(
    file_name,
    column_name = 'obs',
    get_oxygen = TRUE,
    check_for_zero = c('A', 'gsw'),
    include_user_remark_column = TRUE,
    remove_NA_rows = TRUE,
    ...
)
{
    # Get the names of sheets in the workbook
    sheet_names <- openxlsx::getSheetNames(file_name)

    if (!'Measurements' %in% sheet_names) {
        stop(paste0(
            'A sheet named `Measurements` could not be found in Excel file `',
            file_name,
            '`'
        ))
    }

    has_remarks <- 'Remarks' %in% sheet_names

    # Read the entire first sheet of the workbook into a single data frame
    rawdata <- openxlsx::readWorkbook(
        file_name,
        sheet = 'Measurements',
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
    licor_variable_names <- make.unique(replace_unicode(rawdata[data_row, ]))

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
    colnames(licor_data)                <- licor_variable_names
    colnames(licor_variable_categories) <- licor_variable_names
    colnames(licor_variable_units)      <- licor_variable_names

    # Remove NA rows if necessary
    if (remove_NA_rows) {
        all_NA <- sapply(seq_len(nrow(licor_data)), function(i) {
            all(is.na(as.list(licor_data[i, ])))
        })
        licor_data <- licor_data[!all_NA, ]
    }

    # Get the raw preamble data
    raw_preamble <- rawdata[seq_len(data_row - 2), ]

    # Get all the preamble rows as exdf objects
    row_exdf_list <- lapply(
        seq(1, nrow(raw_preamble), by = 2),
        function(i) {extract_6800_excel_preamble_row(raw_preamble, i)}
    )

    # Combine preamble rows into one exdf object
    licor_preamble <- do.call(cbind, row_exdf_list)

    # Get the remarks, if possible
    if (has_remarks) {
        # Read the entire second sheet of the workbook to get the remarks
        rawdata_remarks <- openxlsx::readWorkbook(
            file_name,
            sheet = 'Remarks',
            colNames = FALSE,
            skipEmptyRows = FALSE,
            skipEmptyCols = FALSE
        )

        # Replace any unicode
        rawdata_remarks[, 1] <- replace_unicode(rawdata_remarks[, 1])
        rawdata_remarks[, 2] <- replace_unicode(rawdata_remarks[, 2])

        # Find the user remark rows, whose first column values are formatted
        # like HH:MM:SS
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
    } else {
        user_remarks <- data.frame(
            remark_time = character(),
            remark_value = character()
        )

        remarks <- data.frame(matrix(nrow = 1, ncol = 0))
    }

    # Create the exdf object
    exdf_obj <- exdf(
        cbind(licor_data,                licor_preamble$main_data),
        cbind(licor_variable_units,      licor_preamble$units),
        cbind(licor_variable_categories, licor_preamble$categories),
        preamble = licor_preamble$main_data,
        remarks = remarks,
        data_row = data_row,
        user_remarks = user_remarks
    )

    # Process the user constant rows
    exdf_obj <- process_user_constant_rows_6800_excel(exdf_obj)

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
