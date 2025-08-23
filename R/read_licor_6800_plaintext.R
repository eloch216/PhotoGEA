# Helping function for processing user constants that were logged as rows
process_user_constant_rows_6800_plaintext <- function(licor_file) {
    # In user constant rows, the first column has a text value formatted like
    # "category:name". This can be used to identify the rows.
    uc_rows <- sapply(seq_len(nrow(licor_file)), function(i) {
        grepl('^[[:alnum:]]+:', licor_file[i, 1])
    })

    # Get the indices of user constant rows
    uc_row_indx <- which(uc_rows)

    # Process the user constant rows if there are any
    if (length(uc_row_indx) > 0) {
        # Update the main table following each user constant row
        for (i in uc_row_indx) {
            # The second column contains the user constant value
            uc_val  <- licor_file[i, 2]

            # Extract the name and category from the first column
            uc_info     <- strsplit(licor_file[i, 1], ':')[[1]]
            uc_name     <- uc_info[2]
            uc_category <- uc_info[1]

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
        licor_file <- licor_file[!uc_rows, , TRUE]

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

# Helping function for adding user remarks to the main data table (also used
# for LI-6800 Excel files)
add_latest_remark <- function(licor_file) {
  # Get the user remarks and order them by time
  user_remarks <- licor_file$user_remarks

  user_remarks[, 'remark_time'] <- as.POSIXct(
    user_remarks[, 'remark_time'],
    format = '%H:%M:%S'
  )

  user_remarks <- user_remarks[order(user_remarks[, 'remark_time']), ]

  # Get the time values from the main data
  data_times <- as.POSIXct(
    licor_file[, 'hhmmss'],
    format = '%H:%M:%S'
  )

  # Initialize the user_remark column and fill it in
  licor_file <- set_variable(
    licor_file,
    'user_remark',
    category = 'add_latest_remark'
  )

  for (i in seq_len(nrow(user_remarks))) {
    time_i   <- user_remarks[i, 'remark_time']
    remark_i <- user_remarks[i, 'remark_value']

    times_to_fill <- !is.na(data_times) & data_times > time_i

    licor_file[times_to_fill, 'user_remark'] <- remark_i
  }

  licor_file
}

# Helping function for extracting user remarks
extract_user_remark_table <- function(file_lines, is_remark) {
  # Get the line contents
  remark_lines <- file_lines[is_remark]

  # Split the remarks into a time and a value
  split_remark_lines <- strsplit(remark_lines, '\t')

  remark_times <- sapply(split_remark_lines, function(x) {
    x[1]
  })

  remark_values <- sapply(split_remark_lines, function(x) {
    paste(x[seq(2, length(x))], collapse = ' ')
  })

  # Replace unicode
  remark_times  <- replace_unicode(remark_times)
  remark_values <- replace_unicode(remark_values)

  # Return as a data frame
  data.frame(
    remark_time = remark_times,
    remark_value = remark_values,
    stringsAsFactors = FALSE
  )
}

# Helping function for converting raw lines to a data frame
licor_6800_lines_to_df <- function(file_lines, rows) {
  split_res <- strsplit(file_lines[rows], '\t')
  lengths   <- sapply(split_res, length)
  chunk_res <- matrix(nrow = length(rows), ncol = max(lengths))

  for (i in seq_along(rows)) {
    rowdata <- split_res[[i]]
    chunk_res[i, seq_along(rowdata)] <- rowdata
  }

  as.data.frame(chunk_res, stringsAsFactors = FALSE)
}

# Helping function for extracting the "preamble" and "remarks" from the header
# info; the preamble is returned as an exdf, and the remarks as a data frame
process_6800_plaintext_header <- function(header_df) {
    # The "preamble" lines begin with "category:name"
    preamble_indx <- grepl('^[[:alnum:]]+:', colnames(header_df))
    preamble      <- header_df[1, preamble_indx]

    # Extract the preamble names and categories
    preamble_split <- strsplit(colnames(preamble), ':')
    preamble_cat   <- sapply(preamble_split, function(x) {x[1]})
    preamble_names <- sapply(preamble_split, function(x) {x[2]})

    # Rename the preamble columns
    colnames(preamble) <- preamble_names

    # Assemble the preamble categories
    preamble_categories      <- preamble
    preamble_categories[1, ] <- preamble_cat

    # Create the preamble exdf
    preamble_exdf <- exdf(preamble, categories = preamble_categories)

    # Remove any row names that may have been added
    rownames(preamble_exdf$main_data)  <- NULL
    rownames(preamble_exdf$units)      <- NULL
    rownames(preamble_exdf$categories) <- NULL

    # The remaining lines in the header are "remarks" lines
    remarks <- header_df[1, !preamble_indx]

    # Return both tables
    list(
      preamble_exdf = preamble_exdf,
      remarks = remarks
    )
}

# Main function for reading plaintext LI-6800 log files
read_licor_6800_plaintext <- function(
    file_name,
    get_oxygen = TRUE,
    include_user_remark_column = TRUE,
    remove_NA_rows = TRUE,
    ...
)
{
    # First read the file as a set of lines
    fconn <- file(file_name)
    file_lines <- readLines(fconn)
    close(fconn)

    # Find the remark lines, which begin with HH:MM:SS followed by a tab
    line_is_remark <- grepl(
      '^[[:digit:]]{2}:[[:digit:]]{2}:[[:digit:]]{2}\t',
      file_lines
    )

    # Get the user remarks as a table
    user_remarks <- extract_user_remark_table(file_lines, line_is_remark)

    # Remove the remark lines so they don't appear in the main data
    file_lines <- file_lines[!line_is_remark]

    # Find the rows where the [Header] and [Data] sections begin
    header_indices <- which(file_lines == '[Header]')
    data_indices <- which(file_lines == '[Data]')

    # Check basic file properties
    if (length(file_lines) < 1) {
        stop(
            'file ', file_name, ' is malformed: must contain at least one line'
        )
    }

    if (length(header_indices) < 1 || length(data_indices) < 1) {
        stop(
            'file ', file_name, ' is malformed: it must contain at least one ',
            '[Header] line and at least one [Data] line'
        )
    }

    if (length(header_indices) != length(data_indices)) {
        stop(
            'file ', file_name, ' is malformed: it contains different numbers ',
            'of [Header] and [Data] lines'
        )
    }

    # Define a helping function for reading a row as a data frame and then
    # replacing any Unicode characters
    get_processed_row_data_frame <- function(skip) {
        row_data      <- licor_6800_lines_to_df(file_lines, skip)
        row_data[1, ] <- replace_unicode(row_data[1, ])
        row_data
    }

    # Get all the data chunks
    data_chunks <- lapply(seq_along(data_indices), function(i) {
        # Get the index of the beginning of this data chunk
        data_indx <- data_indices[i]

        # Get the number of data points in this chunk
        nlines <- if (i < length(data_indices)) {
            header_indices[i + 1] - data_indx - 4
        } else {
            length(file_lines) - data_indx - 3
        }

        # Read the column names, categories, units, and values from this data
        # chunk
        licor_variable_names      <- get_processed_row_data_frame(data_indx + 2)
        licor_variable_units      <- get_processed_row_data_frame(data_indx + 3)
        licor_variable_categories <- get_processed_row_data_frame(data_indx + 1)

        licor_data <- licor_6800_lines_to_df(
          file_lines,
          seq(data_indx + 4, length.out = nlines)
        )

        # Convert the data to numeric values whenever possible
        licor_data <- as.data.frame(
            lapply(licor_data, try_as_numeric),
            stringsAsFactors = FALSE
        )

        # Remove NA rows if necessary
        if (remove_NA_rows) {
            all_NA <- sapply(seq_len(nrow(licor_data)), function(i) {
                all(is.na(as.list(licor_data[i, ])))
            })
            licor_data <- licor_data[!all_NA, ]
        }

        # Get the column names as a vector
        licor_variable_names <-
            make.unique(as.character(licor_variable_names[1, ]))

        # Apply column names
        colnames(licor_variable_units) <- licor_variable_names
        colnames(licor_variable_categories) <- licor_variable_names
        colnames(licor_data) <- licor_variable_names

        # Return this chunk as an exdf object
        exdf(
            licor_data,
            licor_variable_units,
            licor_variable_categories
        )
    })

    # Get all the header chunks
    header_chunks <- lapply(seq_along(header_indices), function(i) {
        # Get the index of the beginning of this header chunk
        header_indx <- header_indices[i]

        # Get the number of points in this header
        nlines <- data_indices[i] - header_indx - 1

        # Read the header information
        preamble_raw <- licor_6800_lines_to_df(
          file_lines,
          seq(header_indx + 1, length.out = nlines)
        )

        # Get just the values (as numeric when possible), and add the names
        licor_preamble <- as.data.frame(
            t(preamble_raw[, seq(2, ncol(preamble_raw))]),
            stringsAsFactors = FALSE
        )

        licor_preamble <- as.data.frame(
            lapply(licor_preamble, try_as_numeric),
            stringsAsFactors = FALSE
        )

        colnames(licor_preamble) <- preamble_raw[, 1]

        # Replace any unicode
        colnames(licor_preamble) <- replace_unicode(colnames(licor_preamble))
        licor_preamble[1, ]      <- replace_unicode(licor_preamble[1, ])

        licor_preamble
    })

    # Get the names of all columns that are present in all of the chunks
    data_columns_to_keep   <- do.call(identify_common_columns, data_chunks)
    header_columns_to_keep <- do.call(identify_common_columns, header_chunks)

    # Extract just these columns
    data_chunks <- lapply(data_chunks, function(x) {
      x[ , data_columns_to_keep, TRUE]
    })

    header_chunks <- lapply(header_chunks, function(x) {
      x[ , header_columns_to_keep]
    })

    # Use `rbind` to combine all the chunks
    exdf_obj    <- do.call(rbind, data_chunks)
    header_part <- do.call(rbind, header_chunks)

    # Process the header
    processed_header <- process_6800_plaintext_header(header_part)
    preamble_exdf    <- processed_header$preamble_exdf
    remarks          <- processed_header$remarks

    # Incorporate preamble info into the main data table
    exdf_obj <- cbind(exdf_obj, preamble_exdf)

    # Store additional information in the data exdf
    exdf_obj$preamble     <- preamble_exdf$main_data
    exdf_obj$remarks      <- remarks
    exdf_obj$user_remarks <- user_remarks

    # Process the user constant rows
    exdf_obj <- process_user_constant_rows_6800_plaintext(exdf_obj)

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
