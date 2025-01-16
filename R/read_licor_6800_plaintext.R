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

    licor_file[data_times > time_i, 'user_remark'] <- remark_i
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

# Main function for reading plaintext LI-6800 log files
read_licor_6800_plaintext <- function(
    file_name,
    get_oxygen = TRUE,
    include_user_remark_column = TRUE,
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

        # Apply column names
        colnames(licor_variable_units) <- licor_variable_names[1, ]
        colnames(licor_variable_categories) <- licor_variable_names[1, ]
        colnames(licor_data) <- licor_variable_names[1, ]

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

        licor_preamble <- stats::setNames(
            as.data.frame(
                t(preamble_raw[, seq(2, ncol(preamble_raw))]),
                stringsAsFactors = FALSE
            ),
            preamble_raw[, 1]
        )

        colnames(licor_preamble) <- replace_unicode(colnames(licor_preamble))
        licor_preamble[1, ] <- replace_unicode(licor_preamble[1, ])

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

    # Store additional information in the data exdf
    exdf_obj$preamble <- header_part
    exdf_obj$user_remarks <- user_remarks

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
