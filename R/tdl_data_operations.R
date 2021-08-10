# combine_tdl_files: a function for combining the information from multiple
# TDL files into a single list. Here, only the filenames, types, units, and data
# are retained (i.e., any parameters specified when reading the files will be
# lost).
#
# ------------------------------------------------------------------------------
#
# INPUTS:
#
# - tdl_files: a list of unnamed lists, each representing the information from
#       a single TDL file (typically produced by a call to the
#       `batch_read_tdl_file` function defined in `read_tdl.R`)
#
# ------------------------------------------------------------------------------
#
# OUTPUT:
#
# a list describing TDL data having the same structure as the output of a call
# to the `read_tdl_file` function.
#
combine_tdl_files <- function(
    tdl_files
)
{
    # Make sure at least one file is defined
    if (length(tdl_files) < 1) {
        stop("The input list is required to contain at least one TDL file")
    }

    # Get the info from the first file to use as a reference
    first_file <- tdl_files[[1]]

    # Set up the combo list
    initial_data_frame <- data.frame(
        matrix(
            ncol = ncol(first_file[['main_data']]),
            nrow = 0
        ),
        stringsAsFactors = FALSE
    )

    colnames(initial_data_frame) <- colnames(first_file[['main_data']])

    combo_info <- list(
        file_name = character(0),
        type = first_file[['type']],
        units = first_file[['units']],
        main_data = initial_data_frame,
        timestamp_colname = first_file[['timestamp_colname']]
    )

    # Get the remaining info
    for (i in seq_along(tdl_files)) {
        # Get the current file
        current_file <- tdl_files[[i]]

        # Check for possible errors
        if (!identical(
            colnames(first_file[['main_data']]),
            colnames(current_file[['main_data']])
        ))
        {
            msg <- paste0(
                "The column names specified in TDL file '",
                current_file[['file_name']],
                "' do not agree with the column names specified in '",
                first_file[['file_name']],
                "', so the two files cannot be combined"
            )
            stop(msg)
        }

        if (!identical(first_file[['types']], current_file[['types']])) {
            msg <- paste0(
                "The types specified in TDL file '",
                current_file[['file_name']],
                "' do not agree with the types specified in '",
                first_file[['file_name']],
                "', so the two files cannot be combined"
            )
            stop(msg)
        }

        if (!identical(first_file[['units']], current_file[['units']])) {
            msg <- paste0(
                "The units specified in TDL file '",
                current_file[['file_name']],
                "' do not agree with the units specified in '",
                first_file[['file_name']],
                "', so the two files cannot be combined"
            )
            stop(msg)
        }

        # Add this file to the combined info
        combo_info[['main_data']] <- rbind(
            combo_info[['main_data']],
            current_file[['main_data']]
        )

        combo_info[['file_name']] <- c(
            combo_info[['file_name']],
            current_file[['file_name']]
        )
    }

    return(combo_info)
}

# identify_tdl_cycles: a function that identifies complete TDL measurement
# cycles, i.e., sets of one measurement from each "site" in the system.
#
# ------------------------------------------------------------------------------
#
# INPUTS:
#
# - tdl_file: a list representing the data from a TDL file (typically produced
#       by a a call to the `read_tdl_file` defined in `read_licor.R`)
#
# - site_column_name: the name of the data column that represents the site
#
# - cycle_start_site_value: the value of the site column that indicates the
#       start of a new cycle
#
# - expected_cycle_length_minutes: the expected length of a full cycle in
#       minutes
#
# - expected_cycle_num_pts: the expected number of points in a full cycle
#
# ------------------------------------------------------------------------------
#
# OUTPUT:
#
# a list describing TDL data having the same structure as the output of a call
# to the `read_tdl_file` function, but now only including full cycles, each of
# which is numbered in a new `cycle` column.
#
identify_tdl_cycles <- function(
    tdl_file,
    site_column_name = 'valve_number',
    cycle_start_site_value = 20,
    expected_cycle_length_minutes = 2.7,
    expected_cycle_num_pts = 9
)
{
    main_data <- tdl_file[['main_data']]
    timestamp_colname <- tdl_file[['timestamp_colname']]

    # Make sure the data isn't empty
    if (nrow(main_data) < 1) {
        stop("TDL cycles cannot be identified because the TDL data is empty")
    }

    # Make sure the data is sorted in chronological order
    main_data <- main_data[order(main_data[[timestamp_colname]]),]

    # Find all the row numbers in the TDL data where the site value is equal to
    # the cycle start site value
    starting_indices <-
        which(main_data[[site_column_name]] == cycle_start_site_value)

    # Check for problems
    if (length(starting_indices) < 1) {
        msg <- paste0(
            "No TDL cycles could be identified because the ",
            site_column_name,
            " column is never equal to the cycle start value of ",
            cycle_start_site_value
        )
        stop(msg)
    }

    # Make a helping function to check whether a possible cycle is valid
    check_cycle <- function(possible_cycle) {
        # Get info about the possible cycle
        n_pts <- nrow(possible_cycle)
        n_unique_pts <- length(unique(possible_cycle[[site_column_name]]))
        tdiff <- difftime(
            max(possible_cycle[[timestamp_colname]]),
            min(possible_cycle[[timestamp_colname]]),
            units = "min"
        )

        # Check for validity
        if (n_pts != n_unique_pts) {
            # The site values should be unique within a cycle
            return(FALSE)
        } else if (n_pts != expected_cycle_num_pts) {
            return(FALSE)
        } else if (abs(tdiff - expected_cycle_length_minutes) > 0.5) {
            # Allow a tolerance away from the expected length
            return(FALSE)
        } else {
            return(TRUE)
        }
    }

    # Set up a new data frame
    new_main_data <- data.frame(
        matrix(
            ncol = ncol(main_data) + 1,
            nrow = 0
        ),
        stringsAsFactors = FALSE
    )
    colnames(new_main_data) <- c(colnames(main_data), 'cycle_num')

    # Fill in the new data frame with only valid TDL cycles
    n_cycles <- 1
    for (i in seq_along(starting_indices)) {
        possible_cycle <- c()
        if (i < length(starting_indices)) {
            possible_cycle <- main_data[seq(starting_indices[i], starting_indices[i+1] - 1),]
        } else {
            possible_cycle <- main_data[seq(starting_indices[i], nrow(main_data)),]
        }

        if (check_cycle(possible_cycle)) {
            possible_cycle[['cycle']] <- n_cycles
            new_main_data <- rbind(new_main_data, possible_cycle)
            n_cycles <- n_cycles + 1
        }
    }

    # Clean up and return the result
    row.names(new_main_data) <- NULL
    tdl_file[['main_data']] <- new_main_data
    return(tdl_file)
}
