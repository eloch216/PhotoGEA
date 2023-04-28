identify_tdl_cycles <- function(
    tdl_exdf,
    valve_column_name,
    cycle_start_valve,
    expected_cycle_length_minutes,
    expected_cycle_num_valves,
    expected_cycle_num_time_pts = expected_cycle_num_valves,
    timestamp_colname
)
{
    if (!is.exdf(tdl_exdf)) {
        stop("identify_tdl_cycles requires an exdf object")
    }

    # Make sure the required variables are defined
    required_variables <- list()
    required_variables[[valve_column_name]] <- NA
    required_variables[[timestamp_colname]] <- NA

    check_required_variables(tdl_exdf, required_variables)

    # Add a new column to the data for the cycle number
    tdl_exdf <- document_variables(
        tdl_exdf,
        c("calculated", "cycle_num", "NA"),
        c("calculated", "elapsed_time", "minutes")
    )

    # Make sure the data isn't empty
    if (nrow(tdl_exdf) < 1) {
        stop("TDL cycles cannot be identified because the TDL data is empty")
    }

    # Make sure the data is sorted in chronological order
    tdl_exdf$main_data <- tdl_exdf[order(tdl_exdf[,timestamp_colname]),]

    # Find all the row numbers in the TDL data where the site value is equal to
    # the cycle start site value
    starting_valve_indices <- which(tdl_exdf[,valve_column_name] == cycle_start_valve)

    # If the valve number is the same for sequential rows, this represents
    # multiple measurements. When this happens, only the first instance of the
    # start valve in a sequential group is potentially the start of a cycle.
    # Here we find just these points, which are the potential starting indices.
    starting_indices <- 1
    for (i in seq(2, length(starting_valve_indices))) {
        if (starting_valve_indices[i] - starting_valve_indices[i-1] > 1) {
            starting_indices <- c(starting_indices, starting_valve_indices[i])
        }
    }

    # Check for problems
    if (length(starting_indices) < 1) {
        msg <- paste0(
            "No TDL cycles could be identified because the ",
            valve_column_name,
            " column is never equal to the cycle start value of ",
            cycle_start_valve
        )
        stop(msg)
    }

    # Make a helping function to check whether a possible cycle is valid
    check_cycle <- function(possible_cycle) {
        # Get info about the possible cycle
        n_time_pts <- nrow(possible_cycle)
        n_valves <- length(unique(possible_cycle[[valve_column_name]]))
        tdiff <- difftime(
            max(possible_cycle[[timestamp_colname]]),
            min(possible_cycle[[timestamp_colname]]),
            units = "min"
        )

        # Check for validity
        if (n_time_pts != expected_cycle_num_time_pts) {
            # We have the wrong number of time points
            return(FALSE)
        }
        else if (n_valves != expected_cycle_num_valves) {
            # We have the wrong number of valves
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
            ncol = ncol(tdl_exdf),
            nrow = 0
        ),
        stringsAsFactors = FALSE
    )
    colnames(new_main_data) <- colnames(tdl_exdf)

    # Get the first time point
    start_time <- tdl_exdf[1, timestamp_colname]

    # Fill in the new data frame with only valid TDL cycles
    n_cycles <- 1
    for (i in seq_along(starting_indices)) {
        possible_cycle <- c()
        if (i < length(starting_indices)) {
            possible_cycle <- tdl_exdf[seq(starting_indices[i], starting_indices[i+1] - 1),]
        } else {
            possible_cycle <- tdl_exdf[seq(starting_indices[i], nrow(tdl_exdf)),]
        }

        if (check_cycle(possible_cycle)) {
            possible_cycle[['cycle_num']] <- n_cycles

            new_main_data <- rbind(new_main_data, possible_cycle)
            n_cycles <- n_cycles + 1
        }
    }

    new_main_data[['elapsed_time']] <- as.double(difftime(
        new_main_data[[timestamp_colname]],
        start_time,
        units = "min"
    ))

    # Clean up and return the result
    row.names(new_main_data) <- NULL
    tdl_exdf$main_data <- new_main_data
    return(tdl_exdf)
}
