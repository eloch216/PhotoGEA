# identify_tdl_cycles: a function that identifies complete TDL measurement
# cycles, i.e., sets of one measurement from each "site" in the system.
#
# ------------------------------------------------------------------------------
#
# INPUTS:
#
# - tdl_exdf: an exdf object representing data from a TDL file
#
# - valve_column_name: the name of the data column that represents the site
#
# - cycle_start_valve: the value of the site column that indicates the
#       start of a new cycle
#
# - expected_cycle_length_minutes: the expected length of a full cycle in
#       minutes
#
# - expected_cycle_num_pts: the expected number of points in a full cycle
#
# - timestamp_colname: the name of the data column that represents the timestamp
#
# ------------------------------------------------------------------------------
#
# OUTPUT:
#
# an exdf object representing data from a TDL file, but now only including full
# cycles, each of which is numbered in a new `num_cycles` column.
#
identify_tdl_cycles <- function(
    tdl_exdf,
    valve_column_name,
    cycle_start_valve,
    expected_cycle_length_minutes,
    expected_cycle_num_pts,
    timestamp_colname
)
{
    # Add a new column to the data for the cycle number
    tdl_exdf <- specify_variables(
        tdl_exdf,
        c("calculated", "cycle_num", "NA")
    )

    # Make sure the data isn't empty
    if (nrow(tdl_exdf) < 1) {
        stop("TDL cycles cannot be identified because the TDL data is empty")
    }

    # Make sure the data is sorted in chronological order
    tdl_exdf[['main_data']] <- tdl_exdf[order(tdl_exdf[,timestamp_colname]),]

    # Find all the row numbers in the TDL data where the site value is equal to
    # the cycle start site value
    starting_indices <- which(tdl_exdf[,valve_column_name] == cycle_start_valve)

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
        n_pts <- nrow(possible_cycle)
        n_unique_pts <- length(unique(possible_cycle[[valve_column_name]]))
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
            ncol = ncol(tdl_exdf),
            nrow = 0
        ),
        stringsAsFactors = FALSE
    )
    colnames(new_main_data) <- colnames(tdl_exdf)

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

    # Clean up and return the result
    row.names(new_main_data) <- NULL

    return(
        exdf(
            new_main_data,
            tdl_exdf[['units']],
            tdl_exdf[['categories']]
        )
    )
}
