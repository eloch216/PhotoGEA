smooth_tdl_data <- function(
    tdl_exdf,
    column_to_be_smoothed,
    valve_column_name,
    valve_number,
    smoothing_function # should be a function with input arguments Y, X, that returns a smoothed version of Y
)
{
    if (!is.exdf(tdl_exdf)) {
        stop("smooth_tdl_data requires an exdf object")
    }

    # Make sure the required variables are defined
    required_variables <- list()
    required_variables[[column_to_be_smoothed]] <- NA
    required_variables[[valve_column_name]] <- NA
    required_variables[['elapsed_time']] <- NA

    check_required_variables(tdl_exdf, required_variables)

    # Extract data corresponding to a single valve
    valve_time_series <- tdl_exdf[tdl_exdf[, valve_column_name] == valve_number, ]

    # Apply the smoothing function to the valve time series
    smoothed <- smoothing_function(valve_time_series[[column_to_be_smoothed]], valve_time_series[['elapsed_time']])

    # Store it in the exdf object
    tdl_exdf$main_data[tdl_exdf[, valve_column_name] == valve_number, column_to_be_smoothed] <- smoothed

    # Return the modified exdf object
    return(tdl_exdf)
}
