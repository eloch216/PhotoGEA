smooth_tdl_data <- function(
    tdl_exdf,
    valve_column_name,
    valve_number,
    smoothing_function # should be a function with input arguments Y, X, that returns a smoothed version of Y
)
{
    # Extract data corresponding to a single valve
    valve_time_series <- tdl_exdf[tdl_exdf[, valve_column_name] == valve_number, ]

    # Apply the smoothing function to the valve time series
    filtered_12C <- smoothing_function(valve_time_series[['Conc12C_Avg']], valve_time_series[['elapsed_time']])
    filtered_13C <- smoothing_function(valve_time_series[['Conc13C_Avg']], valve_time_series[['elapsed_time']])

    # Store it in the exdf object
    tdl_data <- tdl_exdf[['main_data']]
    tdl_data[tdl_data[[valve_column_name]] == valve_number, 'Conc12C_Avg'] <- filtered_12C
    tdl_data[tdl_data[[valve_column_name]] == valve_number, 'Conc13C_Avg'] <- filtered_13C
    tdl_exdf[['main_data']] <- tdl_data

    # Return the modified exdf object
    return(tdl_exdf)
}
