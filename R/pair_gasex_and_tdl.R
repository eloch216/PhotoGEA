pair_gasex_and_tdl <- function(
    gasex_exdf,
    tdl_exdf,
    max_allowed_time_difference = 1,  # minutes
    gasex_timestamp_column_name = 'time',
    tdl_timestamp_column_name = 'TIMESTAMP'
)
{
    # Check the gas exchange data
    if (!is.exdf(gasex_exdf)) {
        stop('gasex_exdf must be an exdf object')
    }

    # Make sure the required variables are defined and have the correct units
    gasex_required_variables <- list()
    gasex_required_variables[['valve_number_s']] <- NA
    gasex_required_variables[['valve_number_r']] <- NA
    gasex_required_variables[[gasex_timestamp_column_name]] <- NA

    check_required_variables(gasex_exdf, gasex_required_variables)

    # Check the TDL data
    if (!is.exdf(tdl_exdf)) {
        stop('tdl_exdf must be an exdf object')
    }

    # Make sure the required variables are defined and have the correct units
    tdl_required_variables <- list()
    tdl_required_variables[['cycle_num']] <- NA
    tdl_required_variables[['valve_number']] <- NA
    tdl_required_variables[['calibrated_12c']] <- 'ppm'
    tdl_required_variables[['calibrated_13c']] <- 'ppm'
    tdl_required_variables[['total_CO2']] <- 'ppm'
    tdl_required_variables[['delta_C13']] <- 'ppt'
    tdl_required_variables[[tdl_timestamp_column_name]] <- NA

    check_required_variables(tdl_exdf, tdl_required_variables)

    # Add some new columns to the Licor file in preparation for adding the TDL
    # data
    gasex_exdf <- document_variables(
        gasex_exdf,
        c('pair_gasex_and_tdl', 'cycle_num',        ''),
        c('pair_gasex_and_tdl', 'tdl_time_s',       ''),
        c('pair_gasex_and_tdl', 'calibrated_12c_s', 'ppm'),
        c('pair_gasex_and_tdl', 'calibrated_13c_s', 'ppm'),
        c('pair_gasex_and_tdl', 'total_CO2_s',      'ppm'),
        c('pair_gasex_and_tdl', 'delta_C13_s',      'ppt'),
        c('pair_gasex_and_tdl', 'tdl_time_r',       ''),
        c('pair_gasex_and_tdl', 'calibrated_12c_r', 'ppm'),
        c('pair_gasex_and_tdl', 'calibrated_13c_r', 'ppm'),
        c('pair_gasex_and_tdl', 'total_CO2_r',      'ppm'),
        c('pair_gasex_and_tdl', 'delta_C13_r',      'ppt')
    )

    # Make sure the time columns have the correct class, or we won't be able to
    # store timestamps in them properly
    gasex_exdf[['main_data']][['tdl_time_s']] <-
        as.POSIXlt(gasex_exdf[['main_data']][['tdl_time_s']])

    gasex_exdf[['main_data']][['tdl_time_r']] <-
        as.POSIXlt(gasex_exdf[['main_data']][['tdl_time_r']])

    # Now extract TDL data for each time point
    for (i in seq_len(nrow(gasex_exdf[['main_data']]))) {
        # Get the sample and reference valve numbers
        sample_valve <- gasex_exdf[i, 'valve_number_s']
        reference_valve <- gasex_exdf[i, 'valve_number_r']

        # Find the TDL cycle that contains the closest time to the Licor data
        # point time
        licor_time <- gasex_exdf[i, gasex_timestamp_column_name]
        time_differences <- abs(difftime(
            tdl_exdf[, tdl_timestamp_column_name],
            licor_time,
            units = 'min'
        ))
        min_time_difference <- min(time_differences)
        indx_of_closest_tdl_pnt <- match(min_time_difference, time_differences)
        cycle_of_closest_tdl_pnt <- tdl_exdf[indx_of_closest_tdl_pnt, 'cycle_num']

        if (min_time_difference > max_allowed_time_difference) {
            msg <- paste0(
                'Could not find a time point in the TDL data corresponding to ',
                format(licor_time),
                ', a point in the Licor data.\nThe nearest time in the TDL data was ',
                format(tdl_exdf[indx_of_closest_tdl_pnt, tdl_timestamp_column_name]),
                ',\nwhich is ', min_time_difference, ' minutes away, exceeding ',
                'the maximum allowed value of ', max_allowed_time_difference,
                ' minutes'
            )
            warning(msg)

            # There is no associated TDL data for this point
            tdl_time_sample <- NA
            tdl_time_reference <- NA
            calibrated_12c_s <- NA
            calibrated_13c_s <- NA
            total_CO2_s <- NA
            delta_C13_s <- NA
            calibrated_12c_r <- NA
            calibrated_13c_r <- NA
            total_CO2_r <- NA
            delta_C13_r <- NA
        } else {
            # Get the closest TDL cycle (as a data frame)
            closest_cycle <-
                tdl_exdf[tdl_exdf[, 'cycle_num'] == cycle_of_closest_tdl_pnt, ]

            # Get the sample and reference valve components of the closest cycle
            closest_cycle_sample <-
                closest_cycle[closest_cycle[, 'valve_number'] == sample_valve, ]

            closest_cycle_reference <-
                closest_cycle[closest_cycle[, 'valve_number'] == reference_valve, ]

            # Find the TDL times for the sample and reference measurements
            tdl_time_sample <- closest_cycle_sample[, tdl_timestamp_column_name]
            tdl_time_reference <- closest_cycle_reference[, tdl_timestamp_column_name]

            # Find the sample CO2 data
            calibrated_12c_s <- closest_cycle_sample[, 'calibrated_12c']
            calibrated_13c_s <- closest_cycle_sample[, 'calibrated_13c']
            total_CO2_s <- closest_cycle_sample[, 'total_CO2']
            delta_C13_s <- closest_cycle_sample[, 'delta_C13']

            # Find the reference CO2 data
            calibrated_12c_r <- closest_cycle_reference[, 'calibrated_12c']
            calibrated_13c_r <- closest_cycle_reference[, 'calibrated_13c']
            total_CO2_r <- closest_cycle_reference[, 'total_CO2']
            delta_C13_r <- closest_cycle_reference[, 'delta_C13']
        }

        # Store results in the Licor file
        gasex_exdf[i, 'cycle_num'] <- cycle_of_closest_tdl_pnt
        gasex_exdf[['main_data']][['tdl_time_s']][i] <- tdl_time_sample
        gasex_exdf[['main_data']][['tdl_time_r']][i] <- tdl_time_reference

        gasex_exdf[i, 'calibrated_12c_s'] <- calibrated_12c_s
        gasex_exdf[i, 'calibrated_13c_s'] <- calibrated_13c_s
        gasex_exdf[i, 'total_CO2_s'] <- total_CO2_s
        gasex_exdf[i, 'delta_C13_s'] <- delta_C13_s

        gasex_exdf[i, 'calibrated_12c_r'] <- calibrated_12c_r
        gasex_exdf[i, 'calibrated_13c_r'] <- calibrated_13c_r
        gasex_exdf[i, 'total_CO2_r'] <- total_CO2_r
        gasex_exdf[i, 'delta_C13_r'] <- delta_C13_r
    }

    return(gasex_exdf)
}
