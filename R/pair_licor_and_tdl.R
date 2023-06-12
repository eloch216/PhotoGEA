pair_licor_and_tdl <- function(
    licor_file,
    tdl_data,
    licor_timestamp_column_name,
    tdl_timestamp_column_name,
    max_allowed_time_difference  # minutes
)
{
    # Add some new columns to the Licor file in preparation for adding the TDL
    # data
    licor_file <- document_variables(
        licor_file,
        c("calibrated TDL",              "cycle_num",              ""),
        c("calibrated TDL (sample)",     "tdl_time_s",             ""),
        c("calibrated TDL (sample)",     "calibrated_12c_s",       "ppm"),
        c("calibrated TDL (sample)",     "calibrated_13c_s",       "ppm"),
        c("calibrated TDL (sample)",     "total_mixing_ratio_s",   "ppm"),
        c("calibrated TDL (sample)",     "total_isotope_ratio_s",  "ppt"),
        c("calibrated TDL (reference)",  "tdl_time_r",             ""),
        c("calibrated TDL (reference)",  "calibrated_12c_r",       "ppm"),
        c("calibrated TDL (reference)",  "calibrated_13c_r",       "ppm"),
        c("calibrated TDL (reference)",  "total_mixing_ratio_r",   "ppm"),
        c("calibrated TDL (reference)",  "total_isotope_ratio_r",  "ppt")
    )

    # Make sure the time columns have the correct class, or we won't be able to
    # store timestamps in them properly
    licor_file[['main_data']][['tdl_time_s']] <-
        as.POSIXlt(licor_file[['main_data']][['tdl_time_s']])

    licor_file[['main_data']][['tdl_time_r']] <-
        as.POSIXlt(licor_file[['main_data']][['tdl_time_r']])

    # Get the sample and reference valve numbers
    sample_valve <- licor_file[1, 'valve_number_s']
    reference_valve <- licor_file[1, 'valve_number_r']

    # Now extract TDL data for each time point
    for (i in seq_len(nrow(licor_file[['main_data']]))) {
        # Find the TDL cycle that contains the closest time to the Licor data
        # point time
        licor_time <- licor_file[['main_data']][[licor_timestamp_column_name]][i]
        time_differences <- abs(difftime(
            tdl_data[[tdl_timestamp_column_name]],
            licor_time,
            units = 'min'
        ))
        min_time_difference <- min(time_differences)
        indx_of_closest_tdl_pnt <- match(min_time_difference, time_differences)
        cycle_of_closest_tdl_pnt <- tdl_data[['cycle_num']][indx_of_closest_tdl_pnt]

        if (min_time_difference > max_allowed_time_difference) {
            msg <- paste0(
                "Could not find a time point in the TDL data corresponding to ",
                format(licor_time),
                ", a point in the Licor data.\nThe nearest time in the TDL data was ",
                format(tdl_data[[tdl_timestamp_column_name]][indx_of_closest_tdl_pnt]),
                ",\nwhich is ", min_time_difference, " minutes away, exceeding ",
                "the maximum allowed value of ", max_allowed_time_difference,
                " minutes"
            )
            warning(msg)

            # There is no associated TDL data for this point
            tdl_time_sample <- NA
            tdl_time_reference <- NA
            calibrated_12c_s <- NA
            calibrated_13c_s <- NA
            total_mixing_ratio_s <- NA
            total_isotope_ratio_s <- NA
            calibrated_12c_r <- NA
            calibrated_13c_r <- NA
            total_mixing_ratio_r <- NA
            total_isotope_ratio_r <- NA
        } else {
            # Find the TDL times for the sample and reference measurements
            tdl_time_sample <-
                tdl_data[tdl_data[['cycle_num']] == cycle_of_closest_tdl_pnt &
                    tdl_data[['valve_number']] == sample_valve, tdl_timestamp_column_name]

            tdl_time_reference <-
                tdl_data[tdl_data[['cycle_num']] == cycle_of_closest_tdl_pnt &
                    tdl_data[['valve_number']] == reference_valve, tdl_timestamp_column_name]

            # Find the sample CO2 data
            calibrated_12c_s <-
                tdl_data[tdl_data[['cycle_num']] == cycle_of_closest_tdl_pnt &
                    tdl_data[['valve_number']] == sample_valve, 'calibrated_12c']

            calibrated_13c_s <-
                tdl_data[tdl_data[['cycle_num']] == cycle_of_closest_tdl_pnt &
                    tdl_data[['valve_number']] == sample_valve, 'calibrated_13c']

            total_mixing_ratio_s <-
                tdl_data[tdl_data[['cycle_num']] == cycle_of_closest_tdl_pnt &
                    tdl_data[['valve_number']] == sample_valve, 'total_mixing_ratio']

            total_isotope_ratio_s <-
                tdl_data[tdl_data[['cycle_num']] == cycle_of_closest_tdl_pnt &
                    tdl_data[['valve_number']] == sample_valve, 'total_isotope_ratio']

            # Find the reference CO2 data
            calibrated_12c_r <-
                tdl_data[tdl_data[['cycle_num']] == cycle_of_closest_tdl_pnt &
                    tdl_data[['valve_number']] == reference_valve, 'calibrated_12c']

            calibrated_13c_r <-
                tdl_data[tdl_data[['cycle_num']] == cycle_of_closest_tdl_pnt &
                    tdl_data[['valve_number']] == reference_valve, 'calibrated_13c']

            total_mixing_ratio_r <-
                tdl_data[tdl_data[['cycle_num']] == cycle_of_closest_tdl_pnt &
                    tdl_data[['valve_number']] == reference_valve, 'total_mixing_ratio']

            total_isotope_ratio_r <-
                tdl_data[tdl_data[['cycle_num']] == cycle_of_closest_tdl_pnt &
                    tdl_data[['valve_number']] == reference_valve, 'total_isotope_ratio']
        }

        # Store results in the Licor file
        licor_file[['main_data']][['cycle_num']][i] <- cycle_of_closest_tdl_pnt
        licor_file[['main_data']][['tdl_time_s']][i] <- tdl_time_sample
        licor_file[['main_data']][['tdl_time_r']][i] <- tdl_time_reference

        licor_file[['main_data']][['calibrated_12c_s']][i] <- calibrated_12c_s
        licor_file[['main_data']][['calibrated_13c_s']][i] <- calibrated_13c_s
        licor_file[['main_data']][['total_mixing_ratio_s']][i] <- total_mixing_ratio_s
        licor_file[['main_data']][['total_isotope_ratio_s']][i] <- total_isotope_ratio_s

        licor_file[['main_data']][['calibrated_12c_r']][i] <- calibrated_12c_r
        licor_file[['main_data']][['calibrated_13c_r']][i] <- calibrated_13c_r
        licor_file[['main_data']][['total_mixing_ratio_r']][i] <- total_mixing_ratio_r
        licor_file[['main_data']][['total_isotope_ratio_r']][i] <- total_isotope_ratio_r
    }

    return(licor_file)
}
