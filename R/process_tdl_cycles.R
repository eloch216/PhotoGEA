process_erml_tdl_cycle <- function(
    tdl_cycle,                        # data frame
    valve_column_name,                # column name string
    noaa_valve,                       # value of valve column
    calibration_0_valve,              # value of valve column
    calibration_1_valve,              # value of valve column
    calibration_2_valve,              # value of valve column
    calibration_3_valve,              # value of valve column
    raw_12c_colname,                  # column name string
    raw_13c_colname,                  # column name string
    noaa_cylinder_co2_concentration,  # ppm
    noaa_cylinder_isotope_ratio,      # ppt
    calibration_isotope_ratio,        # ppt
    f_other,                          # fraction of CO2 that is not 13C16O16O or 12C16O16O
    R_VPDB
)
{
    # Make a correction by subtracting the carbon concentrations measured in a
    # line that should have no carbon
    zero_carbon_reference <-
        tdl_cycle[tdl_cycle[[valve_column_name]] == calibration_0_valve,]

    tdl_cycle[['zero_corrected_12c']] <-
        tdl_cycle[[raw_12c_colname]] - zero_carbon_reference[[raw_12c_colname]]

    tdl_cycle[['zero_corrected_13c']] <-
        tdl_cycle[[raw_13c_colname]] - zero_carbon_reference[[raw_13c_colname]]

    # Make a data frame to store the zero information
    calibration_zero <- data.frame(
        cycle_num = tdl_cycle[['cycle_num']][1],
        offset_12c = zero_carbon_reference[[raw_12c_colname]],
        offset_13c = zero_carbon_reference[[raw_13c_colname]]
    )

    # Make adjustments for 12C calibration using a standard cylinder from NOAA.
    # Here we determine a multiplicative factor for converting the measured
    # zero-corrected 12C concentrations to true concentrations.
    total_mixing_ratio_noaa <- noaa_cylinder_co2_concentration * (1 - f_other)
    R_noaa <- R_VPDB * (1 + noaa_cylinder_isotope_ratio / 1000)
    noaa_12C16O16O <- total_mixing_ratio_noaa / (1 + R_noaa)

    noaa_reference <-
        tdl_cycle[tdl_cycle[[valve_column_name]] == noaa_valve,]

    gain_12CO2 <- noaa_12C16O16O / noaa_reference[['zero_corrected_12c']]

    tdl_cycle[['calibrated_12c']] <-
        tdl_cycle[['zero_corrected_12c']] * gain_12CO2

    # Make a data frame containing the info used to calibrate the 12CO2 data
    calibration_12CO2 <- data.frame(
        cycle_num = tdl_cycle[['cycle_num']][1],
        total_mixing_ratio_noaa = total_mixing_ratio_noaa,
        R_noaa = R_noaa,
        noaa_12C16O16O = noaa_12C16O16O,
        gain_12CO2 = gain_12CO2
    )

    # Make adjustments for 13C calibration using another reference that has been
    # characterized using the Long lab's mass spectrometer. This tank is mixed
    # with varying amounts of N2 and measured as three separate lines. Here we
    # determine a polynomial function for converting the measured zero-corrected
    # 13C concentrations to true concentrations. The expected values are
    # determined from the calibrated 12C values using the reference cylinder's
    # known isotope ratio.
    conversion_factor <- R_VPDB * (calibration_isotope_ratio / 1000 + 1)
    expected_13c_values <- c(
        tdl_cycle[tdl_cycle[[valve_column_name]] == calibration_0_valve, 'calibrated_12c'] * conversion_factor,
        tdl_cycle[tdl_cycle[[valve_column_name]] == calibration_1_valve, 'calibrated_12c'] * conversion_factor,
        tdl_cycle[tdl_cycle[[valve_column_name]] == calibration_2_valve, 'calibrated_12c'] * conversion_factor,
        tdl_cycle[tdl_cycle[[valve_column_name]] == calibration_3_valve, 'calibrated_12c'] * conversion_factor
    )

    measured_13c_values <- c(
        tdl_cycle[tdl_cycle[[valve_column_name]] == calibration_0_valve, 'zero_corrected_13c'],
        tdl_cycle[tdl_cycle[[valve_column_name]] == calibration_1_valve, 'zero_corrected_13c'],
        tdl_cycle[tdl_cycle[[valve_column_name]] == calibration_2_valve, 'zero_corrected_13c'],
        tdl_cycle[tdl_cycle[[valve_column_name]] == calibration_3_valve, 'zero_corrected_13c']
    )

    # Fit a quadratic model
    quadratic_fit <-
        lm(expected_13c_values ~ poly(measured_13c_values, 2, raw = TRUE))

    # Get the coefficients
    fit_summary <- summary(quadratic_fit)
    fit_coeff <- fit_summary[['coefficients']]

    a0 <- fit_coeff[[1]]  # intercept
    a1 <- fit_coeff[[2]]  # coefficient of x
    a2 <- fit_coeff[[3]]  # coefficient of x^2

    r_squared <- fit_summary[['r.squared']]

    # Make a data frame containing the conversion factor and fit parameters
    calibration_13CO2_fit <- data.frame(
        cycle_num = tdl_cycle[['cycle_num']][1],
        calibrated_12c_to_13c_conversion_factor = conversion_factor,
        a0 = a0,
        a1 = a1,
        a2 = a2,
        r_squared = r_squared
    )

    # Determine calibrated 13C values
    fitted_13c_values <-
        a0 + a1 * measured_13c_values + a2 * measured_13c_values^2

    tdl_cycle[['calibrated_13c']] <-
        a0 + a1 * tdl_cycle[['zero_corrected_13c']] + a2 * tdl_cycle[['zero_corrected_13c']]^2

    # Make a data frame containing the points used for the fit
    calibration_13CO2_data <- data.frame(
        cycle_num = tdl_cycle[['cycle_num']][1],
        measured_13c_values = measured_13c_values,
        expected_13c_values = expected_13c_values,
        fitted_13c_values = fitted_13c_values
    )

    # Determine the total mixing and isotope ratios from the calibrated 12C and
    # 13C concentrations using Equations A.1 and A.2 from Griffis et al.
    # Agricultural and Forest Meteorology 124, 15-29 (2004)
    # (https://doi.org/10.1016/j.agrformet.2004.01.009).
    tdl_cycle[['total_mixing_ratio']] <-
        tdl_cycle[['calibrated_13c']] + tdl_cycle[['calibrated_12c']]

    tdl_cycle[['total_isotope_ratio']] <-
        1000 * (tdl_cycle[['calibrated_13c']] / tdl_cycle[['calibrated_12c']] / R_VPDB - 1)

    return(list(
        tdl_data = tdl_cycle,
        calibration_zero = calibration_zero,
        calibration_12CO2 = calibration_12CO2,
        calibration_13CO2_data = calibration_13CO2_data,
        calibration_13CO2_fit = calibration_13CO2_fit
    ))
}

process_tdl_cycles <- function(
    tdl_data,
    cycle_processing_function,
    ...
)
{
    # Make sure the cycle column is defined
    if (!'cycle_num' %in% colnames(tdl_data)) {
        stop("The TDL must have a column called 'cycle_num', which is typically added via the 'identify_tdl_cycles' function")
    }

    # Make sure there is at least one TDL cycle in the data
    cycles <- unique(tdl_data[['cycle_num']])
    if (length(cycles) < 1) {
        stop("The TDL data must contain at least one cycle")
    }

    # Process each cycle, creating a list where each element is a list created
    # by `cycle_processing_function`
    cycle_result_list <- by(
        tdl_data,
        tdl_data[['cycle_num']],
        function(x) {cycle_processing_function(x, ...)}
    )

    # Get the names of the elements of the return value from `process_tdl_cycle`
    element_names <- names(cycle_result_list[[1]])

    # Convert the result list into a list of data frames and return it
    return(
        setNames(
            lapply(
                element_names,
                function(n) {
                    smart_rbind(lapply(cycle_result_list, function(x) {x[[n]]}))
                }
            ),
            element_names
        )
    )
}
