source("licor_data_operations.R")

process_tdl_cycle <- function(
    tdl_cycle,
    valve_column_name = 'valve_number',
    zero_carbon_reference_valve = 20,
    noaa_reference_valve = 2,
    ll_reference_valve_1 = 21,
    ll_reference_valve_2 = 23,
    ll_reference_valve_3 = 26,
    raw_12c_colname = 'Conc12C_Avg',
    raw_13c_colname = 'Conc13C_Avg',
    noaa_reference_co2_concentration = 294.996,     # ppm
    noaa_reference_isotope_ratio = -8.40,           # ppt
    ll_reference_isotope_ratio = -11.505,           # ppt
    f_other = 0.00474,                              # fraction of CO2 that is not 13C16O16O or 12C16O16O
    R_VPDB = 0.0111797
)
{
    # Make a correction by subtracting the carbon concentrations measured in a
    # line that should have no carbon
    zero_carbon_reference <-
        tdl_cycle[which(tdl_cycle[[valve_column_name]] == zero_carbon_reference_valve),]

    tdl_cycle[['zero_corrected_12c']] <-
        tdl_cycle[[raw_12c_colname]] - zero_carbon_reference[[raw_12c_colname]]

    tdl_cycle[['zero_corrected_13c']] <-
        tdl_cycle[[raw_13c_colname]] - zero_carbon_reference[[raw_13c_colname]]

    # Make adjustments for 12C calibration using a standard cylinder from NOAA.
    # Here we determine a multiplicative factor for converting the measured
    # zero-corrected 12C concentrations to true concentrations.
    total_mixing_ratio_noaa <- noaa_reference_co2_concentration * (1 - f_other)
    R_noaa <- R_VPDB * (1 + noaa_reference_isotope_ratio / 1000)
    noaa_12C16O16O <- total_mixing_ratio_noaa / (1 + R_noaa)

    noaa_reference <-
        tdl_cycle[which(tdl_cycle[[valve_column_name]] == noaa_reference_valve),]

    gain_12CO2 <- noaa_12C16O16O / noaa_reference[['zero_corrected_12c']]

    tdl_cycle[['calibrated_12c']] <-
        tdl_cycle[['zero_corrected_12c']] * gain_12CO2

    # Make adjustments for 13C calibration using another reference that has been
    # characterized using the Long lab's mass spectrometer. This tank is mixed
    # with varying amounts of N2 and measured as three separate lines. Here we
    # determine a polynomial function for converting the measured zero-corrected
    # 13C concentrations to true concentrations. The expected values are
    # determined from the calibrated 12C values using the reference cylinder's
    # known isotope ratio.
    conversion_factor <- R_VPDB * (ll_reference_isotope_ratio / 1000 + 1)
    expected_13c_values <- c(
        tdl_cycle[which(tdl_cycle[[valve_column_name]] == zero_carbon_reference_valve), 'calibrated_12c'] * conversion_factor,
        tdl_cycle[which(tdl_cycle[[valve_column_name]] == ll_reference_valve_1), 'calibrated_12c'] * conversion_factor,
        tdl_cycle[which(tdl_cycle[[valve_column_name]] == ll_reference_valve_2), 'calibrated_12c'] * conversion_factor,
        tdl_cycle[which(tdl_cycle[[valve_column_name]] == ll_reference_valve_3), 'calibrated_12c'] * conversion_factor
    )

    measured_13c_values <- c(
        tdl_cycle[which(tdl_cycle[[valve_column_name]] == zero_carbon_reference_valve), 'zero_corrected_13c'],
        tdl_cycle[which(tdl_cycle[[valve_column_name]] == ll_reference_valve_1), 'zero_corrected_13c'],
        tdl_cycle[which(tdl_cycle[[valve_column_name]] == ll_reference_valve_2), 'zero_corrected_13c'],
        tdl_cycle[which(tdl_cycle[[valve_column_name]] == ll_reference_valve_3), 'zero_corrected_13c']
    )

    # Fit a quadratic model
    quadratic_fit <- lm(expected_13c_values ~ poly(measured_13c_values, 2, raw=TRUE))

    # Get the coefficients
    fit_summary <- summary(quadratic_fit)
    fit_coeff <- fit_summary[['coefficients']]

    a0 <- fit_coeff[[1]]  # intercept
    a1 <- fit_coeff[[2]]  # coefficient of x
    a2 <- fit_coeff[[3]]  # coefficient of x^2

    r_squared <- fit_summary[['r.squared']]

    # Determine calibrated 13C values
    tdl_cycle[['calibrated_13c']] <-
        a0 + a1 * tdl_cycle[['zero_corrected_13c']] + a2 * tdl_cycle[['zero_corrected_13c']]^2

    # Determine the total mixing and isotope ratios from the calibrated 12C and
    # 13C concentrations using Equations A.1 and A.2 from Griffis et al.
    # Agricultural and Forest Meteorology 124, 15-29 (2004)
    # (https://doi.org/10.1016/j.agrformet.2004.01.009).
    tdl_cycle[['total_mixing_ratio']] <-
        tdl_cycle[['calibrated_13c']] + tdl_cycle[['calibrated_12c']]

    tdl_cycle[['total_isotope_ratio']] <-
        1000 * (tdl_cycle[['calibrated_13c']] / tdl_cycle[['calibrated_12c']] / R_VPDB - 1)

    return(tdl_cycle)

    return(list(
        tdl_cycle = tdl_cycle,
        gain_12CO2 = gain_12CO2,
        calibration_13CO2 = list(
            a0 = a0,
            a1 = a1,
            a2 = a2,
            r_squared = r_squared,
            measured_13c_values = measured_13c_values,
            expected_13c_values = expected_13c_values
        )
    )))
}
