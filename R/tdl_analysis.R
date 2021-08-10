source("read_tdl.R")
source("tdl_data_operations.R")
source("tdl_calculations.R")

MAKE_TDL_PLOTS <- TRUE

tdl_files <- batch_read_tdl_file(choose_input_tdl_files())

tdl_files <- batch_extract_licor_variables(
    tdl_files,
    c(
        'TIMESTAMP',
        'valve_number',
        'Conc12C_Avg',
        'Conc13C_Avg'
    )
)

tdl_files <- combine_tdl_files(tdl_files)

tdl_files <- identify_tdl_cycles(tdl_files)

processed_tdl_data <- process_tdl_cycles(tdl_files[['main_data']])

if (MAKE_TDL_PLOTS) {
    # Make a plot of all the fits from the processing
    tdl_fitting <- xyplot(
        expected_13c_values + fitted_13c_values ~ measured_13c_values | factor(cycle_num),
        data = processed_tdl_data[['calibration_13CO2_data']],
        type = 'b',
        pch = 20,
        xlab = "Measured 13CO2 mixing ratio (ppm)",
        ylab = "True 13CO2 mixing ratio (ppm)",
        auto = TRUE
    )
    x11(width = 12, height = 6)
    print(tdl_fitting)

    # Make a plot showing how the zeroes drift with time
    tdl_12CO2_zero_drift <- xyplot(
        offset_12c ~ cycle_num,
        data = processed_tdl_data[['calibration_zero']],
        type = 'l',
        xlab = "TDL cycle number",
        ylab = "12CO2 zero offset"
    )
    x11()
    print(tdl_12CO2_zero_drift)

    tdl_13CO2_zero_drift <- xyplot(
        offset_13c ~ cycle_num,
        data = processed_tdl_data[['calibration_zero']],
        type = 'l',
        xlab = "TDL cycle number",
        ylab = "13CO2 zero offset"
    )
    x11()
    print(tdl_13CO2_zero_drift)

    # Make a plot showing how the 12CO2 calibration factor drifts with time
    tdl_12CO2_calibration_drift <- xyplot(
        gain_12CO2 ~ cycle_num,
        data = processed_tdl_data[['calibration_12CO2']],
        type = 'l',
        xlab = "TDL cycle number",
        ylab = "12CO2 gain factor"
    )
    x11()
    print(tdl_12CO2_calibration_drift)

    # Make plots showing how the 13CO2 calibration parameters drift with time
    tdl_13CO2_calibration_drift_a0 <- xyplot(
        a0 ~ cycle_num,
        data = processed_tdl_data[['calibration_13CO2_fit']],
        type = 'l',
        xlab = "TDL cycle number",
        ylab = "13CO2 polynomial calibration parameter a0",
        auto = TRUE,
        main = "true_13CO2 =\na0 +\na1 * measured_13CO2 +\na2 * measured_13CO2^2"
    )
    x11()
    print(tdl_13CO2_calibration_drift_a0)

    tdl_13CO2_calibration_drift_a1 <- xyplot(
        a1 ~ cycle_num,
        data = processed_tdl_data[['calibration_13CO2_fit']],
        type = 'l',
        xlab = "TDL cycle number",
        ylab = "13CO2 polynomial calibration parameter a1",
        auto = TRUE,
        main = "true_13CO2 =\na0 +\na1 * measured_13CO2 +\na2 * measured_13CO2^2"
    )
    x11()
    print(tdl_13CO2_calibration_drift_a1)


    tdl_13CO2_calibration_drift_a2 <- xyplot(
        a2 ~ cycle_num,
        data = processed_tdl_data[['calibration_13CO2_fit']],
        type = 'l',
        xlab = "TDL cycle number",
        ylab = "13CO2 polynomial calibration parameter a2",
        auto = TRUE,
        main = "true_13CO2 =\na0 +\na1 * measured_13CO2 +\na2 * measured_13CO2^2"
    )
    x11()
    print(tdl_13CO2_calibration_drift_a2)
}
