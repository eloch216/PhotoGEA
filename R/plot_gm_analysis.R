source("gm_analysis.R")

MAKE_TDL_PLOTS <- FALSE

MAKE_GM_PLOTS <- TRUE

# Make plots, if desired
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

if (MAKE_GM_PLOTS) {
    # Make gmc boxplots for each event
    gmc_event_boxplot <- bwplot(
        gmc ~ event | genotype,
        data = licor_files_no_outliers[['main_data']],
        ylim = c(0, MAX_GM),
        ylab = "Mesophyll conductance to CO2 (mol / m^2 / s / bar)",
    )
    x11()
    print(gmc_event_boxplot)

    # Make gmc boxplots for each rep
    gmc_rep_boxplot <- bwplot(
        gmc ~ event_rep | genotype,
        data = licor_files_no_outliers[['main_data']],
        ylim = c(0, MAX_GM),
        ylab = "Mesophyll conductance to CO2 (mol / m^2 / s / bar)"
    )
    x11()
    print(gmc_rep_boxplot)

    # Make Cc boxplots for each event
    cc_event_boxplot <- bwplot(
        Cc ~ event | genotype,
        data = licor_files_no_outliers[['main_data']],
        #ylim = c(0, MAX_GM),
        ylab = "CO2 concentration in chloroplast (micromol / mol)"
    )
    x11()
    print(cc_event_boxplot)

    # Make Cc boxplots for each rep
    cc_rep_boxplot <- bwplot(
        Cc ~ event_rep | genotype,
        data = licor_files_no_outliers[['main_data']],
        #ylim = c(0, MAX_GM),
        ylab = "CO2 concentration in chloroplast (micromol / mol)"
    )
    x11()
    print(cc_rep_boxplot)

    # Make Cc boxplots for each event
    cc_event_boxplot <- bwplot(
        Cc ~ event | genotype,
        data = licor_files_no_outliers[['main_data']],
        ylim = c(0, 275),
        ylab = "CO2 concentration in chloroplast (micromol / mol)"
    )
    x11()
    print(cc_event_boxplot)

    # Make Cc boxplots for each rep
    cc_rep_boxplot <- bwplot(
        Cc ~ event_rep | genotype,
        data = licor_files_no_outliers[['main_data']],
        ylim = c(0, 275),
        ylab = "CO2 concentration in chloroplast (micromol / mol)"
    )
    x11()
    print(cc_rep_boxplot)

    # Make drawdown boxplots for each event
    drawdown_event_boxplot <- bwplot(
        Ci-Cc ~ event | genotype,
        data = licor_files_no_outliers[['main_data']],
        ylim = c(0, 150),
        ylab = "CO2 drawdown to chloroplast (Ci - Cc) (micromol / mol)"
    )
    x11()
    print(drawdown_event_boxplot)

    # Make drawdown boxplots for each rep
    drawdown_boxplot <- bwplot(
        Ci-Cc ~ event_rep | genotype,
        data = licor_files_no_outliers[['main_data']],
        ylim = c(0, 150),
        ylab = "CO2 drawdown to chloroplast (Ci - Cc) (micromol / mol)"
    )
    x11()
    print(drawdown_boxplot)
}
