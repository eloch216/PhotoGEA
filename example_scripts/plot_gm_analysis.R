# This script generates several plots from the calculations performed in
# `gm_analysis.R`. In fact, this script calls that one before making any plots.
#
# ------------------------------------------------------------------------------
#
# This script requires the `lattice` library, which can be installed using the
# following command if it is not already installed:
#
# install.packages('lattice')
#
# ------------------------------------------------------------------------------
#
# To run the script, set the R working directory to the directory that contains
# this script and type:
#
# source('plot_gm_analysis.R')

library(lattice)

source("gm_analysis.R")

###                            ###
### DECIDE WHICH PLOTS TO MAKE ###
###                            ###

MAKE_TDL_PLOTS <- TRUE

MAKE_GM_PLOTS <- TRUE

###                            ###
### MAKE TDL PLOTS, IF DESIRED ###
###                            ###

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

###                           ###
### MAKE GM PLOTS, IF DESIRED ###
###                           ###

if (MAKE_GM_PLOTS) {
    # Make gmc boxplots for each event
    gmc_event_boxplot <- bwplot(
        gmc ~ event | genotype,
        data = licor_files_no_outliers[['main_data']],
        ylim = c(0, MAX_GM),
        ylab = "Mesophyll conductance to CO2 (mol / m^2 / s / bar)"
    )
    x11()
    print(gmc_event_boxplot)

    # Make gmc boxplots for each rep
    gmc_rep_boxplot <- bwplot(
        gmc ~ event_replicate | genotype,
        data = licor_files_no_outliers[['main_data']],
        ylim = c(0, MAX_GM),
        ylab = "Mesophyll conductance to CO2 (mol / m^2 / s / bar)"
    )
    x11()
    print(gmc_rep_boxplot)

    # Make gmc bars for each event
    gmc_barchart <- barchart(
        gmc_avg ~ event,
        data = event_stats,
        ylim = c(0, MAX_GM),
        ylab = "Mesophyll conductance to CO2 (mol / m^2 / s / bar)",
        panel = function(x, y, ..., subscripts) {
            panel.barchart(x, y, subscripts = subscripts, ...)
            panel.arrows(x, y, x, event_stats[['gmc_upper']], length = 0.2, angle = 90, col = "black", lwd = 1)
            panel.arrows(x, y, x, event_stats[['gmc_lower']], length = 0.2, angle = 90, col = "black", lwd = 1)
        }
    )
    x11()
    print(gmc_barchart)

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
        Cc ~ event_replicate | genotype,
        data = licor_files_no_outliers[['main_data']],
        ylim = c(0, 275),
        ylab = "CO2 concentration in chloroplast (micromol / mol)"
    )
    x11()
    print(cc_rep_boxplot)

    # Make Cc bars for each event
    cc_barchart <- barchart(
        Cc_avg ~ event,
        data = event_stats,
        ylim = c(0, 275),
        ylab = "CO2 concentration in chloroplast (micromol / mol)",
        panel = function(x, y, ..., subscripts) {
            panel.barchart(x, y, subscripts = subscripts, ...)
            panel.arrows(x, y, x, event_stats[['Cc_upper']], length = 0.2, angle = 90, col = "black", lwd = 1)
            panel.arrows(x, y, x, event_stats[['Cc_lower']], length = 0.2, angle = 90, col = "black", lwd = 1)
        }
    )
    x11()
    print(cc_barchart)

    # Make drawdown boxplots for each event
    drawdown_event_boxplot <- bwplot(
        `Ci-Cc` ~ event | genotype,
        data = licor_files_no_outliers[['main_data']],
        ylim = c(0, 150),
        ylab = "CO2 drawdown to chloroplast (Ci - Cc) (micromol / mol)"
    )
    x11()
    print(drawdown_event_boxplot)

    # Make drawdown boxplots for each rep
    drawdown_boxplot <- bwplot(
        `Ci-Cc` ~ event_replicate | genotype,
        data = licor_files_no_outliers[['main_data']],
        ylim = c(0, 150),
        ylab = "CO2 drawdown to chloroplast (Ci - Cc) (micromol / mol)"
    )
    x11()
    print(drawdown_boxplot)

    # Make drawdown bars for each event
    drawdown_barchart <- barchart(
        `Ci-Cc_avg` ~ event,
        data = event_stats,
        ylim = c(0, 150),
        ylab = "CO2 drawdown to chloroplast (Ci - Cc) (micromol / mol)",
        panel = function(x, y, ..., subscripts) {
            panel.barchart(x, y, subscripts = subscripts, ...)
            panel.arrows(x, y, x, event_stats[['Ci-Cc_upper']], length = 0.2, angle = 90, col = "black", lwd = 1)
            panel.arrows(x, y, x, event_stats[['Ci-Cc_lower']], length = 0.2, angle = 90, col = "black", lwd = 1)
        }
    )
    x11()
    print(drawdown_barchart)

    # Make assimilation boxplots for each event
    assimilation_event_boxplot <- bwplot(
        A ~ event | genotype,
        data = licor_files_no_outliers[['main_data']],
        ylim = c(0, 50),
        ylab = "Net CO2 assimilation rate (micromol / m^2 / s)"
    )
    x11()
    print(assimilation_event_boxplot)

    # Make assimilation boxplots for each rep
    assimilation_boxplot <- bwplot(
        A ~ event_replicate | genotype,
        data = licor_files_no_outliers[['main_data']],
        ylim = c(0, 50),
        ylab = "Net CO2 assimilation rate (micromol / m^2 / s)"
    )
    x11()
    print(assimilation_boxplot)

    # Make assimilation bars for each event
    assimilation_barchart <- barchart(
        A_avg ~ event,
        data = event_stats,
        ylim = c(0, 50),
        ylab = "Net CO2 assimilation rate (micromol / m^2 / s)",
        panel = function(x, y, ..., subscripts) {
            panel.barchart(x, y, subscripts = subscripts, ...)
            panel.arrows(x, y, x, event_stats[['A_upper']], length = 0.2, angle = 90, col = "black", lwd = 1)
            panel.arrows(x, y, x, event_stats[['A_lower']], length = 0.2, angle = 90, col = "black", lwd = 1)
        }
    )
    x11()
    print(assimilation_barchart)

}
