# This script generates several plots from the calculations performed in
# response_curve_analysis.R. In fact, this script calls that one before making
# any plots.
#
# ------------------------------------------------------------------------------
#
# To run the script, set the R working directory to the directory that contains
# this script and type:
#
# source('plot_response_curve_analysis.R')
#
# ------------------------------------------------------------------------------
#
# For questions or comments, please contact Ed Lochocki (eloch@illinois.edu)

source('response_curve_analysis.R')

# Decide whether to view data frames along with the plots (can be useful for
# inspection to make sure the results look reasonable)
VIEW_DATA_FRAMES <- FALSE

if (VIEW_DATA_FRAMES) {
    View(all_samples)
    View(all_stats)
}

###                            ###
### PLOT RESPONSE CURVES TO CI ###
###                            ###

rc_caption <- "Average response curves for each genotype"

# Choose colors for the different genotypes to use when plotting average A-Ci
# curves. To see other available palettes, use one of the following commands:
#  display.brewer.all(colorblindFriendly = TRUE)
#  display.brewer.all(colorblindFriendly = FALSE)
rc_cols <- brewer.pal(9, "Set1")
rc_cols <- c("#000000", rc_cols[c(1:5,7:9)])
rc_cols <- rc_cols[1:length(unique(all_stats_subset[[GENOTYPE_COLUMN_NAME]]))]
rc_cols <- rev(rc_cols)

# Make a slightly different version of the color specification to use for the
# error bars
rc_error_cols <- rep(rc_cols, each=length(MEASUREMENT_NUMBERS))

# Set the line width to use for plotting average response curves. (This will
# only apply if type is 'l' or 'b' in the calls to xyplot below.)
line_width <- 1

# Plot the average A-Ci curves
aci_curves <- xyplot(
    A_avg ~ Ci_avg,
    group = genotype,
    data = all_stats_subset,
    type = 'b',
    pch = 16,
    lwd = line_width,
    auto = TRUE,
    grid = TRUE,
    main = rc_caption,
    xlab = "Intercellular [CO2] (ppm)",
    #xlab = "Intercellular [CO2] (ppm)\n(error bars: standard error of the mean)",
    ylab = "Net CO2 assimilation rate (micromol / m^2 / s)\n(error bars: standard error of the mean for same CO2 setpoint)",
    ylim = c(-10, 45),
    xlim = c(-50, 1300),
    par.settings=list(
        superpose.line=list(col=rc_cols),
        superpose.symbol=list(col=rc_cols)
    ),
    panel = function(x, y, ...) {
        panel.arrows(x, y, x, all_stats_subset[['A_upper']], length = 0.05, angle = 90, col = rc_error_cols, lwd = line_width)
        panel.arrows(x, y, x, all_stats_subset[['A_lower']], length = 0.05, angle = 90, col = rc_error_cols, lwd = line_width)
        #panel.arrows(x, y, all_stats_subset[['Ci_upper']], y, length = 0.05, angle = 90, col = rc_error_cols, lwd = line_width)
        #panel.arrows(x, y, all_stats_subset[['Ci_lower']], y, length = 0.05, angle = 90, col = rc_error_cols, lwd = line_width)
        panel.xyplot(x, y, ...)
    }
)

x11(width = 8, height = 6)
print(aci_curves)

# Plot the average ETR-Ci curves
eci_curves <- xyplot(
    ETR_avg ~ Ci_avg,
    group = genotype,
    data = all_stats_subset,
    type = 'b',
    pch = 16,
    lwd = line_width,
    auto = TRUE,
    grid = TRUE,
    main = rc_caption,
    xlab = "Intercellular [CO2] (ppm)",
    #xlab = "Intercellular [CO2] (ppm)\n(error bars: standard error of the mean)",
    ylab = "Electron transport rate (micromol / m^2 / s)\n(error bars: standard error of the mean for same CO2 setpoint)",
    ylim = c(0, 325),
    xlim = c(-50, 1300),
    par.settings=list(
        superpose.line=list(col=rc_cols),
        superpose.symbol=list(col=rc_cols)
    ),
    panel = function(x, y, ...) {
        panel.arrows(x, y, x, all_stats_subset[['ETR_upper']], length = 0.05, angle = 90, col = rc_error_cols, lwd = line_width)
        panel.arrows(x, y, x, all_stats_subset[['ETR_lower']], length = 0.05, angle = 90, col = rc_error_cols, lwd = line_width)
        #panel.arrows(x, y, all_stats_subset[['Ci_upper']], y, length = 0.05, angle = 90, col = rc_error_cols, lwd = line_width)
        #panel.arrows(x, y, all_stats_subset[['Ci_lower']], y, length = 0.05, angle = 90, col = rc_error_cols, lwd = line_width)
        panel.xyplot(x, y, ...)
    }
)

x11(width = 8, height = 6)
print(eci_curves)

###                                     ###
### PLOT ALL INDIVIDUAL RESPONSE CURVES ###
###                                     ###

ind_caption <- "Individual response curves for each genotype and rep"

num_reps <- length(unique(all_samples_subset[[REP_COLUMN_NAME]]))

# Choose colors for the different reps to use when plotting individual response
# curves. To see other available palettes, use one of the following commands:
#  display.brewer.all(colorblindFriendly = TRUE)
#  display.brewer.all(colorblindFriendly = FALSE)
ind_cols <- c(
    "#000000",
    brewer.pal(12, "Paired"),
    brewer.pal(8, "Set2"),
    brewer.pal(8, "Dark2")
)

# Plot each individual A-Ci curve, where each genotype will have multiple traces
# corresponding to different plants
multi_aci_curves <- xyplot(
    A ~ Ci | genotype,
    group = rep,
    data = all_samples_subset,
    type = 'b',
    pch = 20,
    auto.key = list(space = "right"),
    grid = TRUE,
    main = ind_caption,
    xlab = "Intercellular [CO2] (ppm)",
    ylab = "Net CO2 assimilation rate (micromol / m^2 / s)",
    ylim = c(-10, 50),
    xlim = c(-100, 1600),
    par.settings=list(
        superpose.line=list(col=ind_cols),
        superpose.symbol=list(col=ind_cols)
    )
)

x11(width = 8, height = 6)
print(multi_aci_curves)

# Plot each individual gsw-Ci curve, where each genotype will have multiple
# traces corresponding to different plants
multi_gsci_curves <- xyplot(
    gsw ~ Ci | genotype,
    group = rep,
    data = all_samples_subset,
    type = 'b',
    pch = 20,
    auto.key = list(space = "right"),
    grid = TRUE,
    main = ind_caption,
    xlab = "Intercellular [CO2] (ppm)",
    ylab = "Stomatal conductance to water (mol / m^2 / s)",
    ylim = c(0, 0.8),
    xlim = c(-100, 1600),
    par.settings=list(
        superpose.line=list(col=ind_cols),
        superpose.symbol=list(col=ind_cols)
    )
)

x11(width = 8, height = 6)
print(multi_gsci_curves)

###                                                    ###
### MAKE BOX-WHISKER PLOTS FOR FIRST MEASUREMENT POINT ###
###                                                    ###

boxplot_caption <- paste0(
    "Quartiles for measurement point ",
    POINT_FOR_BOX_PLOTS,
    "\n(where CO2_r_sp = ",
    all_samples_one_point[['CO2_r_sp']][1],
    ")"
)

a_boxplot <- bwplot(
    A ~ genotype,
    data = all_samples_one_point,
    ylab = "Net CO2 assimilation rate (micromol / m^2 / s)",
    ylim = c(0, 35),
    main = boxplot_caption,
    xlab = "Genotype"
)

x11(width = 6, height = 6)
print(a_boxplot)

phips2_boxplot <- bwplot(
    PhiPS2 ~ genotype,
    data = all_samples_one_point,
    ylab = "Photosystem II operating efficiency (dimensionless)",
    ylim = c(0, 0.4),
    main = boxplot_caption,
    xlab = "Genotype"
)

x11(width = 6, height = 6)
print(phips2_boxplot)


etr_boxplot <- bwplot(
    ETR ~ genotype,
    data = all_samples_one_point,
    ylab = "Electron transport rate (micromol / m^2 / s)",
    ylim = c(0, 275),
    main = boxplot_caption,
    xlab = "Genotype"
)

x11(width = 6, height = 6)
print(etr_boxplot)
