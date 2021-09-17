# This script generates several plots from the calculations performed in
# `c4_response_curve_analysis.R`. In fact, this script calls that one before
# making any plots.
#
# ------------------------------------------------------------------------------
#
# This script requires the `lattice`, and `RColorBrewer` libraries, which can be
# installed using the following commands if they are not already installed:
#
# install.packages('lattice')
# install.packages('RColorBrewer')
#
# ------------------------------------------------------------------------------
#
# To run the script, set the R working directory to the directory that contains
# this script and type:
#
# source('plot_c4_response_curve_analysis.R')

library(lattice)
library(RColorBrewer)

source('c4_induction_curve_analysis.R')

###                            ###
### PLOT RESPONSE CURVES TO CI ###
###                            ###

rc_caption <- "Average response curves for each event"

# Choose colors for the different events to use when plotting average A-Ci
# curves. To see other available palettes, use one of the following commands:
#  display.brewer.all(colorblindFriendly = TRUE)
#  display.brewer.all(colorblindFriendly = FALSE)
rc_cols <- brewer.pal(9, "Set1")
rc_cols <- c("#000000", rc_cols[c(1:5,7:9)])
rc_cols <- rc_cols[1:length(unique(all_stats_subset[[EVENT_COLUMN_NAME]]))]
rc_cols <- rev(rc_cols)

# Make a slightly different version of the color specification to use for the
# error bars
rc_error_cols <- rep(rc_cols, each=length(MEASUREMENT_NUMBERS))

# Set the line width to use for plotting average response curves. (This will
# only apply if type is 'l' or 'b' in the calls to xyplot below.)
line_width <- 2

# Plot the average induction curves
induction_curves <- xyplot(
    all_stats_subset[['A_avg']] ~ all_stats_subset[['elapsed_time']],
    group = all_stats_subset[[EVENT_COLUMN_NAME]],
    type = 'l',
    pch = 16,
    lwd = line_width,
    auto = TRUE,
    grid = TRUE,
    main = rc_caption,
    xlab = "Elapsed time (minutes)",
    #xlab = "Intercellular [CO2] (ppm)\n(error bars: standard error of the mean)",
    ylab = "Net CO2 assimilation rate (micromol / m^2 / s)\n(error bars: standard error of the mean for same CO2 setpoint)",
    ylim = c(-10, 60),
    #xlim = c(-50, 1300),
    par.settings=list(
        superpose.line=list(col=rc_cols),
        superpose.symbol=list(col=rc_cols)
    ),
    panel = function(x, y, ...) {
        #panel.arrows(x, y, x, all_stats_subset[['A_upper']], length = 0.05, angle = 90, col = rc_error_cols, lwd = line_width)
        #panel.arrows(x, y, x, all_stats_subset[['A_lower']], length = 0.05, angle = 90, col = rc_error_cols, lwd = line_width)
        panel.xyplot(x, y, ...)
    }
)

x11(width = 8, height = 6)
print(induction_curves)


###                                     ###
### PLOT ALL INDIVIDUAL INDUCTION CURVES ###
###                                     ###

ind_caption <- "Individual induction curves for each event and rep"

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

# Plot each individual response curve, where each event will have multiple traces
# corresponding to different plants
multi_induction_curves <- xyplot(
    all_samples_subset[['A']] ~ all_samples_subset[['elapsed_time']] | all_samples_subset[[EVENT_COLUMN_NAME]],
    group = all_samples_subset[[REP_COLUMN_NAME]],
    type = 'b',
    pch = 20,
    auto.key = list(space = "right"),
    grid = TRUE,
    main = ind_caption,
    xlab = "Elapsed time (minutes)",
    ylab = "Net CO2 assimilation rate (micromol / m^2 / s)",
    ylim = c(-10, 60),
    #xlim = c(-100, 1600),
    par.settings=list(
        superpose.line=list(col=ind_cols),
        superpose.symbol=list(col=ind_cols)
    )
)

x11(width = 8, height = 6)
print(multi_induction_curves)



###                                                    ###
### MAKE BOX-WHISKER PLOTS FOR last MEASUREMENT POINT ###
###                                                    ###

boxplot_caption <- paste0(
    "Quartiles for measurement point ",
    POINT_FOR_BOX_PLOTS,
    "\n(where CO2_r_sp = ",
    all_samples_one_point[['CO2_r_sp']][1],
    ")"
)

a_boxplot <- bwplot(
    all_samples_one_point[['A']] ~ all_samples_one_point[[EVENT_COLUMN_NAME]],
    ylab = "Net CO2 assimilation rate (micromol / m^2 / s)",
    ylim = c(0, 70),
    main = boxplot_caption,
    xlab = "Genotype"
)

x11(width = 6, height = 6)
print(a_boxplot)

phips2_boxplot <- bwplot(
    all_samples_one_point[['PhiPS2']] ~ all_samples_one_point[[EVENT_COLUMN_NAME]],
    ylab = "Photosystem II operating efficiency (dimensionless)",
    ylim = c(0, 0.4),
    main = boxplot_caption,
    xlab = "Genotype"
)

x11(width = 6, height = 6)
print(phips2_boxplot)


etr_boxplot <- bwplot(
    all_samples_one_point[['ETR']] ~ all_samples_one_point[[EVENT_COLUMN_NAME]],
    ylab = "Electron transport rate (micromol / m^2 / s)",
    ylim = c(0, 275),
    main = boxplot_caption,
    xlab = "Genotype"
)

x11(width = 6, height = 6)
print(etr_boxplot)

###                                             ###
### MAKE BAR CHARTS FOR FIRST MEASUREMENT POINT ###
###                                             ###

barchart_caption <- paste0(
    "Averages for measurement point ",
    POINT_FOR_BOX_PLOTS,
    "\n(where CO2_r_sp = ",
    all_stats_one_point[['CO2_r_sp_avg']][1],
    ")"
)

assimilation_barchart <- barchart(
    all_stats_one_point[['A_avg']] ~ all_stats_one_point[[EVENT_COLUMN_NAME]],
    ylim = c(0, 60),
    ylab = "Net CO2 assimilation rate (micromol / m^2 / s)",
    main = barchart_caption,
    panel = function(x, y, ..., subscripts) {
        panel.barchart(x, y, subscripts = subscripts, ...)
        panel.arrows(x, y, x, all_stats_one_point[['A_upper']], length = 0.2, angle = 90, col = "black", lwd = 1)
        panel.arrows(x, y, x, all_stats_one_point[['A_lower']], length = 0.2, angle = 90, col = "black", lwd = 1)
    }
)

x11(width = 6, height = 6)
print(assimilation_barchart)
