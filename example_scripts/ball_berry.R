# Load required libraries
library(PhotoGEA)
library(lattice)
library(RColorBrewer)

# Define some important column names
A_COLUMN_NAME <- 'A'
BB_INDEX_COLUMN_NAME <- 'bb_index'
CA_COLUMN_NAME <- 'Ca'
CSURFACE_COLUMN_NAME <- 'Csurface'
DELTAPCHAM_COLUMN_NAME <- 'DeltaPcham'
E_COLUMN_NAME <- 'E'
GBW_COLUMN_NAME <- 'gbw'
GSW_COLUMN_NAME <- 'gsw'
H2O_S_COLUMN_NAME <- 'H2O_s'
PA_COLUMN_NAME <- 'Pa'
RHLEAF_COLUMN_NAME <- 'RHleaf'
TLEAF_COLUMN_NAME <- 'TleafCnd'

# Read Licor files, extract important columns, and combine the results into one
# exdf object
multi_file_info <- batch_read_licor_file(
    choose_input_licor_files(),
    preamble_data_rows = c(3, 5, 7, 9, 11, 13),
    variable_category_row = 14,
    variable_name_row = 15,
    variable_unit_row = 16,
    data_start_row = 17,
    timestamp_colname = 'time'
)

extracted_multi_file_info <- batch_extract_variables(
    multi_file_info,
    c(
        'obs',
        'time',
        'elapsed',
        'date',
        'hhmmss',
        'species',
        'plot',
        'instrument',
        A_COLUMN_NAME,
        CA_COLUMN_NAME,
        DELTAPCHAM_COLUMN_NAME,
        E_COLUMN_NAME,
        GBW_COLUMN_NAME,
        GSW_COLUMN_NAME,
        H2O_S_COLUMN_NAME,
        PA_COLUMN_NAME,
        TLEAF_COLUMN_NAME
    )
)

combined_info <- combine_exdf(extracted_multi_file_info)

# Calculate gas properties and use the result to determine the Ball-Berry index
combined_info <- calculate_gas_properties(
    combined_info,
    A_COLUMN_NAME,
    CA_COLUMN_NAME,
    DELTAPCHAM_COLUMN_NAME,
    E_COLUMN_NAME,
    GBW_COLUMN_NAME,
    GSW_COLUMN_NAME,
    H2O_S_COLUMN_NAME,
    PA_COLUMN_NAME,
    TLEAF_COLUMN_NAME
)

combined_info <- calculate_ball_berry_index(
    combined_info,
    A_COLUMN_NAME,
    RHLEAF_COLUMN_NAME,
    CSURFACE_COLUMN_NAME
)

# Make a new identifier based on the plot and instrument names
combined_info[,'plot_instrument'] <- paste(
    'plot',
    combined_info[,'plot'],
    '-',
    combined_info[,'instrument']
)

# Make a new identifier based on the species, plot, and instrument names
combined_info[,'species_plot_instrument'] <- paste(
    combined_info[,'species'],
    '- plot',
    combined_info[,'plot'],
    '-',
    combined_info[,'instrument']
)

# Check the data for any issues before proceeding with additional analysis
check_response_curve_data(
    combined_info[['main_data']],
    'species',
    'plot_instrument'
)

# Do Ball-Berry fitting
bb_results <- fit_ball_berry(
    combined_info,
    'species_plot_instrument',
    GSW_COLUMN_NAME,
    BB_INDEX_COLUMN_NAME
)

bb_parameters <- bb_results[['parameters']]
bb_fits <- bb_results[['fits']]

# Get averages, standard deviations, and standard errors for the Ball-Berry
# slope and intercept for each crop
bb_intercept_avg <- tapply(bb_parameters[['bb_intercept']], bb_parameters[['species']], mean)
bb_slope_avg <- tapply(bb_parameters[['bb_slope']], bb_parameters[['species']], mean)

bb_intercept_sd <- tapply(bb_parameters[['bb_intercept']], bb_parameters[['species']], sd)
bb_slope_sd <- tapply(bb_parameters[['bb_slope']], bb_parameters[['species']], sd)

bb_intercept_num <- tapply(bb_parameters[['bb_intercept']], bb_parameters[['species']], length)
bb_slope_num <- tapply(bb_parameters[['bb_slope']], bb_parameters[['species']], length)

bb_parameters_stats <- data.frame(
    bb_intercept_avg,
    bb_intercept_sd,
    bb_intercept_num,
    bb_intercept_stderr = bb_intercept_sd / sqrt(bb_intercept_num),
    bb_slope_avg,
    bb_slope_sd,
    bb_slope_num,
    bb_slope_stderr = bb_slope_sd / sqrt(bb_slope_num)
)
bb_parameters_stats[['species']] <- rownames(bb_parameters_stats)
rownames(bb_parameters_stats) <- NULL

# Choose limits
gsw_range <- c(0, 0.70)
bb_index_range <- c(0, 0.12)

# Plot the individual Ball-Berry curves
multi_bb_curves <- xyplot(
    gsw ~ bb_index | species,
    group = plot_instrument,
    data = combined_info[['main_data']],
    type = 'b',
    pch = 20,
    auto.key = list(space = 'right'),
    grid = TRUE,
    main = 'Individual Ball-Berry curves for each species and plot / instrument',
    xlab = 'Ball-Berry index (mmol / m^2 / s)',
    ylab = 'Stomatal conductance to H2O (mmol / m^2 / s)',
    ylim = gsw_range,
    xlim = bb_index_range,
    par.settings = list(
        superpose.line = list(col = default_colors),
        superpose.symbol = list(col = default_colors)
    )
)

x11(width = 8, height = 6)
print(multi_bb_curves)

# Plot the individual Ball-Berry curves and fits for each species
for (sp in unique(bb_fits[['species']])) {
    bb_intercept <- bb_parameters_stats[which(bb_parameters_stats[['species']]==sp), "bb_intercept_avg"]
    bb_slope <- bb_parameters_stats[which(bb_parameters_stats[['species']]==sp), "bb_slope_avg"]
    bb_intercept_err <- bb_parameters_stats[which(bb_parameters_stats[['species']]==sp), "bb_intercept_stderr"]
    bb_slope_err <- bb_parameters_stats[which(bb_parameters_stats[['species']]==sp), "bb_slope_stderr"]

    plot_obj <- xyplot(
        gsw + gsw_fit ~ bb_index | plot_instrument,
        data = bb_fits[which(bb_fits[['species']] == sp),],
        auto = TRUE,
        grid = TRUE,
        type = 'b',
        pch = 16,
        ylim = gsw_range,
        xlim = bb_index_range,
        xlab = 'Ball-Berry index (mmol / m^2 / s)',
        ylab = 'Stomatal conductance to H2O (mmol / m^2 / s)',
        main = paste0(
            "Ball-Berry fits for ", sp,
            ":\naverage intercept = ", bb_intercept, " +/- ", bb_intercept_err,
            "\naverage slope = ", bb_slope, " +/- ", bb_slope_err
        )
    )
    x11(width = 10, height = 6)
    print(plot_obj)
}
