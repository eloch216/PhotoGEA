# Load required libraries
library(PhotoGEA)
library(lattice)
library(RColorBrewer)

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
        'A',
        'Ca',
        'DeltaPcham',
        'E',
        'gbw',
        'gsw',
        'H2O_s',
        'Pa',
        'TleafCnd'
    )
)

combined_info <- combine_exdf(extracted_multi_file_info)

# Calculate gas properties and use the result to determine the Ball-Berry index
combined_info <- calculate_gas_properties(combined_info)

combined_info[,'bb_index'] <-
    0.01 * combined_info[,'A'] * combined_info[,'RHleaf'] / combined_info[,'Csurface']

combined_info <- specify_variables(
        combined_info,
        c('calculated', 'bb_index', 'mmol m^(-2) s^(-1)')
)

# Make a new identifier based on the plot and instrument names
combined_info[,'plot_instrument'] <- paste(
    'plot',
    combined_info[,'plot'],
    '-',
    combined_info[,'instrument']
)

# Choose colors for the different reps to use when plotting individual response
# curves. To see other available palettes, use one of the following commands:
#  display.brewer.all(colorblindFriendly = TRUE)
#  display.brewer.all(colorblindFriendly = FALSE)
ind_cols <- c(
    '#000000',
    brewer.pal(12, 'Paired'),
    brewer.pal(8, 'Set2'),
    brewer.pal(8, 'Dark2')
)

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
    ylim = c(0, 0.70),
    xlim = c(0, 0.12),
    par.settings=list(
        superpose.line=list(col=ind_cols),
        superpose.symbol=list(col=ind_cols)
    )
)

x11(width = 8, height = 6)
print(multi_bb_curves)
