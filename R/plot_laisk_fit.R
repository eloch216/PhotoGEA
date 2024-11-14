plot_laisk_fit <- function(
    fit_results,
    identifier_column_name,
    plot_type,
    cols = multi_curve_colors(),
    a_column_name = 'A',
    ci_column_name = 'Ci',
    ppfd_column_name = 'PPFD',
    ...
)
{
    if (length(fit_results) != 4 ||
        !is.exdf(fit_results$first_fit_parameters) ||
        !is.exdf(fit_results$first_fits) ||
        !is.exdf(fit_results$second_fit_parameters) ||
        !is.exdf(fit_results$second_fits)) {
        stop('fit_results must be the output from calculate_RL_laisk')
    }

    # Make sure the plot type is acceptable
    if (!tolower(plot_type) %in% c('first', 'second')) {
        stop('`plot_type` must be either `first` or `second`')
    }

    # Make sure the required variables are defined and have the correct units
    required_variables <- list()
    required_variables[[ci_column_name]]         <- unit_dictionary[['Ci']]
    required_variables[[identifier_column_name]] <- NA
    required_variables[[ppfd_column_name]]       <- NA

    # Don't throw an error if some columns are all NA
    check_required_variables(fit_results$first_fits, required_variables, check_NA = FALSE)

    required_variables[[a_column_name]] <- unit_dictionary[['A']]

    # Don't throw an error if some columns are all NA
    check_required_variables(fit_results$first_fits, required_variables, check_NA = FALSE)

    if (plot_type == 'first') {
        # Get the name of the fitted column
        a_fit_column_name <- paste0(a_column_name, '_fit')

        # Try converting light levels to numeric values
        if (!is.factor(laisk_results$first_fits[, ppfd_column_name])) {
            laisk_results$first_fits[, ppfd_column_name] <- sapply(
                laisk_results$first_fits[, ppfd_column_name],
                try_as_numeric
            )
        }

        # Make sure light levels are ordered
        laisk_results$first_fits[, ppfd_column_name] <- factor(
            laisk_results$first_fits[, ppfd_column_name],
            levels = sort(unique(laisk_results$first_fits[, ppfd_column_name]))
        )

        # Create the y-axis label
        y_label <- paste0(
            'Net CO2 assimilation rate [ ', laisk_results$first_fits$units[[a_column_name]], ' ]',
            '\n(open red circle: estimated intersection point)'
        )

        # Plot the individual fits
        lattice::xyplot(
            laisk_results$first_fits[, a_column_name] ~ laisk_results$first_fits[, ci_column_name] | laisk_results$first_fits[, identifier_column_name],
            group = laisk_results$first_fits[, ppfd_column_name],
            type = 'p',
            par.settings = list(superpose.symbol = list(col = cols, pch = 16)),
            auto.key = list(space = 'right'),
            curve_ids = laisk_results$first_fits[, identifier_column_name],
            xlab = paste('Intercellular CO2 concentration [', laisk_results$first_fits$units[[ci_column_name]], ']'),
            ylab = y_label,
            panel = function(...) {
                # Get info about this curve
                args <- list(...)
                curve_id <- args$curve_ids[args$subscripts][1]

                curve_data <-
                    laisk_results$first_fits[laisk_results$first_fits[, identifier_column_name] == curve_id, ]

                curve_parameters <-
                    laisk_results$second_fit_parameters[laisk_results$second_fit_parameters[, identifier_column_name] == curve_id, ]

                ppfd_vals <- unique(curve_data[, ppfd_column_name])

                # Plot fits at each light level
                for (qin in ppfd_vals) {
                    # Get data for this light level
                    ppfd_curve_data <- curve_data[curve_data[, ppfd_column_name] == qin, ]

                    # Plot the fit
                    lattice::panel.lines(
                        ppfd_curve_data[, a_fit_column_name] ~ ppfd_curve_data[, ci_column_name],
                        col = 'black'
                    )
                }

                # Plot measured points
                lattice::panel.xyplot(...)

                # Plot RL and Ci_star
                lattice::panel.points(
                    -curve_parameters[1, 'RL'] ~ curve_parameters[1, 'Ci_star'],
                    col = 'red',
                    pch = 1
                )
            },
            ...
        )
    } else {
        # Try converting light levels to numeric values
        if (!is.factor(laisk_results$second_fits[, ppfd_column_name])) {
            laisk_results$second_fits[, ppfd_column_name] <- sapply(
                laisk_results$second_fits[, ppfd_column_name],
                try_as_numeric
            )
        }

        # Make sure light levels are ordered
        laisk_results$second_fits[, ppfd_column_name] <- factor(
            laisk_results$second_fits[, ppfd_column_name],
            levels = sort(unique(laisk_results$second_fits[, ppfd_column_name]))
        )

        # Plot the intercepts and slopes
        lattice::xyplot(
            laisk_results$second_fits[, 'laisk_intercept'] ~ laisk_results$second_fits[, 'laisk_slope'] | laisk_results$second_fits[, identifier_column_name],
            group = laisk_results$second_fits[, ppfd_column_name],
            type = 'p',
            par.settings = list(superpose.symbol = list(col = cols, pch = 16)),
            auto.key = list(space = 'right'),
            curve_ids = laisk_results$second_fits[, identifier_column_name],
            xlab = paste('Laisk slope [', laisk_results$second_fits$units[['laisk_slope']], ']'),
            ylab = paste('Laisk intercept [', laisk_results$second_fits$units[['laisk_intercept']], ']'),
            panel = function(...) {
                # Get info about this curve
                args <- list(...)
                curve_id <- args$curve_ids[args$subscripts][1]

                curve_data <-
                    laisk_results$second_fits[laisk_results$second_fits[, identifier_column_name] == curve_id, ]

                lattice::panel.lines(curve_data[['laisk_intercept_fit']] ~ curve_data[['laisk_slope']], col = 'black')

                lattice::panel.xyplot(...)
            },
            ...
        )
    }
}
