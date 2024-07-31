plot_c3_aci_fit <- function(
    fit_results,
    identifier_column_name,
    x_name,
    plot_operating_point = TRUE,
    a_column_name = 'A',
    cc_column_name = 'Cc',
    ci_column_name = 'Ci',
    ...
)
{
    if (length(fit_results) != 3 || !is.exdf(fit_results$fits) ||
        !is.exdf(fit_results$parameters) || !is.exdf(fit_results$fits_interpolated)) {
        stop('fit_results must be the output from fit_c3_aci or fit_c3_variable_j')
    }

    # Make sure the x variable is acceptable
    if (!x_name %in% c(ci_column_name, cc_column_name)) {
        stop('The x_name must be the Ci or Cc column name')
    }

    # Get the appropriate operating point name
    operating_x <- if (x_name == ci_column_name) {
        'operating_Ci'
    } else {
        'operating_Cc'
    }

    # Make sure the required variables are defined and have the correct units
    required_variables <- list()
    required_variables[[cc_column_name]]         <- 'micromol mol^(-1)'
    required_variables[[ci_column_name]]         <- 'micromol mol^(-1)'
    required_variables[[identifier_column_name]] <- NA

    check_required_variables(fit_results$fits_interpolated, required_variables)
    
    required_variables[[a_column_name]] <- unit_dictionary[['A']]

    check_required_variables(fit_results$fits, required_variables)

    required_variables <- list()
    required_variables[['operating_Ci']]       <- 'micromol mol^(-1)'
    required_variables[['operating_Cc']]       <- 'micromol mol^(-1)'
    required_variables[['operating_An_model']] <- unit_dictionary[['A']]

    check_required_variables(fit_results$parameters, required_variables)

    # Choose line settings
    assim_cols <- multi_curve_colors()[1:4]
    assim_cols[1] <- '#676767'

    line_settings <- list(
        col = assim_cols,
        lwd = c(4, 2, 2, 2),
        lty = c(1, 5, 5, 5)
    )

    # Plot the fits, operating point, and raw data
    lattice::xyplot(
        An + Ac + Aj + Ap ~ fit_results$fits_interpolated[, x_name] | fit_results$fits_interpolated[, identifier_column_name],
        data = fit_results$fits_interpolated$main_data,
        type = 'l',
        par.settings = list(superpose.line = line_settings),
        auto.key = list(space = 'right', lines = TRUE, points = FALSE),
        xlab = paste(x_name, '[', fit_results$fits_interpolated$units[[x_name]], ']'),
        ylab = paste(
            'Net CO2 assimilation rate [', fit_results$fits_interpolated$units[['An']],
            ']\n(black circles: measured data; red circle: estimated operating point)'
        ),
        curve_ids = fit_results$fits_interpolated[, identifier_column_name],
        panel = function(...) {
            # Plot the fit lines
            lattice::panel.xyplot(...)

            # Get info about this curve
            args <- list(...)
            curve_id <- args$curve_ids[args$subscripts][1]

            curve_parameters <-
                fit_results$parameters[fit_results$parameters[, identifier_column_name] == curve_id, ]

            curve_data <-
                fit_results$fits[fit_results$fits[, identifier_column_name] == curve_id, ]

            # Plot the operating point, if desired
            if (plot_operating_point) {
                lattice::panel.points(
                    curve_parameters[1, 'operating_An_model'] ~ curve_parameters[1, operating_x],
                    col = 'red',
                    pch = 16
                )
            }

            # Plot the measured data points
            lattice::panel.points(
                curve_data[, a_column_name] ~ curve_data[, x_name],
                col = 'black',
                pch = 16
            )
        },
        ...
    )
}
