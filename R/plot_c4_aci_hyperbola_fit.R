plot_c4_aci_hyperbola_fit <- function(
    fit_results,
    identifier_column_name,
    a_column_name = 'A',
    ci_column_name = 'Ci',
    ...
)
{
    if (length(fit_results) != 3 || !is.exdf(fit_results$fits) ||
        !is.exdf(fit_results$parameters) || !(is.exdf(fit_results$fits_interpolated))) {
        stop('fit_results must be the output from fit_c4_aci_hyperbola')
    }

    # Make sure the required variables are defined and have the correct units
    required_variables <- list()
    required_variables[[ci_column_name]]         <- 'micromol mol^(-1)'
    required_variables[[identifier_column_name]] <- NA

    check_required_variables(fit_results$fits_interpolated, required_variables)

    required_variables[[a_column_name]] <- unit_dictionary[['A']]

    check_required_variables(fit_results$fits, required_variables)

    # Choose line settings
    assim_cols <- multi_curve_colors()[1:4]
    assim_cols[1] <- '#676767'

    line_settings <- list(
        col = assim_cols,
        lwd = c(4, 2, 2, 2),
        lty = c(1, 5, 5, 5)
    )

    # Plot the fits and raw data
    lattice::xyplot(
        An + Ainitial + Amax ~ fit_results$fits_interpolated[, ci_column_name] | fit_results$fits_interpolated[, identifier_column_name],
        data = fit_results$fits_interpolated$main_data,
        type = 'l',
        par.settings = list(superpose.line = line_settings),
        auto.key = list(space = 'right', lines = TRUE, points = FALSE),
        xlab = paste(ci_column_name, '[', fit_results$fits_interpolated$units[[ci_column_name]], ']'),
        ylab = paste(
            'Net CO2 assimilation rate [', fit_results$fits_interpolated$units[['An']],
            ']\n(filled black circles: measured data used for fits',
            '\nopen black circles: measured data excluded from fits)'
        ),
        curve_ids = fit_results$fits_interpolated[, identifier_column_name],
        panel = function(...) {

            # Get info about this curve
            args <- list(...)
            curve_id <- args$curve_ids[args$subscripts][1]

            curve_data <-
                fit_results$fits[fit_results$fits[, identifier_column_name] == curve_id, ]

            used_for_fit <- points_for_fitting(curve_data)

            # Plot the fit lines
            lattice::panel.xyplot(...)

            # Plot the measured data points
            lattice::panel.points(
                curve_data[used_for_fit, a_column_name] ~ curve_data[used_for_fit, x_name],
                col = 'black',
                pch = 16
            )

            lattice::panel.points(
                curve_data[!used_for_fit, a_column_name] ~ curve_data[!used_for_fit, x_name],
                col = 'black',
                pch = 1
            )
        },
        ...
    )
}
