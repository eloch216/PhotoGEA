plot_ball_berry_fit <- function(
    fit_results,
    identifier_column_name,
    bb_index_column_name = 'bb_index',
    gsw_column_name = 'gsw',
    ...
)
{
    if (length(fit_results) != 2 || !is.exdf(fit_results$fits) || !is.exdf(fit_results$parameters)) {
        stop('fit_results must be the output from fit_ball_berry')
    }

    gsw_fit_column_name <- paste0(gsw_column_name, '_fit')

    # Make sure the required variables are defined and have the correct units
    required_variables <- list()
    required_variables[[gsw_column_name]]        <- unit_dictionary[['gsw']]
    required_variables[[gsw_fit_column_name]]    <- unit_dictionary[['gsw']]
    required_variables[[bb_index_column_name]]   <- unit_dictionary[['bb_index']]
    required_variables[[identifier_column_name]] <- NA

    check_required_variables(fit_results$fits, required_variables)

    # Choose line settings
    line_settings <- list(
        col = '#676767',
        lwd = c(2, 0),
        lty = c(1, 0)
    )

    # Plot the fits, operating point, and raw data
    lattice::xyplot(
        fit_results$fits[, gsw_fit_column_name] + fit_results$fits[, gsw_column_name] ~
            fit_results$fits[, bb_index_column_name] | fit_results$fits[, identifier_column_name],
        type = 'l',
        par.settings = list(superpose.line = line_settings),
        auto.key = list(space = 'right', lines = TRUE, points = FALSE, text = 'Fit'),
        xlab = paste('Ball-Berry index [', fit_results$fits$units[[bb_index_column_name]], ']'),
        ylab = paste(
            'Stomatal conductance to H2O [', fit_results$fits$units[[gsw_column_name]],
            ']\n(filled black circles: measured data used for fits',
            '\nopen black circles: measured data excluded from fits)'
        ),
        curve_ids = fit_results$fits[, identifier_column_name],
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
                curve_data[used_for_fit, gsw_column_name] ~ curve_data[used_for_fit, bb_index_column_name],
                col = 'black',
                pch = 16
            )

            lattice::panel.points(
                curve_data[!used_for_fit, gsw_column_name] ~ curve_data[!used_for_fit, bb_index_column_name],
                col = 'black',
                pch = 1
            )
        },
        ...
    )
}
