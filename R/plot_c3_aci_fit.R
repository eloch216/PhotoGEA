plot_c3_aci_fit <- function(
    fit_results,
    identifier_column_name,
    x_name,
    plot_operating_point = TRUE,
    plot_Ad = FALSE,
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

    # Don't throw an error if some columns are all NA
    check_required_variables(fit_results$fits_interpolated, required_variables, check_NA = FALSE)

    required_variables[[a_column_name]] <- unit_dictionary[['A']]

    # Don't throw an error if some columns are all NA
    check_required_variables(fit_results$fits, required_variables, check_NA = FALSE)

    required_variables <- list()
    required_variables[['operating_Ci']]       <- 'micromol mol^(-1)'
    required_variables[['operating_Cc']]       <- 'micromol mol^(-1)'
    required_variables[['operating_An_model']] <- unit_dictionary[['A']]

    # Don't throw an error if some columns are all NA
    check_required_variables(fit_results$parameters, required_variables, check_NA = FALSE)

    # Choose line settings
    assim_cols <- multi_curve_colors()[1:5]
    assim_cols[1] <- '#676767'

    line_settings <- list(
        col = assim_cols,
        lwd = c(4, 2, 2, 2, 2),
        lty = c(1, 5, 5, 5, 2)
    )

    # Specify the y-axis label
    y_label <- paste0(
        'Net CO2 assimilation rate [ ', fit_results$fits_interpolated$units[['An']],
        ' ]\n(filled black circles: measured data used for fits',
        '\nopen black circles: measured data excluded from fits'
    )

    y_label <- paste0(
        y_label,
        if (plot_operating_point) {
            '\nopen red circle: estimated operating point)'
        } else {
            ')'
        }
    )

    # Plot the fits, operating point, and raw data
    lattice::xyplot(
        if (plot_Ad) {
            An + Ac + Aj + Ap + Ad ~ fit_results$fits_interpolated[, x_name] | fit_results$fits_interpolated[, identifier_column_name]
        } else {
            An + Ac + Aj + Ap ~ fit_results$fits_interpolated[, x_name] | fit_results$fits_interpolated[, identifier_column_name]
        },
        data = fit_results$fits_interpolated$main_data,
        type = 'l',
        par.settings = list(superpose.line = line_settings),
        auto.key = list(space = 'right', lines = TRUE, points = FALSE),
        xlab = paste(x_name, '[', fit_results$fits_interpolated$units[[x_name]], ']'),
        ylab = y_label,
        curve_ids = fit_results$fits_interpolated[, identifier_column_name],
        panel = function(...) {
            # Get info about this curve
            args <- list(...)
            curve_id <- args$curve_ids[args$subscripts][1]

            curve_parameters <-
                fit_results$parameters[fit_results$parameters[, identifier_column_name] == curve_id, ]

            curve_data <-
                fit_results$fits[fit_results$fits[, identifier_column_name] == curve_id, ]

            curve_data_interpolated <-
                fit_results$fits_interpolated[fit_results$fits_interpolated[, identifier_column_name] == curve_id, ]

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

            # Plot the operating point, if desired
            if (plot_operating_point) {
                lattice::panel.points(
                    curve_parameters[1, 'operating_An_model'] ~ curve_parameters[1, operating_x],
                    col = 'red',
                    pch = 1
                )
            }
        },
        ...
    )
}
