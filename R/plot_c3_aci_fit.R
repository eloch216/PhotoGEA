plot_c3_aci_fit <- function(
    fit_results,
    identifier_column_name,
    x_name,
    c_step = 1,
    plot_operating_point = TRUE,
    a_column_name = 'A',
    cc_column_name = 'Cc',
    ci_column_name = 'Ci',
    j_norm_column_name = 'J_norm',
    kc_column_name = 'Kc',
    ko_column_name = 'Ko',
    oxygen_column_name = 'oxygen',
    rd_norm_column_name = 'Rd_norm',
    total_pressure_column_name = 'total_pressure',
    vcmax_norm_column_name = 'Vcmax_norm',
    ...
)
{
    if (length(fit_results) != 2 || !is.exdf(fit_results$fits) || !is.exdf(fit_results$parameters)) {
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
    required_variables[['alpha_g']]                  <- unit_dictionary[['alpha_g']]
    required_variables[['Gamma_star']]               <- unit_dictionary[['Gamma_star']]
    required_variables[['J_at_25']]                  <- unit_dictionary[['J_at_25']]
    required_variables[['Rd_at_25']]                 <- unit_dictionary[['Rd_at_25']]
    required_variables[['Tp']]                       <- unit_dictionary[['Tp']]
    required_variables[['Vcmax_at_25']]              <- unit_dictionary[['Vcmax_at_25']]
    required_variables[[a_column_name]]              <- unit_dictionary[['A']]
    required_variables[[cc_column_name]]             <- 'micromol mol^(-1)'
    required_variables[[ci_column_name]]             <- 'micromol mol^(-1)'
    required_variables[[identifier_column_name]]     <- NA
    required_variables[[j_norm_column_name]]         <- 'normalized to J at 25 degrees C'
    required_variables[[kc_column_name]]             <- 'micromol mol^(-1)'
    required_variables[[ko_column_name]]             <- 'mmol mol^(-1)'
    required_variables[[oxygen_column_name]]         <- unit_dictionary[['oxygen']]
    required_variables[[rd_norm_column_name]]        <- 'normalized to Rd at 25 degrees C'
    required_variables[[total_pressure_column_name]] <- 'bar'
    required_variables[[vcmax_norm_column_name]]     <- 'normalized to Vcmax at 25 degrees C'

    check_required_variables(fit_results$fits, required_variables)

    required_variables <- list()
    required_variables[['operating_Ci']]       <- 'micromol mol^(-1)'
    required_variables[['operating_Cc']]       <- 'micromol mol^(-1)'
    required_variables[['operating_An_model']] <- unit_dictionary[['A']]

    check_required_variables(fit_results$parameters, required_variables)

    # Get FvCB model outputs with a finer Cc spacing
    new_assim <- do.call(rbind, by(fit_results$fits, fit_results$fits[, identifier_column_name], function(x) {
        # Get constant inputs
        atp_use       <- x[1, 'atp_use']
        nadph_use     <- x[1, 'nadph_use']
        curvature_cj  <- x[1, 'curvature_cj']
        curvature_cjp <- x[1, 'curvature_cjp']

        # Get Ci sequence to use
        cc_seq <- seq(
            min(x[, cc_column_name]),
            max(x[, cc_column_name]),
            by = c_step
        )

        # Restrict to key columns
        assim_input_col <- c(
            'alpha_g',
            'Gamma_star',
            'J_at_25',
            'Rd_at_25',
            'Tp',
            'Vcmax_at_25',
            cc_column_name,
            ci_column_name,
            j_norm_column_name,
            kc_column_name,
            ko_column_name,
            oxygen_column_name,
            rd_norm_column_name,
            total_pressure_column_name,
            vcmax_norm_column_name
        )

        assim_input <- x[, c(identifier_column_name, assim_input_col), TRUE]

        # Change the number of rows
        assim_input$main_data[seq_along(cc_seq), ] <- NA
        assim_input[, identifier_column_name] <- x[1, identifier_column_name]

        # Interpolate each column
        for (i in seq_along(assim_input_col)) {
            cn <- assim_input_col[i]

            assim_input[, cn] <- if (all(is.na(x[, cn]))) {
                NA
            } else {
                approx(x[, cc_column_name], x[, cn], cc_seq)[['y']]
            }
        }

        # Calculate assimilation rates
        assim <- calculate_c3_assimilation(
            assim_input,
            '', # alpha_g
            '', # Gamma_star
            '', # J_at_25
            '', # Rd_at_25
            '', # Tp
            '', # Vcmax_at_25
            atp_use,
            nadph_use,
            curvature_cj,
            curvature_cjp,
            cc_column_name,
            j_norm_column_name,
            kc_column_name,
            ko_column_name,
            oxygen_column_name,
            rd_norm_column_name,
            total_pressure_column_name,
            vcmax_norm_column_name
        )

        cbind(
            assim_input[, c(identifier_column_name, ci_column_name, cc_column_name), TRUE],
            assim
        )
    }))

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
        An + Ac + Aj + Ap ~ new_assim[, x_name] | new_assim[, identifier_column_name],
        data = new_assim$main_data,
        type = 'l',
        par.settings = list(superpose.line = line_settings),
        auto.key = list(space = 'right'),
        xlab = paste(x_name, '[', new_assim$units[[x_name]], ']'),
        ylab = paste(
            'Net CO2 assimilation rate [', new_assim$units[['An']],
            ']\n(black circles: measured data; red circle: estimated operating point)'
        ),
        curve_ids = new_assim[, identifier_column_name],
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
