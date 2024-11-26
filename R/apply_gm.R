apply_gm <- function(
    licor_exdf,
    photosynthesis_type = 'C3',
    calculate_drawdown = TRUE,
    a_column_name = 'A',
    ca_column_name = 'Ca',
    ci_column_name = 'Ci',
    gmc_column_name = 'gmc',
    total_pressure_column_name = 'total_pressure',
    perform_checks = TRUE,
    return_exdf = TRUE
)
{
    if (perform_checks) {
        if (!is.exdf(licor_exdf)) {
            stop('apply_gm requires an exdf object')
        }

        # Make sure the required variables are defined and have the correct units
        required_variables <- list()
        required_variables[[a_column_name]]              <- 'micromol m^(-2) s^(-1)'
        required_variables[[ci_column_name]]             <- 'micromol mol^(-1)'
        required_variables[[gmc_column_name]]            <- 'mol m^(-2) s^(-1) bar^(-1)'
        required_variables[[total_pressure_column_name]] <- 'bar'

        if (calculate_drawdown) {
            required_variables[[ca_column_name]] <- 'micromol mol^(-1)'
        }

        check_required_variables(licor_exdf, required_variables)
    }

    # Calculate internal CO2 concentration and partial pressure
    internal_c <-
        licor_exdf[, ci_column_name] - licor_exdf[, a_column_name] /
            (licor_exdf[, gmc_column_name] * licor_exdf[, total_pressure_column_name]) # micromol / mol

    internal_c_pressure <- internal_c * licor_exdf[, total_pressure_column_name] # microbar

    if (return_exdf) {
        # Define new column names
        drawdown_cm_column_name <- 'drawdown_cm'
        drawdown_cs_column_name <- 'drawdown_cs'

        c_column_name <- if(toupper(photosynthesis_type) == 'C3') {
            'Cc'
        } else if (toupper(photosynthesis_type) == 'C4') {
            'Cm'
        } else {
            stop(paste('Unsupported photosynthesis_type:', photosynthesis_type))
        }

        pc_column_name <- paste0('P', c_column_name)

        # Store results and document new columns
        licor_exdf[, c_column_name] <- internal_c

        licor_exdf[, pc_column_name] <- internal_c_pressure

        licor_exdf <- document_variables(
            licor_exdf,
            c('apply_gm', c_column_name,          'micromol mol^(-1)'),
            c('apply_gm', pc_column_name,         'microbar')
        )

        if (calculate_drawdown) {
            licor_exdf[, drawdown_cm_column_name] <-
                licor_exdf[, ci_column_name] - licor_exdf[, c_column_name]

            licor_exdf[, drawdown_cs_column_name] <-
                licor_exdf[, ca_column_name] - licor_exdf[, ci_column_name]

            licor_exdf <- document_variables(
                licor_exdf,
                c('apply_gm', drawdown_cm_column_name, 'micromol mol^(-1)'),
                c('apply_gm', drawdown_cs_column_name, 'micromol mol^(-1)')
            )
        }

        return(licor_exdf)
    } else {
        return(list(internal_c = internal_c))
    }
}
