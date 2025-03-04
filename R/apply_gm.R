apply_gm <- function(
    exdf_obj,
    gmc_at_25 = '', # mol / m^2 / s / bar (this value is sometimes being fitted)
    photosynthesis_type = 'C3',
    calculate_drawdown = TRUE,
    a_column_name = 'A',
    ca_column_name = 'Ca',
    ci_column_name = 'Ci',
    gmc_norm_column_name = 'gmc_norm',
    total_pressure_column_name = 'total_pressure',
    perform_checks = TRUE,
    return_exdf = TRUE
)
{
    if (perform_checks) {
        if (!is.exdf(exdf_obj)) {
            stop('apply_gm requires an exdf object')
        }

        # Make sure the required variables are defined and have the correct units
        required_variables <- list()
        required_variables[[a_column_name]]              <- 'micromol m^(-2) s^(-1)'
        required_variables[[ci_column_name]]             <- 'micromol mol^(-1)'
        required_variables[[gmc_norm_column_name]]       <- unit_dictionary(gmc_norm_column_name)
        required_variables[[total_pressure_column_name]] <- 'bar'

        if (calculate_drawdown) {
            required_variables[[ca_column_name]] <- 'micromol mol^(-1)'
        }

        flexible_param <- list(
            gmc_at_25 = gmc_at_25
        )

        required_variables <-
            require_flexible_param(required_variables, flexible_param)

        check_required_variables(exdf_obj, required_variables)
    }

    # Retrieve values of flexible parameters as necessary
    if (!value_set(gmc_at_25)) {gmc_at_25 <- exdf_obj[, 'gmc_at_25']}

    # Calculate gmc at leaf temperature
    gmc_tl <- gmc_at_25 * exdf_obj[, gmc_norm_column_name] # mol / m^2 / s / bar

    # Calculate internal CO2 concentration and partial pressure
    internal_c <- exdf_obj[, ci_column_name] - exdf_obj[, a_column_name] /
                    (gmc_tl * exdf_obj[, total_pressure_column_name]) # micromol / mol

    internal_c_pressure <- internal_c * exdf_obj[, total_pressure_column_name] # microbar

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
        exdf_obj[, c_column_name]  <- internal_c
        exdf_obj[, pc_column_name] <- internal_c_pressure
        exdf_obj[, 'gmc_at_25']    <- gmc_at_25
        exdf_obj[, 'gmc_tl']       <- gmc_tl

        exdf_obj <- document_variables(
            exdf_obj,
            c('apply_gm', c_column_name,          'micromol mol^(-1)'),
            c('apply_gm', pc_column_name,         'microbar'),
            c('apply_gm', 'gmc_at_25',            'mol m^(-2) s^(-1) bar^(-1)'),
            c('apply_gm', 'gmc_tl',               'mol m^(-2) s^(-1) bar^(-1)')
        )

        if (calculate_drawdown) {
            exdf_obj[, drawdown_cm_column_name] <-
                exdf_obj[, ci_column_name] - exdf_obj[, c_column_name]

            exdf_obj[, drawdown_cs_column_name] <-
                exdf_obj[, ca_column_name] - exdf_obj[, ci_column_name]

            exdf_obj <- document_variables(
                exdf_obj,
                c('apply_gm', drawdown_cm_column_name, 'micromol mol^(-1)'),
                c('apply_gm', drawdown_cs_column_name, 'micromol mol^(-1)')
            )
        }

        return(exdf_obj)
    } else {
        return(list(internal_c = internal_c))
    }
}
