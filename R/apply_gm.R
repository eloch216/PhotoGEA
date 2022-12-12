apply_gm <- function(
    licor_exdf,
    photosynthesis_type = 'C3',
    a_column_name = 'A',
    ca_column_name = 'Ca',
    ci_column_name = 'Ci',
    gmc_column_name = 'gmc',
    pa_column_name = 'Pa',
    deltapcham_column_name = 'DeltaPcham'
)
{
    if (!is.exdf(licor_exdf)) {
        stop('apply_gm requires an exdf object')
    }

    # Make sure the required variables are defined and have the correct units
    required_variables <- list()
    required_variables[[a_column_name]] <- 'micromol m^(-2) s^(-1)'
    required_variables[[ca_column_name]] <- 'micromol mol^(-1)'
    required_variables[[ci_column_name]] <- 'micromol mol^(-1)'
    required_variables[[gmc_column_name]] <- 'mol m^(-2) s^(-1) bar^(-1)'
    required_variables[[pa_column_name]] <- 'kPa'
    required_variables[[deltapcham_column_name]] <- 'kPa'

    check_required_variables(licor_exdf, required_variables)

    # Define new column names
    drawdown_m_column_name <- 'drawdown_m'
    drawdown_s_column_name <- 'drawdown_s'

    c_column_name <- if(toupper(photosynthesis_type) == 'C3') {
        'Cc'
    } else if (toupper(photosynthesis_type) == 'C4') {
        'Cm'
    } else {
        stop(paste('Unsupported photosynthesis_type:', photosynthesis_type))
    }

    # Make calculations
    licor_exdf[,c_column_name] <-
        licor_exdf[,ci_column_name] - licor_exdf[,a_column_name] /
            (licor_exdf[,gmc_column_name] * (licor_exdf[,pa_column_name] + licor_exdf[,deltapcham_column_name]) / 100)

    licor_exdf[,drawdown_m_column_name] <-
        licor_exdf[,ci_column_name] - licor_exdf[,c_column_name]

    licor_exdf[,drawdown_s_column_name] <-
        licor_exdf[,ca_column_name] - licor_exdf[,ci_column_name]

    # Document the columns that were added
    licor_exdf <- document_variables(
        licor_exdf,
        c('apply_gm', c_column_name,          'micromol mol^(-1)'),
        c('apply_gm', drawdown_m_column_name, 'micromol mol^(-1)'),
        c('apply_gm', drawdown_s_column_name, 'micromol mol^(-1)')
    )

    return(licor_exdf)
}
