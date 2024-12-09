calculate_wue <- function(
    exdf_obj,
    calculate_c3 = FALSE,
    a_column_name = 'A',
    ca_column_name = 'Ca',
    cc_column_name = 'Cc',
    ci_column_name = 'Ci',
    e_column_name = 'E',
    gmc_column_name = 'gmc_tl',
    gsw_column_name = 'gsw',
    h2o_a_column_name = 'H2O_s',
    h2o_i_column_name = 'H2O_i',
    total_pressure_column_name = 'total_pressure'
)
{
    if (!is.exdf(exdf_obj)) {
        stop('calculate_wue requires an exdf object')
    }

    # Make sure the required variables are defined and have the correct units
    required_variables <- list()
    required_variables[[a_column_name]]     <- 'micromol m^(-2) s^(-1)'
    required_variables[[ca_column_name]]    <- 'micromol mol^(-1)'
    required_variables[[ci_column_name]]    <- 'micromol mol^(-1)'
    required_variables[[e_column_name]]     <- 'mol m^(-2) s^(-1)'
    required_variables[[gsw_column_name]]   <- 'mol m^(-2) s^(-1)'
    required_variables[[h2o_a_column_name]] <- 'mmol mol^(-1)'
    required_variables[[h2o_i_column_name]] <- 'mmol mol^(-1)'

    if (calculate_c3) {
        required_variables[[gmc_column_name]]            <- 'mol m^(-2) s^(-1) bar^(-1)'
        required_variables[[cc_column_name]]             <- 'micromol mol^(-1)'
        required_variables[[total_pressure_column_name]] <- 'bar'
    }

    check_required_variables(exdf_obj, required_variables)

    # Calculate new columns
    exdf_obj[, 'lWUE'] <- exdf_obj[, a_column_name] / exdf_obj[, e_column_name]                # micromol CO2 / mol H2O
    exdf_obj[, 'iWUE'] <- exdf_obj[, a_column_name] / exdf_obj[, gsw_column_name]              # micromol CO2 / mol H2O
    exdf_obj[, 'Cia_ratio'] <- exdf_obj[, ci_column_name] / exdf_obj[, ca_column_name]         # dimensionless
    exdf_obj[, 'drawdown_sw'] <- exdf_obj[, h2o_i_column_name] - exdf_obj[, h2o_a_column_name] # mmol / mol

    if (calculate_c3) {
        exdf_obj[, 'g_ratio'] <- exdf_obj[, gmc_column_name] *
            exdf_obj[, total_pressure_column_name] / exdf_obj[, gsw_column_name] # dimensionless

        exdf_obj[, 'drawdown_ct'] <- exdf_obj[, ca_column_name] - exdf_obj[, cc_column_name] # micromol / mol
    }

    # Document the columns that were added and return the exdf object
    if (calculate_c3) {
        exdf_obj <- document_variables(
            exdf_obj,
            c('calculate_wue', 'g_ratio',     'dimensionless'),
            c('calculate_wue', 'drawdown_ct', 'micromol mol^(-1)')
        )
    }

    document_variables(
        exdf_obj,
        c('calculate_wue', 'lWUE',        'micromol CO2 / mol H2O'),
        c('calculate_wue', 'iWUE',        'micromol CO2 / mol H2O'),
        c('calculate_wue', 'Cia_ratio',   'dimensionless'),
        c('calculate_wue', 'drawdown_sw', 'mmol mol^(-1)')
    )
}
