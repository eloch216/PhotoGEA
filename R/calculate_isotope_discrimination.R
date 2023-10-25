calculate_isotope_discrimination <- function(
    exdf_obj,
    co2_r_column_name = 'CO2_r',
    co2_s_column_name = 'CO2_s',
    delta_C13_r_column_name = 'delta_C13_r',
    delta_C13_s_column_name = 'delta_C13_s',
    h2o_r_column_name = 'H2O_r',
    h2o_s_column_name = 'H2O_s',
    tdl_12C_r_column_name = 'calibrated_12c_r',
    tdl_12C_s_column_name = 'calibrated_12c_s'
)
{
    if (!is.exdf(exdf_obj)) {
        stop('calculate_isotope_discrimination requires an exdf object')
    }

    # Make sure the required variables are defined and have the correct units
    required_variables <- list()
    required_variables[[co2_r_column_name]] <- 'micromol mol^(-1)'
    required_variables[[co2_s_column_name]] <- 'micromol mol^(-1)'
    required_variables[[delta_C13_r_column_name]] <- 'ppt'
    required_variables[[delta_C13_s_column_name]] <- 'ppt'
    required_variables[[h2o_r_column_name]] <- 'mmol mol^(-1)'
    required_variables[[h2o_s_column_name]] <- 'mmol mol^(-1)'
    required_variables[[tdl_12C_r_column_name]] <- 'ppm'
    required_variables[[tdl_12C_s_column_name]] <- 'ppm'

    check_required_variables(exdf_obj, required_variables)

    # Extract some important columns
    CO2_r <- exdf_obj[, co2_r_column_name]             # micromol / mol
    CO2_s <- exdf_obj[, co2_s_column_name]             # micromol / mol
    delta_C13_r <- exdf_obj[, delta_C13_r_column_name] # ppt
    delta_C13_s <- exdf_obj[, delta_C13_s_column_name] # ppt
    H2O_r <- exdf_obj[, h2o_r_column_name]             # mmol / mol
    H2O_s <- exdf_obj[, h2o_s_column_name]             # mmol / mol
    tdl_12C_r <- exdf_obj[, tdl_12C_r_column_name]     # ppm
    tdl_12C_s <- exdf_obj[, tdl_12C_s_column_name]     # ppm

    # Correct CO2_s and CO2_r to account for the "dilution effect"
    CO2_r_corrected <- CO2_r / (1 - H2O_r * 1e-3) # micromol / mol
    CO2_s_corrected <- CO2_s / (1 - H2O_s * 1e-3) # micromol / mol

    # Calculate xsi using either the TDL or gas exchange CO2 concentrations
    xsi_gasex <- CO2_r_corrected / (CO2_r_corrected - CO2_s_corrected) # dimensionless
    xsi_tdl <- tdl_12C_r / (tdl_12C_r - tdl_12C_s)                     # dimensionless

    # Calculate Delta_obs using either the TDL or gas exchange CO2
    # concentrations. The factors of 1e3 in the code below are necessary because
    # we express isotope ratios in units of ppt.
    Delta_obs_gasex <-
        1e3 * xsi_gasex * (delta_C13_s - delta_C13_r) /
        (1e3 + delta_C13_s - xsi_gasex * (delta_C13_s - delta_C13_r)) # ppt

    Delta_obs_tdl <-
        1e3 * xsi_tdl * (delta_C13_s - delta_C13_r) /
        (1e3 + delta_C13_s - xsi_tdl * (delta_C13_s - delta_C13_r)) # ppt

    # Store the calculated quantities in the exdf object
    exdf_obj[, 'CO2_r_corrected'] <- CO2_r_corrected
    exdf_obj[, 'CO2_s_corrected'] <- CO2_s_corrected
    exdf_obj[, 'Delta_obs_gasex'] <- Delta_obs_gasex
    exdf_obj[, 'Delta_obs_tdl'] <- Delta_obs_tdl
    exdf_obj[, 'xsi_gasex'] <- xsi_gasex
    exdf_obj[, 'xsi_tdl'] <- xsi_tdl

    # Document the columns that were added and return the exdf
    document_variables(
        exdf_obj,
        c('calculate_isotope_discrimination', 'CO2_r_corrected', 'micromol mol^(-1)'),
        c('calculate_isotope_discrimination', 'CO2_s_corrected', 'micromol mol^(-1)'),
        c('calculate_isotope_discrimination', 'Delta_obs_gasex', 'ppt'),
        c('calculate_isotope_discrimination', 'Delta_obs_tdl',   'ppt'),
        c('calculate_isotope_discrimination', 'xsi_gasex',       'dimensionless'),
        c('calculate_isotope_discrimination', 'xsi_tdl',         'dimensionless')
    )
}
