estimate_licor_variance <- function(
    exdf_obj,
    sd_CO2_r, # micromol / mol
    sd_CO2_s, # micromol / mol
    sd_flow,  # micromol / s
    sd_H2O_r, # mmol / mol
    sd_H2O_s, # mmol / mol
    a_column_name = 'A',
    co2_r_column_name = 'CO2_r',
    co2_s_column_name = 'CO2_s',
    corrfact_column_name = 'CorrFact',
    flow_column_name = 'Flow',
    h2o_r_column_name = 'H2O_r',
    h2o_s_column_name = 'H2O_s',
    s_column_name = 'S'
)
{
    # Check inputs
    if (!is.exdf(exdf_obj)) {
        stop('estimate_licor_variance requires an exdf object')
    }

    # Make sure the required variables are defined and have the correct units.
    required_variables <- list()
    required_variables[[a_column_name]]        <- unit_dictionary('A')
    required_variables[[co2_r_column_name]]    <- unit_dictionary('CO2_r')
    required_variables[[co2_s_column_name]]    <- unit_dictionary('CO2_s')
    required_variables[[corrfact_column_name]] <- unit_dictionary('CorrFact')
    required_variables[[flow_column_name]]     <- unit_dictionary('Flow')
    required_variables[[h2o_r_column_name]]    <- unit_dictionary('H2O_r')
    required_variables[[h2o_s_column_name]]    <- unit_dictionary('H2O_s')
    required_variables[[s_column_name]]        <- unit_dictionary('S')

    check_required_variables(exdf_obj, required_variables)

    # Extract key variables to make the following equations simpler; also
    # convert some units along the way
    A        <- exdf_obj[, a_column_name]            # micromol / m^2 / s
    CO2_r    <- exdf_obj[, co2_r_column_name] * 1e-6 # mol / mol
    CO2_s    <- exdf_obj[, co2_s_column_name] * 1e-6 # mol / mol
    corrfact <- exdf_obj[, corrfact_column_name]     # dimensionless
    flow     <- exdf_obj[, flow_column_name]         # micromol / s
    H2O_r    <- exdf_obj[, h2o_r_column_name] * 1e-3 # mol / mol
    H2O_s    <- exdf_obj[, h2o_s_column_name] * 1e-3 # mol / mol
    s        <- exdf_obj[, s_column_name] * 1e-4     # m^2

    # Calculate a few factors that will be useful later
    front_f     <- corrfact * flow / s  # micromol / m^2 / s
    water_in_f  <- 1 - corrfact * H2O_r # dimensionless
    water_out_f <- 1 - corrfact * H2O_s # dimensionless
    ca_f        <- corrfact * CO2_s     # dimensionless

    # Calculate variance terms
    var_CO2_r <- (front_f * sd_CO2_r * 1e-6)^2                                      # (micromol / m^2 / s)^2
    var_CO2_s <- (-front_f * water_in_f / water_out_f * sd_CO2_s * 1e-6)^2          # (micromol / m^2 / s)^2
    var_flow  <- (A / flow * sd_flow)^2                                             # (micromol / m^2 / s)^2
    var_H2O_r <- (front_f * ca_f / water_out_f * sd_H2O_r * 1e-3)^2                 # (micromol / m^2 / s)^2
    var_H2O_s <- (-front_f * ca_f * water_in_f / water_out_f^2 * sd_H2O_s * 1e-3)^2 # (micromol / m^2 / s)^2

    # Calculate variance and stdev of the net CO2 assimilation rate
    var_A <- var_CO2_r + var_CO2_s + var_flow + var_H2O_r + var_H2O_s
    sd_A <- sqrt(var_A)

    # Store variances in the exdf object
    exdf_obj[, 'sd_A']      <- sd_A
    exdf_obj[, 'var_A']     <- var_A
    exdf_obj[, 'var_CO2_r'] <- var_CO2_r
    exdf_obj[, 'var_CO2_s'] <- var_CO2_s
    exdf_obj[, 'var_flow']  <- var_flow
    exdf_obj[, 'var_H2O_r'] <- var_H2O_r
    exdf_obj[, 'var_H2O_s'] <- var_H2O_s

    # Document the columns that were added and return the exdf
    document_variables(
        exdf_obj,
        c('estimate_licor_variance', 'sd_A',      'micromol m^(-2) s^(-1)'),
        c('estimate_licor_variance', 'var_A',     'micromol^2 m^(-4) s^(-2)'),
        c('estimate_licor_variance', 'var_CO2_r', 'micromol^2 m^(-4) s^(-2)'),
        c('estimate_licor_variance', 'var_CO2_s', 'micromol^2 m^(-4) s^(-2)'),
        c('estimate_licor_variance', 'var_flow',  'micromol^2 m^(-4) s^(-2)'),
        c('estimate_licor_variance', 'var_H2O_r', 'micromol^2 m^(-4) s^(-2)'),
        c('estimate_licor_variance', 'var_H2O_s', 'micromol^2 m^(-4) s^(-2)')
    )
}
