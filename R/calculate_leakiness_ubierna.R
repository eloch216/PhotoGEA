calculate_leakiness_ubierna <- function(
    exdf_obj,
    e = -3, # ppt; isotopic fractionation during day respiration
    a_bar_column_name = 'a_bar',
    a_column_name = 'A',
    ci_column_name = 'Ci',
    co2_s_column_name = 'CO2_s',
    csurface_column_name = 'Csurface',
    delta_c13_r_column_name = 'delta_C13_r',
    delta_obs_tdl_column_name = 'Delta_obs_tdl',
    rl_column_name = 'RL',
    t_column_name = 't'
)
{
    if (!is.exdf(exdf_obj)) {
        stop('calculate_leakiness_ubierna requires an exdf object')
    }

    # Make sure the required variables are defined and have the correct units
    required_variables <- list()
    required_variables[[a_bar_column_name]]         <- 'ppt'
    required_variables[[a_column_name]]             <- 'micromol m^(-2) s^(-1)'
    required_variables[[ci_column_name]]            <- 'micromol mol^(-1)'
    required_variables[[co2_s_column_name]]         <- 'micromol mol^(-1)'
    required_variables[[csurface_column_name]]      <- 'micromol mol^(-1)'
    required_variables[[delta_c13_r_column_name]]   <- 'ppt'
    required_variables[[delta_obs_tdl_column_name]] <- 'ppt'
    required_variables[[rl_column_name]]            <- 'micromol m^(-2) s^(-1)'
    required_variables[[t_column_name]]             <- 'dimensionless'

    check_required_variables(exdf_obj, required_variables)

    # Extract some important columns
    A             <- exdf_obj[, a_column_name]             # micromol / m^2 / s
    a_bar         <- exdf_obj[, a_bar_column_name]         # ppt
    Ci            <- exdf_obj[, ci_column_name]            # micromol / mol
    CL            <- exdf_obj[, csurface_column_name]      # micromol / mol
    CO2_s         <- exdf_obj[, co2_s_column_name]         # micromol / mol
    delta_C13_r   <- exdf_obj[, delta_c13_r_column_name]   # ppt
    Delta_obs_tdl <- exdf_obj[, delta_obs_tdl_column_name] # ppt
    RL            <- exdf_obj[, rl_column_name]            # micromol / m^2 / s
    t             <- exdf_obj[, t_column_name]             # dimensionless

    # Get the values of some constants that are defined in `constants.R`
    b_prime_3        <- ISOTOPE_CONSTANTS$b               # ppt
    b_prime_4        <- ISOTOPE_CONSTANTS$b_prime_4       # ppt
    delta_13c_growth <- ISOTOPE_CONSTANTS$delta_Ca_growth # ppt
    s                <- ISOTOPE_CONSTANTS$s               # ppt

    # RLm is assumed to be proportional to RL
    Rm_frac <- 0.5
    RLm <- RL * Rm_frac # micromol / m^2 / s

    # Equation 21 from Ubierna et al. (2013). e_prime is the fractionation
    # during decarboxylation including measurement artefacts
    e_prime <- e + delta_C13_r - delta_13c_growth # ppt

    # The high-light version of Equation 16 from Ubierna et al. (2013). phi_i is
    # the leakiness estimated with the isotope method including respiratory and
    # photorespiratory fractionations
    phi_i_top <-
        ((1 - t) * Delta_obs_tdl * CO2_s - a_bar * (CO2_s - Ci)) / ((1 + t) * Ci) -
        b_prime_4 + e_prime * RLm / (A + 0.5 * RL) # ppt

    phi_i_bottom <- b_prime_3 - s + e_prime * (RLm / (A + 0.5 * RL) - RL / (A + RL)) # ppt

    phi_i <- phi_i_top / phi_i_bottom # dimensionless

    # The high-light version of Equation 15 from Ubierna et al. (2013). phi_is is
    # the leakiness estimated with the isotope method including respiratory and
    # photorespiratory fractionations and Cs
    phi_is <- phi_i # dimensionless

    # Equation 17 from Ubierna et al. (2013). phi_sim is the leakiness estimated
    # with the isotope method ignoring respiratory and photorespiratory
    # fractionations and Cs
    phi_sim_top <- Delta_obs_tdl * (1 - t) * CO2_s -
        a_bar * (CO2_s - Ci) - (1 + t) * Ci * b_prime_4 # ppt * micromol / mol

    phi_sim_bottom <- Ci * (1 + t) * (b_prime_3 - s) # ppt * micromol / mol

    phi_sim <- phi_sim_top / phi_sim_bottom # dimensionless

    # Store the calculated quantities in the exdf object
    exdf_obj[, 'e_prime'] <- e_prime
    exdf_obj[, 'phi_i'] <- phi_i
    exdf_obj[, 'phi_is'] <- phi_is
    exdf_obj[, 'phi_sim'] <- phi_sim

    # Document the columns that were added and return the exdf
    document_variables(
        exdf_obj,
        c('calculate_leakiness', 'e_prime', 'ppt'),
        c('calculate_leakiness', 'phi_i',   'dimensionless'),
        c('calculate_leakiness', 'phi_is',  'dimensionless'),
        c('calculate_leakiness', 'phi_sim', 'dimensionless')
    )
}
