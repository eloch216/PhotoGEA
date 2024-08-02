calculate_gm_ubierna <- function(
    exdf_obj,
    e = -3, # ppt; isotopic fractionation during day respiration
    f = 11, # ppt; isotopic fractionation during photorespiration
    a_bar_column_name = 'a_bar',
    a_column_name = 'A',
    ci_column_name = 'Ci',
    co2_s_column_name = 'CO2_s',
    csurface_column_name = 'Csurface',
    delta_c13_r_column_name = 'delta_C13_r',
    delta_obs_tdl_column_name = 'Delta_obs_tdl',
    gamma_star_column_name = 'Gamma_star',
    rl_column_name = 'RL',
    total_pressure_column_name = 'total_pressure',
    t_column_name = 't'
)
{
    ##
    ## Preliminaries
    ##

    if (!is.exdf(exdf_obj)) {
        stop('calculate_gm_ubierna requires an exdf object')
    }

    # Make sure the required variables are defined and have the correct units
    required_variables <- list()
    required_variables[[a_column_name]]              <- 'micromol m^(-2) s^(-1)'
    required_variables[[a_bar_column_name]]          <- 'ppt'
    required_variables[[ci_column_name]]             <- 'micromol mol^(-1)'
    required_variables[[co2_s_column_name]]          <- 'micromol mol^(-1)'
    required_variables[[csurface_column_name]]       <- 'micromol mol^(-1)'
    required_variables[[delta_c13_r_column_name]]    <- 'ppt'
    required_variables[[delta_obs_tdl_column_name]]  <- 'ppt'
    required_variables[[gamma_star_column_name]]     <- 'micromol mol^(-1)'
    required_variables[[total_pressure_column_name]] <- 'bar'
    required_variables[[rl_column_name]]             <- 'micromol m^(-2) s^(-1)'
    required_variables[[t_column_name]]              <- 'dimensionless'

    check_required_variables(exdf_obj, required_variables)

    # Extract some important columns
    A              <- exdf_obj[, a_column_name]              # micromol / m^2 / s
    a_bar          <- exdf_obj[, a_bar_column_name]          # ppt
    Ci             <- exdf_obj[, ci_column_name]             # micromol / mol
    CO2_s          <- exdf_obj[, co2_s_column_name]          # micromol / mol
    Csurface       <- exdf_obj[, csurface_column_name]       # micromol / mol
    delta_C13_r    <- exdf_obj[, delta_c13_r_column_name]    # ppt
    Delta_obs_tdl  <- exdf_obj[, delta_obs_tdl_column_name]  # ppt
    gstar          <- exdf_obj[, gamma_star_column_name]     # micromol / mol
    total_pressure <- exdf_obj[, total_pressure_column_name] # bar
    RL             <- exdf_obj[, rl_column_name]             # micromol / m^2 / s
    t              <- exdf_obj[, t_column_name]              # dimensionless

    # Get the values of some constants that are defined in `constants.R`
    b_prime_3    <- ISOTOPE_CONSTANTS$b               # ppt
    a_m          <- ISOTOPE_CONSTANTS$a_m             # ppt
    delta_growth <- ISOTOPE_CONSTANTS$delta_Ca_growth # ppt

    # Convert some quantities from ppm to partial pressures
    ppCO2_s <- (CO2_s * 1e-6) * total_pressure          # bar
    ppCO2_surface <- (Csurface * 1e-6) * total_pressure # bar
    ppCO2_i <- (Ci * 1e-6) * total_pressure             # bar
    Gamma_star <- (gstar * 1e-6) * total_pressure       # bar

    ##
    ## Calculate isotopic fractionation values
    ##

    # Equation 30
    e_star <- delta_C13_r - delta_growth # ppt

    # Equation 28 for the total isotopic fractionation due to day respiration
    e_prime <- e + e_star # ppt

    ##
    ## Calculate fractionation factors
    ##

    # Using un-numbered equations just after Equation 19
    alpha_b <- 1 + b_prime_3 * 1e-3 # dimensionless; fractionation factor during Rubisco carboxylation
    alpha_e <- 1 + e_prime * 1e-3   # dimensionless; fractionation factor during day respiration
    alpha_f <- 1 + f * 1e-3         # dimensionless; fractionation factor during photorespiration

    # Calculate some factors from t that will be used in later calculations
    t_factor_1 <- 1 + t
    t_factor_2 <- 1 / (1 - t)
    t_factor_3 <- (1 + t) / (1 - t)
    t_factor_4 <- (1 - t) / (1 + t)

    ##
    ## Calculate discrimination and mesophyll conductance
    ##

    # Equation 31 for discrimination that would occur if Ci = Cc in the absence
    # of any respiratory fractionation
    Delta_i <-
        t_factor_2 * a_bar +
        t_factor_2 * (t_factor_1 * b_prime_3 - a_bar) * (ppCO2_i / ppCO2_s) # ppt

    # Equation 33 for discrimination associated with day respiration
    Delta_e <-
        t_factor_3 * (alpha_b / alpha_e) * e_prime * (RL / (A + RL)) *
        (ppCO2_i - Gamma_star) / ppCO2_s # ppt

    # Equation 34 for discrimination associated with photorespiration
    Delta_f <- t_factor_3 * (alpha_b / alpha_f) * f * (Gamma_star / ppCO2_s) # ppt

    # Equation 44 for mesophyll conductance (using Evans and von Caemmerer
    # notation)
    Delta_difference <- Delta_i - Delta_e - Delta_f - Delta_obs_tdl # ppt

    equation_top <-
        t_factor_3 *
        (b_prime_3 - a_m - (alpha_b / alpha_e) * e_prime * (RL / (A + RL))) *
        (A / ppCO2_s) # ppt * micromol / m^2 / s / bar

    gmc <- equation_top / Delta_difference * 1e-6 # mol / m^2 / s / bar

    # Store the calculated quantities in the exdf object
    exdf_obj[, 'Delta_difference'] <- Delta_difference
    exdf_obj[, 'Delta_e'] <- Delta_e
    exdf_obj[, 'Delta_f'] <- Delta_f
    exdf_obj[, 'Delta_i'] <- Delta_i
    exdf_obj[, 'e_prime'] <- e_prime
    exdf_obj[, 'equation_top'] <- equation_top
    exdf_obj[, 'gmc'] <- gmc

    # # Document the columns that were added and return the exdf
    document_variables(
        exdf_obj,
        c('calculate_gm_ubierna', 'Delta_difference', 'ppt'),
        c('calculate_gm_ubierna', 'Delta_e',          'ppt'),
        c('calculate_gm_ubierna', 'Delta_f',          'ppt'),
        c('calculate_gm_ubierna', 'Delta_i',          'ppt'),
        c('calculate_gm_ubierna', 'e_prime',          'ppt'),
        c('calculate_gm_ubierna', 'equation_top',     'ppt * micromol / m^2 / s / bar'),
        c('calculate_gm_ubierna', 'gmc',              'mol m^(-2) s^(-1) bar^(-1)')
    )
}
