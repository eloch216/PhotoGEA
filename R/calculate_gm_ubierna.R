calculate_gm_ubierna <- function(
    licor_exdf,
    e = -3, # ppt; isotopic fractionation during day respiration
    f = 11  # ppt; isotopic fractionation during photorespiration
)
{
    if (!is.exdf(licor_exdf)) {
        stop('calculate_gm_ubierna requires an exdf object')
    }

    # Make sure the required variables are defined and have the correct units
    required_variables <- list()
    required_variables$A             <- 'micromol m^(-2) s^(-1)'
    required_variables$a_bar         <- 'ppt'
    required_variables$Ci            <- 'micromol mol^(-1)'
    required_variables$CO2_s         <- 'micromol mol^(-1)'
    required_variables$Csurface      <- 'micromol mol^(-1)'
    required_variables$delta_C13_r   <- 'ppt'
    required_variables$Delta_obs_tdl <- 'ppt'
    required_variables$Gamma_star    <- 'micromol mol^(-1)'
    required_variables$Pa            <- 'kPa'
    required_variables$Rd            <- 'micromol m^(-2) s^(-1)'
    required_variables$t             <- 'dimensionless'

    check_required_variables(licor_exdf, required_variables)

    # Extract some important columns
    A             <- licor_exdf[, 'A']             # micromol / m^2 / s
    a_bar         <- licor_exdf[, 'a_bar']         # ppt
    Ci            <- licor_exdf[, 'Ci']            # micromol / mol
    CO2_s         <- licor_exdf[, 'CO2_s']         # micromol / mol
    Csurface      <- licor_exdf[, 'Csurface']      # micromol / mol
    delta_C13_r   <- licor_exdf[, 'delta_C13_r']   # ppt
    Delta_obs_tdl <- licor_exdf[, 'Delta_obs_tdl'] # ppt
    Pa            <- licor_exdf[, 'Pa']            # kPa
    Rd            <- licor_exdf[, 'Rd']            # micromol / m^2 / s
    t             <- licor_exdf[, 't']             # dimensionless

    # Get the values of some constants that are defined in `constants.R`
    b_prime_3    <- ISOTOPE_CONSTANTS$b               # ppt
    ai           <- ISOTOPE_CONSTANTS$a_m             # ppt
    delta_growth <- ISOTOPE_CONSTANTS$delta_Ca_growth # ppt

    # Make the calculations; here we make use of 1 kPa = 0.01 bar
    ppCO2_s <- (CO2_s * 1e-6) * (Pa * 1e-2)                         # bar
    ppCO2_surface <- (Csurface * 1e-6) * (Pa * 1e-2)                # bar
    ppCO2_i <- (Ci * 1e-6) * (Pa * 1e-2)                            # bar
    Gamma_star <- (licor_exdf[, 'Gamma_star'] * 1e-6) * (Pa * 1e-2) # bar

    # Total isotopic fractionation due to day respiration
    e_prime <- e + delta_C13_r - delta_growth # ppt

    # Calculate some factors from t that will be used in later calculations
    t_factor_1 <- 1 + t
    t_factor_2 <- 1 / (1 - t)
    t_factor_3 <- (1 + t) / (1 - t)
    t_factor_4 <- (1 - t) / (1 + t)

    # Delta_i is determined from Equation 10 in Farquhar et al. "On the
    # Relationship Between Carbon Isotope Discrimination and the
    # Intercellular Carbon Dioxide Concentration in Leaves" Functional
    # Plant Biol. 9, 121–137 (1982)
    # (https://doi.org/10.1071/pp9820121)
    #
    # where a modification of a -> a_bar / (1 - t) and
    # b -> b_prime_3 * (1 + t) / (1 - t)
    #
    # Does this make sense?
    Delta_i <-
        t_factor_2 * a_bar +
        t_factor_2 * (t_factor_1 * b_prime_3 - a_bar) * ppCO2_i / ppCO2_s

    # Delta_e is determined by ?
    Delta_e <-
        t_factor_3 * e_prime * Rd * (ppCO2_i - Gamma_star) /
        ((A + Rd) * ppCO2_s)

    # Delta_f is determined by ?
    Delta_f <- t_factor_3 * (f * (Gamma_star / ppCO2_s))

    # calculate gm!
    delta_difference <- Delta_i - Delta_obs_tdl - Delta_e - Delta_f

    equation_top <-
        t_factor_3 *
        (b_prime_3 - ai - (e_prime * Rd) / (A + Rd)) * (A / ppCO2_s)

    gmc <- equation_top / delta_difference / 1e6

    # Store the calculated quantities in the exdf object
    licor_exdf[, 'delta_difference'] <- delta_difference
    licor_exdf[, 'Delta_e'] <- Delta_e
    licor_exdf[, 'Delta_f'] <- Delta_f
    licor_exdf[, 'Delta_i'] <- Delta_i
    licor_exdf[, 'e_prime'] <- e_prime
    licor_exdf[, 'equation_top'] <- equation_top
    licor_exdf[, 'gmc'] <- gmc

    # # Document the columns that were added and return the exdf
    document_variables(
        licor_exdf,
        c('calculate_gm_ubierna', 'delta_difference', 'ppt'),
        c('calculate_gm_ubierna', 'Delta_e',          'ppt'),
        c('calculate_gm_ubierna', 'Delta_f',          'ppt'),
        c('calculate_gm_ubierna', 'Delta_i',          'ppt'),
        c('calculate_gm_ubierna', 'e_prime',          'ppt'),
        c('calculate_gm_ubierna', 'equation_top',     '?'),
        c('calculate_gm_ubierna', 'gmc',              'mol m^(-2) s^(-1) bar^(-1)')
    )
}
