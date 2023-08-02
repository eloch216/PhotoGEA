calculate_gm_ubierna <- function(licor_exdf)
{
    if (!is.exdf(licor_exdf)) {
        stop('calculate_gm_ubierna requires an exdf object')
    }

    # Make sure the required variables are defined and have the correct units
    required_variables <- list()
    required_variables$A             <- 'micromol m^(-2) s^(-1)'
    required_variables$Ci            <- 'micromol mol^(-1)'
    required_variables$CO2_s         <- 'micromol mol^(-1)'
    required_variables$Csurface      <- 'micromol mol^(-1)'
    required_variables$delta_C13_r   <- 'ppt'
    required_variables$E             <- 'mol m^(-2) s^(-1)'
    required_variables$gtc           <- 'mol m^(-2) s^(-1)'
    required_variables$Oxygen        <- '%'
    required_variables$Pa            <- 'kPa'
    required_variables$Rd            <- 'micromol m^(-2) s^(-1)'
    required_variables$Delta_obs_tdl <- 'ppt'

    check_required_variables(licor_exdf, required_variables)

    # Extract some important columns
    A             <- licor_exdf[, 'A']             # micromol / m^2 / s
    Ci            <- licor_exdf[, 'Ci']            # micromol / mol
    CO2_s         <- licor_exdf[, 'CO2_s']         # micromol / mol
    Csurface      <- licor_exdf[, 'Csurface']      # micromol / mol
    delta_C13_r   <- licor_exdf[, 'delta_C13_r']   # ppt
    E             <- licor_exdf[, 'E']             # mol / m^2 / s
    gtc           <- licor_exdf[, 'gtc']           # mol / m^2 / s
    Oxygen        <- licor_exdf[, 'Oxygen']        # percent
    Pa            <- licor_exdf[, 'Pa']            # kPa
    Rd            <- licor_exdf[, 'Rd']            # micromol / m^2 / s
    Delta_obs_tdl <- licor_exdf[, 'Delta_obs_tdl'] # ppt

    # Define some constants to avoid magic numbers in the equations
    ubierna_a_b <- 2.9  # 13C fractionation during diffusion through the leaf boundary layer (Ubierna, 2017)
    ubierna_a_s <- 4.4  # 13C fractionation due to diffusion in air (Ubierna, 2017)
    b_prime_3 <- 29     # where does this value come from? Ubierna (2017) uses 30, as do other sources
    delta_growth <- -8  # where does this value come from?
    alpha <- 0.5        # where does this value come from?
    specificity <- 2115 # where does this value come from?
    f <- 11.6           # fractionation factor during photorespiration; value from https://doi.org/10.1104/pp.108.130153
    ai <- 1.8           # check https://doi.org/10.1111/j.1365-3040.2012.02591.x

    # Make the calculations
    ppO2 <- (Oxygen * 1e-2) * (Pa * 1e-2)
    ppCO2_s <- (CO2_s * 1e-6) * (Pa * 1e-2)
    ppCO2_surface <- (Csurface * 1e-6) * (Pa * 1e-2)
    ppCO2_i <- (Ci * 1e-6) * (Pa * 1e-2)

    # `t` is a "ternary correction factor" calculated according to Equation
    # 7 of Farquhar & Cernusak, "Ternary effects on the gas exchange of
    # isotopologues of carbon dioxide" Plant, Cell & Environment 35,
    # 1221–1231 (2012)
    # (https://doi.org/10.1111/j.1365-3040.2012.02484.x).
    #
    # Equation 7 reads `t = alpha_ac * E / (2 * g_ac)`, where E is the
    # transpiration rate, `g_ac` is the conductance to diffusion of air in
    # CO2, and `alpha_ac` is determined by `1 / alpha_ac = 1 / (1 + a_bar)`,
    # where `a_bar` is the weighted fractionation across the boundary layer
    # and stomata in series.
    #
    # We can therefore write alpha_ac as 1 + a_bar.
    #
    # In turn, a_bar (13C) is calculated using the equation from Table 4 in
    # Ubierna et al. "Temperature response of mesophyll conductance in three
    # C4 species calculated with two methods: 18O discrimination and in
    # vitro Vpmax" New Phytologist 214, 66–80 (2017)
    # (https://doi.org/10.1111/nph.14359)

    a_bar <-
        (ubierna_a_b * (ppCO2_s - ppCO2_surface) + ubierna_a_s * (ppCO2_surface - ppCO2_i)) /
        (ppCO2_s - ppCO2_i)

    t = (1 + a_bar * 1e-3) * E / gtc / 2

    # What are these?
    e <- delta_C13_r - delta_growth
    Gamma_star <- alpha * ppO2 / specificity

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
        t_factor_3 * e * Rd * (ppCO2_i - Gamma_star) /
        ((A + Rd) * ppCO2_s)

    # Delta_f is determined by ?
    Delta_f <- t_factor_3 * (f * (Gamma_star / ppCO2_s))

    # calculate gm!
    delta_difference <- Delta_i - Delta_obs_tdl - Delta_e - Delta_f

    equation_top <-
        t_factor_3 *
        (b_prime_3 - ai - (e * Rd) / (A + Rd)) * (A / ppCO2_s)

    gmc <- equation_top / delta_difference / 1e6

    # Store the calculated quantities in the exdf object
    licor_exdf[, 'a_bar'] <- a_bar
    licor_exdf[, 'delta_difference'] <- delta_difference
    licor_exdf[, 'Delta_e'] <- Delta_e
    licor_exdf[, 'Delta_f'] <- Delta_f
    licor_exdf[, 'Delta_i'] <- Delta_i
    licor_exdf[, 'e'] <- e
    licor_exdf[, 'equation_top'] <- equation_top
    licor_exdf[, 'Gamma_star'] <- Gamma_star
    licor_exdf[, 'gmc'] <- gmc
    licor_exdf[, 't'] <- t

    # # Document the columns that were added and return the exdf
    document_variables(
        licor_exdf,
        c('calculate_gm_ubierna', 'a_bar',            'ppt'),
        c('calculate_gm_ubierna', 'delta_difference', 'ppt'),
        c('calculate_gm_ubierna', 'Delta_e',          'ppt'),
        c('calculate_gm_ubierna', 'Delta_f',          'ppt'),
        c('calculate_gm_ubierna', 'Delta_i',          'ppt'),
        c('calculate_gm_ubierna', 'e',                'ppt'),
        c('calculate_gm_ubierna', 'equation_top',     '?'),
        c('calculate_gm_ubierna', 'Gamma_star',       'bar'),
        c('calculate_gm_ubierna', 'gmc',              'mol m^(-2) s^(-1) bar^(-1)'),
        c('calculate_gm_ubierna', 't',                'dimensionless')
    )
}
