calculate_gm_ubierna <- function(licor_exdf)
{
    if (!is.exdf(licor_exdf)) {
        stop("calculate_gm_ubierna requires an exdf object")
    }

    # Make sure the required variables are defined and have the correct units
    required_variables <- list()
    required_variables$CO2_s <- "micromol mol^(-1)"
    required_variables$H2O_s <- "mmol mol^(-1)"
    required_variables$CO2_r <- "micromol mol^(-1)"
    required_variables$H2O_r <- "mmol mol^(-1)"
    required_variables$Oxygen <- "%"
    required_variables$Pa <- "kPa"
    required_variables$Csurface <- "micromol mol^(-1)"
    required_variables$Ci <- "micromol mol^(-1)"
    required_variables$total_mixing_ratio_r <- "ppm"
    required_variables$total_mixing_ratio_s <- "ppm"
    required_variables$E <- "mol m^(-2) s^(-1)"
    required_variables$gtc <- "mol m^(-2) s^(-1)"
    required_variables$total_isotope_ratio_r <- "ppt"
    required_variables$total_isotope_ratio_s <- "ppt"
    required_variables$A <- "micromol m^(-2) s^(-1)"
    required_variables$respiration <- "micromol m^(-2) s^(-1)"

    check_required_variables(licor_exdf, required_variables)

    # Add some new columns to the Licor file in preparation for calculating
    # mesophyll conductance
    licor_exdf <- document_variables(
        licor_exdf,
        c("calculated", "Cs_licor",         "micromol mol^(-1)"),
        c("calculated", "Ce_licor",         "micromol mol^(-1)"),
        c("calculated", "ppO2",             "bar"),
        c("calculated", "gsc",              "mol m^(-2) s^(-1)"),
        c("calculated", "gbc",              "mol m^(-2) s^(-1)"),
        c("calculated", "ppCO2_s",          "bar"),
        c("calculated", "ppCO2_r",          "bar"),
        c("calculated", "ppCO2_surface",    "bar"),
        c("calculated", "ppCO2_i",          "bar"),
        c("calculated", "xsi_LICOR",        "?"),
        c("calculated", "xsi_TDL",          "?"),
        c("calculated", "a_bar",            "?"),
        c("calculated", "t",                "?"),
        c("calculated", "gamma_star",       "bar"),
        c("calculated", "e",                "?"),
        c("calculated", "delta_tdl",        "ppt"),
        c("calculated", "delta_i",          "ppt"),
        c("calculated", "delta_e",          "ppt"),
        c("calculated", "delta_f",          "ppt"),
        c("calculated", "equation_top",     "?"),
        c("calculated", "delta_difference", "ppt"),
        c("calculated", "gmc",              "mol m^(-2) s^(-1) bar^(-1)")
    )

    # Define some constants to avoid magic numbers in the equations
    ubierna_a_b <- 2.9                   # 13C fractionation during diffusion through the leaf boundary layer (Ubierna, 2017)
    ubierna_a_s <- 4.4                   # 13C fractionation due to diffusion in air (Ubierna, 2017)
    b_prime_3 <- 29                      # where does this value come from? Ubierna (2017) uses 30, as do other sources
    delta_growth <- -8                   # where does this value come from?
    alpha <- 0.5                         # where does this value come from?
    specificity <- 2115                  # where does this value come from?
    f <- 11.6                            # fractionation factor during photorespiration; value from https://doi.org/10.1104/pp.108.130153
    ai <- 1.8                            # check https://doi.org/10.1111/j.1365-3040.2012.02591.x

    # Make the calculations
    licor_exdf[['main_data']] <- within(licor_exdf[['main_data']], {
        Cs_licor <- 1e6 * CO2_s / (1e6 - H2O_s * 1e3)
        Ce_licor <- 1e6 * CO2_r / (1e6 - H2O_r * 1e3)
        ppO2 <- (Oxygen * 1e-2) * (Pa * 1e-2)

        ppCO2_s <- (CO2_s * 1e-6) * (Pa * 1e-2)
        ppCO2_r <- (CO2_r * 1e-6) * (Pa * 1e-2)
        ppCO2_surface <- (Csurface * 1e-6) * (Pa * 1e-2)
        ppCO2_i <- (Ci * 1e-6) * (Pa * 1e-2)

        xsi_LICOR <- CO2_r / (CO2_r - CO2_s)
        xsi_TDL <- total_mixing_ratio_r / (total_mixing_ratio_r - total_mixing_ratio_s)

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
        e <- total_isotope_ratio_r - delta_growth
        gamma_star <- alpha * ppO2 / specificity

        # Calculate some factors from t that will be used in later calculations
        t_factor_1 <- 1 + t
        t_factor_2 <- 1 / (1 - t)
        t_factor_3 <- (1 + t) / (1 - t)
        t_factor_4 <- (1 - t) / (1 + t)

        # delta_i is determined from Equation 10 in Farquhar et al. "On the
        # Relationship Between Carbon Isotope Discrimination and the
        # Intercellular Carbon Dioxide Concentration in Leaves" Functional
        # Plant Biol. 9, 121–137 (1982)
        # (https://doi.org/10.1071/pp9820121)
        #
        # where a modification of a -> a_bar / (1 - t) and
        # b -> b_prime_3 * (1 + t) / (1 - t)
        #
        # Does this make sense?
        delta_i <-
            t_factor_2 * a_bar +
            t_factor_2 * (t_factor_1 * b_prime_3 - a_bar) * ppCO2_i / ppCO2_s

        # delta_tdl is determined by ?
        delta_tdl <-
            1e3 * xsi_TDL * (total_isotope_ratio_s - total_isotope_ratio_r) /
            (1e3 + total_isotope_ratio_s - xsi_TDL * (total_isotope_ratio_s - total_isotope_ratio_r))

        # delta_e is determined by ?
        delta_e <-
            t_factor_3 * e * respiration * (ppCO2_i - gamma_star) /
            ((A + respiration) * ppCO2_s)

        # delta_f is determined by ?
        delta_f <- t_factor_3 * (f * (gamma_star / ppCO2_s))

        # calculate gm!
        delta_difference <- delta_i - delta_tdl - delta_e - delta_f

        equation_top <-
            t_factor_3 *
            (b_prime_3 - ai - (e * respiration) / (A + respiration)) * (A / ppCO2_s)

        gmc <- equation_top / delta_difference / 1e6
    })

    # remove some extra columns that we don't need to keep
    licor_exdf[,'t_factor_1'] <- NULL
    licor_exdf[,'t_factor_2'] <- NULL
    licor_exdf[,'t_factor_3'] <- NULL
    licor_exdf[,'t_factor_4'] <- NULL

    return(licor_exdf)
}
