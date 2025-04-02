# Equations referenced below come from the following sources:
# 1. Busch et al. Nat. Plants 6, 245–258 (2020) [https://doi.org/10.1038/s41477-020-0606-6]
calculate_gm_busch <- function(
    exdf_obj,
    e = -3, # ppt; isotopic fractionation during day respiration
    f = 11, # ppt; isotopic fractionation during photorespiration
    e_star_equation = 20,
    gm_type = 'dis',
    a_bar_column_name = 'a_bar',
    a_column_name = 'A',
    ci_column_name = 'Ci',
    co2_s_column_name = 'CO2_s',
    csurface_column_name = 'Csurface',
    delta_c13_r_column_name = 'delta_C13_r',
    delta_obs_growth_column_name = 'Delta_obs_growth',
    delta_obs_tdl_column_name = 'Delta_obs_tdl',
    gamma_star_column_name = 'Gamma_star_tl',
    rl_column_name = 'RL',
    total_pressure_column_name = 'total_pressure',
    t_column_name = 't'
)
{
    ##
    ## Preliminaries
    ##

    if (!is.exdf(exdf_obj)) {
        stop('exdf_obj must be an exdf object')
    }

    # Check inputs
    if (!e_star_equation %in% c(19, 20)) {
        stop('e_star_equation must be 19 or 20')
    }

    if (!gm_type %in% c('dis', 'con')) {
        stop('gm_type must be "dis" or "con"')
    }

    # Make sure the required variables are defined and have the correct units
    required_variables <- list()
    required_variables[[a_bar_column_name]]          <- 'ppt'
    required_variables[[a_column_name]]              <- 'micromol m^(-2) s^(-1)'
    required_variables[[ci_column_name]]             <- 'micromol mol^(-1)'
    required_variables[[co2_s_column_name]]          <- 'micromol mol^(-1)'
    required_variables[[csurface_column_name]]       <- 'micromol mol^(-1)'
    required_variables[[delta_c13_r_column_name]]    <- 'ppt'
    required_variables[[delta_obs_tdl_column_name]]  <- 'ppt'
    required_variables[[gamma_star_column_name]]     <- 'micromol mol^(-1)'
    required_variables[[rl_column_name]]             <- 'micromol m^(-2) s^(-1)'
    required_variables[[total_pressure_column_name]] <- 'bar'
    required_variables[[t_column_name]]              <- 'dimensionless'

    if (e_star_equation == 20) {
        required_variables[[delta_obs_growth_column_name]] <- 'ppt'
    }

    check_required_variables(exdf_obj, required_variables)

    # Extract some important columns
    A              <- exdf_obj[, a_column_name]              # micromol / m^2 / s
    a_bar          <- exdf_obj[, a_bar_column_name]          # ppt
    Ca             <- exdf_obj[, co2_s_column_name]          # micromol / mol
    Ci             <- exdf_obj[, ci_column_name]             # micromol / mol
    Cs             <- exdf_obj[, csurface_column_name]       # micromol / mol
    delta_Ca_meas  <- exdf_obj[, delta_c13_r_column_name]    # ppt
    Delta_obs_meas <- exdf_obj[, delta_obs_tdl_column_name]  # ppt
    Gamma_star     <- exdf_obj[, gamma_star_column_name]     # micromol / mol
    total_pressure <- exdf_obj[, total_pressure_column_name] # bar
    RL             <- exdf_obj[, rl_column_name]             # micromol / m^2 / s
    t              <- exdf_obj[, t_column_name]              # dimensionless

    Delta_obs_growth <- if (e_star_equation == 20) {
        exdf_obj[, delta_obs_growth_column_name] # ppt
    } else {
        NULL
    }

    # Get the values of some constants that are defined in `constants.R`
    a_b             <- ISOTOPE_CONSTANTS$a_b             # ppt
    a_s             <- ISOTOPE_CONSTANTS$a_s             # ppt
    a_m             <- ISOTOPE_CONSTANTS$a_m             # ppt
    b               <- ISOTOPE_CONSTANTS$b               # ppt
    delta_Ca_growth <- ISOTOPE_CONSTANTS$delta_Ca_growth # ppt

    ##
    ## Calculate isotopic fractionation values
    ##

    # Equations 19 and 20 from Busch et al. (2020) calculate the apparent
    # isotopic fractionation during day respiration under different assumptions
    e_star <- if (e_star_equation == 19) {
        # Equation 19 from Busch et al. (2020) uses a common simplification
        delta_Ca_meas - delta_Ca_growth # ppt
    } else {
        # Equation 20 from Busch et al. (2020) uses the difference in isotopic
        # composition between new and old assimilates
        (delta_Ca_meas - Delta_obs_meas) - (delta_Ca_growth - Delta_obs_growth) # ppt
    }

    # The total isotopic fractionation due to day respiration, as discussed just
    # before Equation 17 in Busch et al. (2020)
    e_prime <- e + e_star # ppt

    ##
    ## Calculate fractionation factors
    ##

    alpha_b <- 1 + b * 1e-3       # dimensionless; fractionation factor during Rubisco carboxylation
    alpha_e <- 1 + e_prime * 1e-3 # dimensionless; fractionation factor during day respiration
    alpha_f <- 1 + f * 1e-3       # dimensionless; fractionation factor during photorespiration

    # Un-numbered equation following Equation 12 in Busch et al. (2020):
    # fractionation factor during ???
    alpha_R <- 1 + (RL / A) * (e_prime * 1e-3 / alpha_e) # dimensionless

    ##
    ## Calculate ternary corrections
    ##

    # Factors used in subsequent calculations
    t_factor_1 <- 1 / (1 - t)       # dimensionless
    t_factor_2 <- (1 + t) / (1 - t) # dimensionless

    ##
    ## Calculate mesophyll conductance
    ##

    # Factor used in subsequent calculations
    Delta_i_term_1 <-
        a_b * (Ca - Cs) / Ca +
        a_s * (Cs - Ci) / Ca # ppt

    # The calculations for mesophyll conductance depend on whether or not day
    # respiration is assumed to be isotopically connected to the CBB cycle. Note
    # about units: For most equations used in calculate_gm_busch, gas
    # concentrations appear as ratios, so the results are the same whether
    # partial pressures or mole fractions are used. The only exceptions are
    # Equations 21 and 22 for the mesophyll conductance, where there is an
    # unmatched `Ca` in the denomenator. Here we convert from mole fraction
    # (micromol / mol) to partial pressure (bar) using the total pressure.
    if (gm_type == 'con') {
        # Here we assume that day respiration is isotopically connected to the
        # CBB cycle.

        # Factors used in subsequent calculations
        rd_a_factor <- RL / (A + RL)                                # dimensionless
        alpha_factor <- alpha_b / alpha_e                           # dimensionless
        rd_a_alpha_e_factor <- rd_a_factor * alpha_factor * e_prime # ppt

        # In this case, Delta_i is calculated using Equation 2 from Busch et al.
        # (2020) with gm assumed to be infinite; in other words, with Cc = Ci.
        # This assumption causes one factor in the equation to become zero:
        # a_m * (Ci - Cc) / Ca = 0. See paragraph just before Equation 21.
        Delta_i_term_2 <-
            b * Ci / Ca -
            rd_a_alpha_e_factor * ((Ci - Gamma_star) / Ca) -
            (Gamma_star / Ca) * (alpha_b / alpha_f) * f # ppt

        Delta_i <- t_factor_1 * Delta_i_term_1 + t_factor_2 * Delta_i_term_2 # ppt

        # Equation 21 from Busch et al. (2020)
        gm_bottom <- (Ca * 1e-6 * total_pressure) * (Delta_i - Delta_obs_meas) # ppt * bar
        gm_top <- t_factor_2 * A * 1e-6 * (b - a_m - rd_a_alpha_e_factor)      # ppt * mol / m^2 / s
        gmc <- gm_top / gm_bottom                                              # mol / m^2 / s / bar
    } else {
        # Here we assume that day respiration is isotopically disconnected from
        # the CBB cycle.

        # Factors used in subsequent calculations
        rd_a_factor <- RL / A                                       # dimensionless
        alpha_factor <- alpha_b / (alpha_e * alpha_R)               # dimensionless
        rd_a_alpha_e_factor <- rd_a_factor * alpha_factor * e_prime # ppt

        # In this case, Delta_i is calculated using Equation 13 from Busch et
        # al. (2020) with gm assumed to be infinite; in other words, with Cc =
        # Ci. This assumption causes one factor in the equation to become zero:
        # a_m * (Ci - Cc) / Ca = 0. Here we also assume h = 0 (see text
        # following Equation 22).
        Delta_i_term_2 <-
            b * Ci / Ca -
            rd_a_alpha_e_factor * (Ci / Ca) -
            (alpha_b / (alpha_f * alpha_R)) * (Gamma_star / Ca) * f # ppt

        Delta_i <- t_factor_1 * Delta_i_term_1 + t_factor_2 * Delta_i_term_2 # ppt

        # Equation 22 from Busch et al. (2020)
        gm_bottom <- (Ca * 1e-6 * total_pressure) * (Delta_i - Delta_obs_meas) # ppt * bar
        gm_top <- t_factor_2 * A * 1e-6 * (b - a_m - rd_a_alpha_e_factor)      # ppt * mol / m^2 / s
        gmc <- gm_top / gm_bottom                                              # mol / m^2 / s / bar
    }

    # Store the calculated quantities in the exdf object
    exdf_obj[, 'alpha_R'] <- alpha_R
    exdf_obj[, 'Delta_i'] <- Delta_i
    exdf_obj[, 'Delta_i_term_1']  <- Delta_i_term_1
    exdf_obj[, 'Delta_i_term_2']  <- Delta_i_term_2
    exdf_obj[, 'e_prime'] <- e_prime
    exdf_obj[, 'e_star'] <- e_star
    exdf_obj[, 'e_star_equation'] <- e_star_equation
    exdf_obj[, 'gmc'] <- gmc
    exdf_obj[, 'gm_bottom'] <- gm_bottom
    exdf_obj[, 'gm_top'] <- gm_top
    exdf_obj[, 'gm_type'] <- gm_type

    # # Document the columns that were added and return the exdf
    document_variables(
        exdf_obj,
        c('calculate_gm_busch', 'alpha_R',         'dimensionless'),
        c('calculate_gm_busch', 'Delta_i',         'ppt'),
        c('calculate_gm_busch', 'Delta_i_term_1',  'ppt'),
        c('calculate_gm_busch', 'Delta_i_term_2',  'ppt'),
        c('calculate_gm_busch', 'e_prime',         'ppt'),
        c('calculate_gm_busch', 'e_star',          'ppt'),
        c('calculate_gm_busch', 'e_star_equation', ''),
        c('calculate_gm_busch', 'gmc',             'mol m^(-2) s^(-1) bar^(-1)'),
        c('calculate_gm_busch', 'gm_bottom',       'ppt * micromol / mol'),
        c('calculate_gm_busch', 'gm_top',          'ppt * micromol / m^2 / s'),
        c('calculate_gm_busch', 'gm_type',         '')
    )
}
