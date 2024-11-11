calculate_c3_variable_j <- function(
    exdf_obj,
    alpha_g,    # dimensionless      (this value is sometimes being fitted)
    alpha_s,    # dimensionless      (this value is sometimes being fitted)
    alpha_t,    # dimensionless      (this value is sometimes being fitted)
    Gamma_star, # micromol / mol     (this value is sometimes being fitted)
    RL_at_25,   # micromol / m^2 / s (at 25 degrees C; typically this value is being fitted)
    tau,        # dimensionless      (typically this value is being fitted)
    atp_use = 4.0,
    nadph_use = 8.0,
    a_column_name = 'A',
    ci_column_name = 'Ci',
    phips2_column_name = 'PhiPS2',
    qin_column_name = 'Qin',
    rl_norm_column_name = 'RL_norm',
    total_pressure_column_name = 'total_pressure',
    hard_constraints = 0,
    perform_checks = TRUE,
    return_exdf = TRUE
)
{
    if (perform_checks) {
        if (!is.exdf(exdf_obj)) {
            stop('calculate_c3_variable_j requires an exdf object')
        }

        # Make sure the required variables are defined and have the correct units
        required_variables <- list()
        required_variables[[a_column_name]]              <- 'micromol m^(-2) s^(-1)'
        required_variables[[ci_column_name]]             <- 'micromol mol^(-1)'
        required_variables[[phips2_column_name]]         <- 'dimensionless'
        required_variables[[qin_column_name]]            <- 'micromol m^(-2) s^(-1)'
        required_variables[[rl_norm_column_name]]        <- 'normalized to RL at 25 degrees C'
        required_variables[[total_pressure_column_name]] <- 'bar'

        flexible_param <- list(
            alpha_g = alpha_g,
            alpha_s = alpha_s,
            alpha_t = alpha_t,
            Gamma_star = Gamma_star,
            RL_at_25 = RL_at_25,
            tau = tau
        )

        required_variables <-
            require_flexible_param(required_variables, flexible_param)

        check_required_variables(exdf_obj, required_variables)
    }

    # Retrieve values of flexible parameters as necessary
    if (!value_set(alpha_g))    {alpha_g     <- exdf_obj[, 'alpha_g']}
    if (!value_set(alpha_s))    {alpha_s     <- exdf_obj[, 'alpha_s']}
    if (!value_set(alpha_t))    {alpha_t     <- exdf_obj[, 'alpha_t']}
    if (!value_set(Gamma_star)) {Gamma_star <- exdf_obj[, 'Gamma_star']}
    if (!value_set(RL_at_25))   {RL_at_25   <- exdf_obj[, 'RL_at_25']}
    if (!value_set(tau))        {tau        <- exdf_obj[, 'tau']}

    # Extract a few columns from the exdf object to make the equations easier to
    # read, converting units as necessary
    pressure   <- exdf_obj[, total_pressure_column_name]        # bar
    Ci         <- exdf_obj[, ci_column_name]                    # micromol / mol
    PCi        <- Ci * pressure                                 # microbar

    An     <- exdf_obj[, a_column_name]                         # micromol / m^2 / s
    PhiPS2 <- exdf_obj[, phips2_column_name]                    # dimensionless
    Qin    <- exdf_obj[, qin_column_name]                       # micromol / m^2 / s

    RL_tl <- RL_at_25 * exdf_obj[, rl_norm_column_name]         # micromol / m^2 / s

    # Make sure key inputs have reasonable values
    msg <- character()

    # Always check parameters that cannot be fit, and make sure we are not
    # mixing models
    mixed_j_coeff <- (abs(atp_use - 4) > 1e-10 || abs(nadph_use - 8) > 1e-10) &&
        (any(alpha_g > 0, na.rm = TRUE) || any(alpha_s > 0, na.rm = TRUE) || any(alpha_t > 0, na.rm = TRUE))

    if (any(PhiPS2 < 0, na.rm = TRUE))   {msg <- append(msg, 'PhiPS2 must be >= 0')}
    if (any(pressure < 0, na.rm = TRUE)) {msg <- append(msg, 'pressure must be >= 0')}
    if (any(Qin < 0, na.rm = TRUE))      {msg <- append(msg, 'Qin must be >= 0')}
    if (mixed_j_coeff)                   {msg <- append(msg, 'atp_use must be 4 and nadph_use must be 8 when alpha_g / alpha_s / alpha_t are nonzero')}

    # Optionally check whether Ci is reasonable
    if (hard_constraints >= 1) {
        if (any(Ci < 0, na.rm = TRUE)) {msg <- append(msg, 'Ci must be >= 0')}
    }

    # Optionally check reasonableness of parameters that can be fit
    if (hard_constraints >= 2) {
        if (any(alpha_g < 0 | alpha_g > 1, na.rm = TRUE))                    {msg <- append(msg, 'alpha_g must be >= 0 and <= 1')}
        if (any(alpha_s < 0 | alpha_s > 1, na.rm = TRUE))                    {msg <- append(msg, 'alpha_s must be >= 0 and <= 1')}
        if (any(alpha_t < 0 | alpha_t > 1, na.rm = TRUE))                    {msg <- append(msg, 'alpha_t must be >= 0 and <= 1')}
        if (any(alpha_g + 2 * alpha_t + 4 * alpha_s / 3 > 1, na.rm = TRUE))  {msg <- append(msg, 'alpha_g + 2 * alpha_t + 4 * alpha_s / 3 must be <= 1')}
        if (any(Gamma_star < 0, na.rm = TRUE))                               {msg <- append(msg, 'Gamma_star must be >= 0')}
        if (any(RL_at_25 < 0, na.rm = TRUE))                                 {msg <- append(msg, 'RL_at_25 must be >= 0')}
        if (any(tau < 0 | tau > 1, na.rm = TRUE))                            {msg <- append(msg, 'tau must be >= 0 and <= 1')}
    }

    msg <- paste(msg, collapse = '. ')

    # We only bypass these checks if !perform_checks && return_exdf
    if (perform_checks || !return_exdf) {
        if (msg != '') {
            stop(msg)
        }
    }

    # Get the effective value of Gamma_star
    Gamma_star_agt <- (1 - alpha_g + 2 * alpha_t) * Gamma_star * pressure # microbar

    # Calculate J_F (actual RuBP regeneration rate as estimated from
    # fluorescence) using Equation 5
    J_F <- tau * Qin * PhiPS2 # micromol / m^2 / s

    # Calculate gmc (mesophyll conductance to CO2) using Equation 7 from Harley
    # et al. (1992)
    AnRd <- An + RL_tl # micromol / m^2 / s

    gmc_top <- Gamma_star_agt * (J_F + (nadph_use + 16 * alpha_g - 8 * alpha_t + 8 * alpha_s) * AnRd) # microbar * micromol / m^2 / s

    gmc_bottom <- J_F - atp_use * AnRd # micromol / m^2 / s

    gmc <- An / (Ci - gmc_top / gmc_bottom) # mol / m^2 / s / bar

    # Calculate Cc
    Cc <- Ci - An / (gmc * pressure) # micromol / mol

    # Calculate the partial derivative of Cc with respect to A using Equation 10
    # from Harley et al. (1992)
    dCcdA <- (atp_use + nadph_use) * Gamma_star_agt * J_F / (J_F - atp_use * AnRd)^2 # bar m^2 s / mol

    # Indicate trust according to the slope criteria from Harley et al. (1992)
    harley_slope_trust <- dCcdA >= 10 & dCcdA <= 50

    if (return_exdf) {
        # Make a new exdf object from the calculated variables and make sure units
        # are included
        output <- exdf(data.frame(
            alpha_g = alpha_g,
            alpha_s = alpha_s,
            alpha_t = alpha_t,
            Gamma_star = Gamma_star,
            Gamma_star_agt = Gamma_star_agt,
            RL_at_25 = RL_at_25,
            tau = tau,
            RL_tl = RL_tl,
            J_F = J_F,
            gmc = gmc,
            Cc = Cc,
            dCcdA = dCcdA,
            harley_slope_trust = harley_slope_trust,
            atp_use = atp_use,
            nadph_use = nadph_use,
            c3_variable_j_msg = msg,
            stringsAsFactors = FALSE
        ))

        document_variables(
            output,
            c('calculate_c3_variable_j', 'alpha_g',            'dimensionless'),
            c('calculate_c3_variable_j', 'alpha_s',            'dimensionless'),
            c('calculate_c3_variable_j', 'alpha_t',            'dimensionless'),
            c('calculate_c3_variable_j', 'Gamma_star',         'micromol mol^(-1)'),
            c('calculate_c3_variable_j', 'Gamma_star_agt',      'microbar'),
            c('calculate_c3_variable_j', 'RL_at_25',           'micromol m^(-2) s^(-1)'),
            c('calculate_c3_variable_j', 'tau',                'dimensionless'),
            c('calculate_c3_variable_j', 'RL_tl',              'micromol m^(-2) s^(-1)'),
            c('calculate_c3_variable_j', 'J_F',                'micromol m^(-2) s^(-1)'),
            c('calculate_c3_variable_j', 'gmc',                'mol m^(-2) s^(-1) bar^(-1)'),
            c('calculate_c3_variable_j', 'Cc',                 'micromol mol^(-1)'),
            c('calculate_c3_variable_j', 'dCcdA',              'bar m^(2) s mol^(-1)'),
            c('calculate_c3_variable_j', 'harley_slope_trust', ''),
            c('calculate_c3_variable_j', 'atp_use',            'dimensionless'),
            c('calculate_c3_variable_j', 'nadph_use',          'dimensionless'),
            c('calculate_c3_variable_j', 'c3_variable_j_msg',  '')
        )
    } else {
        return(list(gmc = gmc, Cc = Cc))
    }
}
