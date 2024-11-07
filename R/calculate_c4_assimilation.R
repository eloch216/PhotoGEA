calculate_c4_assimilation <- function(
    exdf_obj,
    alpha_psii,  # dimensionless        (typically this value is fixed to 0)
    gbs,         # mol / m^2 / s / bar  (typically this value is fixed to 0.003)
    Jmax_at_25,  # micromol / m^2 / s   (at 25 degrees C; typically this value is being fitted)
    RL_at_25,    # micromol / m^2 / s   (at 25 degrees C; typically this value is being fitted)
    Rm_frac,     # dimensionless        (typically this value is fixed to 0.5 or 1.0)
    Vcmax_at_25, # micromol / m^2 / s   (at 25 degrees C; typically this value is being fitted)
    Vpmax_at_25, # micromol / m^2 / s   (at 25 degrees C; typically this value is being fitted)
    Vpr,         # micromol / m^2 / s   (typically this value is being fitted)
    absorptance = 0.85,
    f_spectral = 0.15,
    rho = 0.5,
    theta = 0.7,
    x_etr = 0.4,
    ao_column_name = 'ao',
    gamma_star_column_name = 'gamma_star',
    jmax_norm_column_name = 'Jmax_norm',
    kc_column_name = 'Kc',
    ko_column_name = 'Ko',
    kp_column_name = 'Kp',
    oxygen_column_name = 'oxygen',
    pcm_column_name = 'PCm',
    qin_column_name = 'Qin',
    rl_norm_column_name = 'RL_norm',
    total_pressure_column_name = 'total_pressure',
    vcmax_norm_column_name = 'Vcmax_norm',
    vpmax_norm_column_name = 'Vpmax_norm',
    hard_constraints = 0,
    perform_checks = TRUE,
    return_exdf = TRUE
)
{
    if (perform_checks) {
        if (!is.exdf(exdf_obj)) {
            stop('calculate_c4_assimilation requires an exdf object')
        }

        # Make sure the required variables are defined and have the correct units
        required_variables <- list()
        required_variables[[ao_column_name]]             <- 'dimensionless'
        required_variables[[gamma_star_column_name]]     <- 'dimensionless'
        required_variables[[jmax_norm_column_name]]      <- unit_dictionary[['Jmax_norm']]
        required_variables[[kc_column_name]]             <- 'microbar'
        required_variables[[ko_column_name]]             <- 'mbar'
        required_variables[[kp_column_name]]             <- 'microbar'
        required_variables[[oxygen_column_name]]         <- unit_dictionary[['oxygen']]
        required_variables[[pcm_column_name]]            <- 'microbar'
        required_variables[[qin_column_name]]            <- 'micromol m^(-2) s^(-1)'
        required_variables[[rl_norm_column_name]]        <- 'normalized to RL at 25 degrees C'
        required_variables[[total_pressure_column_name]] <- 'bar'
        required_variables[[vcmax_norm_column_name]]     <- 'normalized to Vcmax at 25 degrees C'
        required_variables[[vpmax_norm_column_name]]     <- 'normalized to Vpmax at 25 degrees C'

        flexible_param <- list(
            alpha_psii = alpha_psii,
            gbs = gbs,
            Jmax_at_25 = Jmax_at_25,
            RL_at_25 = RL_at_25,
            Rm_frac = Rm_frac,
            Vcmax_at_25 = Vcmax_at_25,
            Vpmax_at_25 = Vpmax_at_25,
            Vpr = Vpr
        )

        required_variables <-
            require_flexible_param(required_variables, flexible_param)

        check_required_variables(exdf_obj, required_variables)
    }

    # Retrieve values of flexible parameters as necessary
    if (!value_set(alpha_psii))  {alpha_psii  <- exdf_obj[, 'alpha_psii']}
    if (!value_set(gbs))         {gbs         <- exdf_obj[, 'gbs']}
    if (!value_set(Jmax_at_25))  {Jmax_at_25  <- exdf_obj[, 'Jmax_at_25']}
    if (!value_set(RL_at_25))    {RL_at_25    <- exdf_obj[, 'RL_at_25']}
    if (!value_set(Rm_frac))     {Rm_frac     <- exdf_obj[, 'Rm_frac']}
    if (!value_set(Vcmax_at_25)) {Vcmax_at_25 <- exdf_obj[, 'Vcmax_at_25']}
    if (!value_set(Vpmax_at_25)) {Vpmax_at_25 <- exdf_obj[, 'Vpmax_at_25']}
    if (!value_set(Vpr))         {Vpr         <- exdf_obj[, 'Vpr']}

    # Extract a few columns from the exdf object to make the equations easier to
    # read, converting units as necessary
    Cm         <- exdf_obj[, pcm_column_name]        # microbar
    Kc         <- exdf_obj[, kc_column_name]         # microbar
    Ko         <- exdf_obj[, ko_column_name] * 1000  # microbar
    Kp         <- exdf_obj[, kp_column_name]         # microbar
    gamma_star <- exdf_obj[, gamma_star_column_name] # dimensionless
    ao         <- exdf_obj[, ao_column_name]         # dimensionless
    Qin        <- exdf_obj[, qin_column_name]        # micromol / m^2 / s

    pressure <- exdf_obj[, total_pressure_column_name] # bar
    oxygen   <- exdf_obj[, oxygen_column_name]         # percent
    POm      <- oxygen * pressure * 1e4                # microbar

    # Make sure key inputs have reasonable values
    msg <- character()

    # Always check parameters that cannot be fit
    if (any(ao < 0, na.rm = TRUE))                {msg <- append(msg, 'ao must be >= 0')}
    if (any(gamma_star < 0, na.rm = TRUE))        {msg <- append(msg, 'gamma_star must be >= 0')}
    if (any(Kc < 0, na.rm = TRUE))                {msg <- append(msg, 'Kc must be >= 0')}
    if (any(Ko < 0, na.rm = TRUE))                {msg <- append(msg, 'Ko must be >= 0')}
    if (any(Kp < 0, na.rm = TRUE))                {msg <- append(msg, 'Kp must be >= 0')}
    if (any(oxygen < 0, na.rm = TRUE))            {msg <- append(msg, 'oxygen must be >= 0')}
    if (any(pressure < 0, na.rm = TRUE))          {msg <- append(msg, 'pressure must be >= 0')}
    if (any(x_etr < 0 | x_etr > 1, na.rm = TRUE)) {msg <- append(msg, 'x_etr must be >= 0 and <= 1')}

    # Optionally check whether PCm is reasonable
    if (hard_constraints >= 1) {
        if (any(Cm < 0, na.rm = TRUE)) {msg <- append(msg, 'PCm must be >= 0')}
    }

    # Optionally check reasonableness of parameters that can be fit
    if (hard_constraints >= 2) {
        if (any(alpha_psii < 0 | alpha_psii > 1, na.rm = TRUE)) {msg <- append(msg, 'alpha_psii must be >= 0 and <= 1')}
        if (any(gbs < 0, na.rm = TRUE))                         {msg <- append(msg, 'gbs must be >= 0')}
        if (any(Jmax_at_25 < 0, na.rm = TRUE))                  {msg <- append(msg, 'Jmax_at_25 must be >= 0')}
        if (any(RL_at_25 < 0, na.rm = TRUE))                    {msg <- append(msg, 'RL_at_25 must be >= 0')}
        if (any(Rm_frac < 0 | Rm_frac > 1, na.rm = TRUE))       {msg <- append(msg, 'Rm_frac must be >= 0 and <= 1')}
        if (any(Vcmax_at_25 < 0, na.rm = TRUE))                 {msg <- append(msg, 'Vcmax_at_25 must be >= 0')}
        if (any(Vpmax_at_25 < 0, na.rm = TRUE))                 {msg <- append(msg, 'Vpmax_at_25 must be >= 0')}
        if (any(Vpr < 0, na.rm = TRUE))                         {msg <- append(msg, 'Vpr must be >= 0')}
    }

    msg <- paste(msg, collapse = '. ')

    # We only bypass these checks if !perform_checks && return_exdf
    if (perform_checks || !return_exdf) {
        if (msg != '') {
            stop(msg)
        }
    }

    # Apply temperature responses to Vcmax, Vpmax, RL, RLm, and Jmax, making use
    # of Table 4.1
    Vcmax_tl <- Vcmax_at_25 * exdf_obj[, vcmax_norm_column_name] # micromol / m^2 / s
    Vpmax_tl <- Vpmax_at_25 * exdf_obj[, vpmax_norm_column_name] # micromol / m^2 / s
    RL_tl <- RL_at_25 * exdf_obj[, rl_norm_column_name]          # micromol / m^2 / s
    RLm_tl <- Rm_frac * RL_tl                                    # micromol / m^2 / s
    Jmax_tl <- Jmax_at_25 * exdf_obj[, jmax_norm_column_name]    # micromol / m^2 / s

    # Equations 4.17 and 4.19
    Vpc <- Cm * Vpmax_tl / (Cm + Kp)   # micromol / m^2 / s
    Vp <- pmin(Vpc, Vpr, na.rm = TRUE) # micromol / m^2 / s

    # Equations 2.14 and 2.15
    J_tl <- sapply(seq_along(Jmax_tl), function(i) {
        j_from_jmax(
            Jmax_tl[i],
            Qin[i],
            absorptance * (1 - f_spectral) * rho,
            theta
        )
    })

    # Calculate individual process-limited assimilation rates. These are not
    # explicitly given by any equations in the textbook, but do appear as terms
    # in some later calculations.
    Apr <- Vpr - RLm_tl + gbs * Cm              # micromol / m^2 / s
    Apc <- Vpc - RLm_tl + gbs * Cm              # micromol / m^2 / s
    Ap <- Vp - RLm_tl + gbs * Cm                # micromol / m^2 / s (Equation 4.25)
    Ar <- Vcmax_tl - RL_tl                      # micromol / m^2 / s (Equation 4.25)
    Ajm <- x_etr * J_tl / 2 - RLm_tl + gbs * Cm # micromol / m^2 / s (Equation 4.45)
    Ajbs <- (1 - x_etr) * J_tl / 3 - RL_tl      # micromol / m^2 / s (Equation 4.45)

    # Calculate terms that appear in several of the next equations
    f1 <- alpha_psii / ao                        # dimensionless
    f2 <- gbs * Kc * (1.0 + POm / Ko)            # micromol / m^2 / s
    f3 <- gamma_star * Vcmax_tl                  # micromol / m^2 / s
    f4 <- Kc / Ko                                # dimensionless
    f5 <- 7 * gamma_star / 3                     # dimensionless
    f6 <- (1 - x_etr) * J_tl / 3 + 7 * RL_tl / 3 # micromol / m^2 / s

    # Equation 4.22 (here we use `ea` rather than `a`, where `e` stands for
    # `enzyme`)
    ea <- 1.0 - f1 * f4  # dimensionless

    # Equation 4.23 (here we use `eb` rather than `b` as in Equation 4.22)
    eb <- -(Ap + Ar + f2 + f1 * (f3 + RL_tl * f4))  # micromol / m^2 / s

    # Equation 4.24 (here we use `ec` rather than `c` as in Equation 4.22)
    ec <- Ar * Ap - (f3 * gbs * POm + RL_tl * f2)  # (micromol / m^2 / s)^2

    # Equation 4.21 for the enzyme-limited assimilation rate
    Ac <- sapply(seq_along(ea), function(i) {
        quadratic_root_minus(ea[i], eb[i], ec[i]) # micromol / m^2 / s
    })

    # Equation 4.42 (here we use `la` rather than `a`, where `l` stands for
    # `light`)
    la <- 1.0 - f1 * f5 # dimensionless

    # Equation 4.43 (here we use `lb` rather than `b` as in Equation 4.43)
    lb <- -(Ajm + Ajbs + gbs * POm * f5 + gamma_star * f1 * f6)  # micromol / m^2 / s

    # Equation 4.45 (here we use `lc` rather than `c` as in Equation 4.44)
    lc <- Ajm * Ajbs - gamma_star * gbs * POm * f6  # (micromol / m^2 / s)^2

    # Equation 4.41 for the light-limited assimilation rate
    Aj <- sapply(seq_along(la), function(i) {
        quadratic_root_minus(la[i], lb[i], lc[i]) # micromol / m^2 / s
    })

    # Equation 4.47 for the overall assimilation rate
    An <- pmin(Ac, Aj, na.rm = TRUE) # micromol / m^2 / s

    if (return_exdf) {
        # Make a new exdf object from the calculated variables and make sure units
        # are included
        output <- exdf(data.frame(
            alpha_psii = alpha_psii,
            gbs = gbs,
            Jmax_at_25 = Jmax_at_25,
            Jmax_tl = Jmax_tl,
            J_tl = J_tl,
            RL_at_25 = RL_at_25,
            RL_tl = RL_tl,
            Rm_frac = Rm_frac,
            RLm_tl = RLm_tl,
            Vcmax_at_25 = Vcmax_at_25,
            Vcmax_tl = Vcmax_tl,
            Vpmax_at_25 = Vpmax_at_25,
            Vpmax_tl = Vpmax_tl,
            Vpr = Vpr,
            Vpc = Vpc,
            Vp = Vp,
            Apc = Apc,
            Apr = Apr,
            Ap = Ap,
            Ar = Ar,
            Ajm = Ajm,
            Ajbs = Ajbs,
            Ac = Ac,
            Aj = Aj,
            An = An,
            c4_assimilation_msg = msg,
            stringsAsFactors = FALSE
        ))

        document_variables(
            output,
            c('calculate_c4_assimilation', 'alpha_psii',          unit_dictionary$alpha_psii),
            c('calculate_c4_assimilation', 'gbs',                 unit_dictionary$gbs),
            c('calculate_c4_assimilation', 'Jmax_at_25',          'micromol m^(-2) s^(-1)'),
            c('calculate_c4_assimilation', 'Jmax_tl',             'micromol m^(-2) s^(-1)'),
            c('calculate_c4_assimilation', 'J_tl',                'micromol m^(-2) s^(-1)'),
            c('calculate_c4_assimilation', 'RL_at_25',            'micromol m^(-2) s^(-1)'),
            c('calculate_c4_assimilation', 'Rm_frac',             unit_dictionary$Rm_frac),
            c('calculate_c4_assimilation', 'Vcmax_at_25',         'micromol m^(-2) s^(-1)'),
            c('calculate_c4_assimilation', 'Vpmax_at_25',         'micromol m^(-2) s^(-1)'),
            c('calculate_c4_assimilation', 'Vpr',                 'micromol m^(-2) s^(-1)'),
            c('calculate_c4_assimilation', 'Vcmax_tl',            'micromol m^(-2) s^(-1)'),
            c('calculate_c4_assimilation', 'Vpmax_tl',            'micromol m^(-2) s^(-1)'),
            c('calculate_c4_assimilation', 'RL_tl',               'micromol m^(-2) s^(-1)'),
            c('calculate_c4_assimilation', 'RLm_tl',              'micromol m^(-2) s^(-1)'),
            c('calculate_c4_assimilation', 'Vpc',                 'micromol m^(-2) s^(-1)'),
            c('calculate_c4_assimilation', 'Vp',                  'micromol m^(-2) s^(-1)'),
            c('calculate_c4_assimilation', 'Apc',                 'micromol m^(-2) s^(-1)'),
            c('calculate_c4_assimilation', 'Apr',                 'micromol m^(-2) s^(-1)'),
            c('calculate_c4_assimilation', 'Ap',                  'micromol m^(-2) s^(-1)'),
            c('calculate_c4_assimilation', 'Ar',                  'micromol m^(-2) s^(-1)'),
            c('calculate_c4_assimilation', 'Ajm',                 'micromol m^(-2) s^(-1)'),
            c('calculate_c4_assimilation', 'Ajbs',                'micromol m^(-2) s^(-1)'),
            c('calculate_c4_assimilation', 'Ac',                  'micromol m^(-2) s^(-1)'),
            c('calculate_c4_assimilation', 'Aj',                  'micromol m^(-2) s^(-1)'),
            c('calculate_c4_assimilation', 'An',                  'micromol m^(-2) s^(-1)'),
            c('calculate_c4_assimilation', 'c4_assimilation_msg', '')
        )
    } else {
        return(An)
    }
}
