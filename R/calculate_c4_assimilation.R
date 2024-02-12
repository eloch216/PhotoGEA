calculate_c4_assimilation <- function(
    exdf_obj,
    Rd_at_25,                  # micromol / m^2 / s   (at 25 degrees C; typically this value is being fitted)
    Vcmax_at_25,               # micromol / m^2 / s   (at 25 degrees C; typically this value is being fitted)
    Vpmax_at_25,               # micromol / m^2 / s   (at 25 degrees C; typically this value is being fitted)
    Vpr,                       # micromol / m^2 / s   (typically this value is being fitted)
    POm = 210000,              # microbar             (typically this value is known from the experimental setup)
    gbs = 0.003,               # mol / m^2 / s / bar  (typically this value is fixed)
    Rm_frac = 0.5,             # dimensionless        (typically this value is fixed)
    alpha_psii = 0,                 # dimensionless        (typically this value is fixed)
    ao_column_name = 'ao',
    gamma_star_column_name = 'gamma_star',
    kc_column_name = 'Kc',
    ko_column_name = 'Ko',
    kp_column_name = 'Kp',
    pcm_column_name = 'PCm',
    rd_norm_column_name = 'Rd_norm',
    vcmax_norm_column_name = 'Vcmax_norm',
    vpmax_norm_column_name = 'Vpmax_norm',
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
        required_variables[[ao_column_name]]         <- 'dimensionless'
        required_variables[[gamma_star_column_name]] <- 'dimensionless'
        required_variables[[kc_column_name]]         <- 'microbar'
        required_variables[[ko_column_name]]         <- 'mbar'
        required_variables[[kp_column_name]]         <- 'microbar'
        required_variables[[pcm_column_name]]        <- 'microbar'
        required_variables[[rd_norm_column_name]]    <- 'normalized to Rd at 25 degrees C'
        required_variables[[vcmax_norm_column_name]] <- 'normalized to Vcmax at 25 degrees C'
        required_variables[[vpmax_norm_column_name]] <- 'normalized to Vpmax at 25 degrees C'

        flexible_param <- list(
            Rd_at_25 = Rd_at_25,
            Vcmax_at_25 = Vcmax_at_25,
            Vpmax_at_25 = Vpmax_at_25,
            Vpr = Vpr
        )

        required_variables <-
            require_flexible_param(required_variables, flexible_param)

        check_required_variables(exdf_obj, required_variables)
    }

    # Retrieve values of flexible parameters as necessary
    if (!value_set(Rd_at_25))    {Rd_at_25    <- exdf_obj[, 'Rd_at_25']}
    if (!value_set(Vcmax_at_25)) {Vcmax_at_25 <- exdf_obj[, 'Vcmax_at_25']}
    if (!value_set(Vpmax_at_25)) {Vpmax_at_25 <- exdf_obj[, 'Vpmax_at_25']}
    if (!value_set(Vpr))         {Vpr         <- exdf_obj[, 'Vpr']}

    # Extract a few columns from the exdf object to make the equations easier to
    # read, converting units as necessary
    Cm <- exdf_obj[, pcm_column_name]                # microbar
    Kc <- exdf_obj[, kc_column_name]                 # microbar
    Ko <- exdf_obj[, ko_column_name] * 1000          # microbar
    Kp <- exdf_obj[, kp_column_name]                 # microbar
    gamma_star <- exdf_obj[, gamma_star_column_name] # dimensionless
    ao <- exdf_obj[, ao_column_name]                 # dimensionless

    # Make sure key inputs have reasonable values
    msg <- character()

    if (any(Cm < 0, na.rm = TRUE))          {msg <- append(msg, 'PCm must be >= 0')}
    if (any(Kc < 0, na.rm = TRUE))          {msg <- append(msg, 'Kc must be >= 0')}
    if (any(Ko < 0, na.rm = TRUE))          {msg <- append(msg, 'Ko must be >= 0')}
    if (any(Kp < 0, na.rm = TRUE))          {msg <- append(msg, 'Kp must be >= 0')}
    if (any(gamma_star < 0, na.rm = TRUE))  {msg <- append(msg, 'gamma_star must be >= 0')}
    if (any(ao < 0, na.rm = TRUE))          {msg <- append(msg, 'ao must be >= 0')}
    if (any(Rd_at_25 < 0, na.rm = TRUE))    {msg <- append(msg, 'Rd_at_25 must be >= 0')}
    if (any(Vcmax_at_25 < 0, na.rm = TRUE)) {msg <- append(msg, 'Vcmax_at_25 must be >= 0')}
    if (any(Vpmax_at_25 < 0, na.rm = TRUE)) {msg <- append(msg, 'Vpmax_at_25 must be >= 0')}
    if (any(Vpr < 0, na.rm = TRUE))         {msg <- append(msg, 'Vpr must be >= 0')}

    msg <- paste(msg, collapse = '. ')

    # We only bypass these checks if !perform_checks && return_exdf
    if (perform_checks || !return_exdf) {
        if (msg != '') {
            stop(msg)
        }
    }

    # Apply temperature responses to Vcmax, Vpmax, Rd, and Rm, making use of
    # Table 4.1
    Vcmax_tl <- Vcmax_at_25 * exdf_obj[, vcmax_norm_column_name] # micromol / m^2 / s
    Vpmax_tl <- Vpmax_at_25 * exdf_obj[, vpmax_norm_column_name] # micromol / m^2 / s
    Rd_tl <- Rd_at_25 * exdf_obj[, rd_norm_column_name]          # micromol / m^2 / s
    Rm_tl <- Rm_frac * Rd_tl                                     # micromol / m^2 / s

    # Equations 4.17 and 4.19
    Vpc <- Cm * Vpmax_tl / (Cm + Kp)  # micromol / m^2 / s
    Vp <- pmin(Vpc, Vpr)              # micromol / m^2 / s

    # Calculate enzyme-limited assimilation rates. These are not explicitly
    # given by any equations in the textbook, but do appear as terms in some
    # later calculations.
    Apr <- Vpr - Rm_tl + gbs * Cm  # micromol / m^2 / s
    Apc <- Vpc - Rm_tl + gbs * Cm  # micromol / m^2 / s
    Ap <- Vp - Rm_tl + gbs * Cm    # micromol / m^2 / s
    Ar <- Vcmax_tl - Rd_tl         # micromol / m^2 / s

    # Calculate terms that appear in several of the next equations
    f1 <- alpha_psii / ao             # dimensionless
    f2 <- gbs * Kc * (1.0 + POm / Ko) # micromol / m^2 / s
    f3 <- gamma_star * Vcmax_tl       # micromol / m^2 / s
    f4 <- Kc / Ko                     # dimensionless

    # Equation 4.22 (here we use `qa` rather than `a`, where `q` stands for
    # `quadratic`)
    qa <- 1.0 - f1 * f4  # dimensionless

    # Equation 4.23 (here we use `qb` rather than `b` as in Equation 4.22)
    qb <- -(Ap + Ar + f2 + f1 * (f3 + Rd_tl * f4))  # micromol / m^2 / s

    # Equation 4.24 (here we use `qc` rather than `c` as in Equation 4.22)
    qc <- Ar * Ap - (f3 * gbs * POm + Rd_tl * f2)  # (micromol / m^2 / s)^2

    # Equation 4.21
    An <- sapply(seq_along(qa), function(i) {
        quadratic_root_minus(qa[i], qb[i], qc[i]) # micromol / m^2 / s
    })

    if (return_exdf) {
        # Make a new exdf object from the calculated variables and make sure units
        # are included
        output <- exdf(data.frame(
            Vcmax_tl = Vcmax_tl,
            Vpmax_tl = Vpmax_tl,
            Rd_tl = Rd_tl,
            Rm_tl = Rm_tl,
            Vpc = Vpc,
            Vpr = Vpr,
            Vp = Vp,
            Apc = Apc,
            Apr = Apr,
            Ap = Ap,
            Ar = Ar,
            An = An,
            c4_assimilation_msg = msg
        ))

        document_variables(
            output,
            c('calculate_c4_assimilation', 'Vcmax_tl',            'micromol m^(-2) s^(-1)'),
            c('calculate_c4_assimilation', 'Vpmax_tl',            'micromol m^(-2) s^(-1)'),
            c('calculate_c4_assimilation', 'Rd_tl',               'micromol m^(-2) s^(-1)'),
            c('calculate_c4_assimilation', 'Rm_tl',               'micromol m^(-2) s^(-1)'),
            c('calculate_c4_assimilation', 'Vpc',                 'micromol m^(-2) s^(-1)'),
            c('calculate_c4_assimilation', 'Vpr',                 'micromol m^(-2) s^(-1)'),
            c('calculate_c4_assimilation', 'Vp',                  'micromol m^(-2) s^(-1)'),
            c('calculate_c4_assimilation', 'Apc',                 'micromol m^(-2) s^(-1)'),
            c('calculate_c4_assimilation', 'Apr',                 'micromol m^(-2) s^(-1)'),
            c('calculate_c4_assimilation', 'Ap',                  'micromol m^(-2) s^(-1)'),
            c('calculate_c4_assimilation', 'Ar',                  'micromol m^(-2) s^(-1)'),
            c('calculate_c4_assimilation', 'An',                  'micromol m^(-2) s^(-1)'),
            c('calculate_c4_assimilation', 'c4_assimilation_msg', '')
        )
    } else {
        return(An)
    }
}
