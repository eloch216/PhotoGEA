calculate_c3_assimilation <- function(
    exdf_obj,
    alpha_g,      # dimensionless      (this value is sometimes being fitted)
    Gamma_star,   # micromol / mol     (this value is sometimes being fitted)
    J_at_25,      # micromol / m^2 / s (at 25 degrees C; typically this value is being fitted)
    Rd_at_25,     # micromol / m^2 / s (at 25 degrees C; typically this value is being fitted)
    Tp,           # micromol / m^2 / s (typically this value is being fitted)
    Vcmax_at_25,  # micromol / m^2 / s (at 25 degrees C; typically this value is being fitted)
    POc = 210000, # microbar           (typically this value is known from the experimental setup)
    atp_use = 4.0,
    nadph_use = 8.0,
    curvature_cj = 1.0,
    curvature_cjp = 1.0,
    cc_column_name = 'Cc',
    j_norm_column_name = 'J_norm',
    kc_column_name = 'Kc',
    ko_column_name = 'Ko',
    rd_norm_column_name = 'Rd_norm',
    total_pressure_column_name = 'total_pressure',
    vcmax_norm_column_name = 'Vcmax_norm',
    perform_checks = TRUE,
    return_exdf = TRUE
)
{
    if (perform_checks) {
        if (!is.exdf(exdf_obj)) {
            stop('calculate_c3_assimilation requires an exdf object')
        }

        # Make sure the required variables are defined and have the correct units
        required_variables <- list()
        required_variables[[cc_column_name]]             <- 'micromol mol^(-1)'
        required_variables[[j_norm_column_name]]         <- 'normalized to J at 25 degrees C'
        required_variables[[kc_column_name]]             <- 'micromol mol^(-1)'
        required_variables[[ko_column_name]]             <- 'mmol mol^(-1)'
        required_variables[[rd_norm_column_name]]        <- 'normalized to Rd at 25 degrees C'
        required_variables[[total_pressure_column_name]] <- 'bar'
        required_variables[[vcmax_norm_column_name]]     <- 'normalized to Vcmax at 25 degrees C'

        flexible_param <- list(
            alpha_g = alpha_g,
            Gamma_star = Gamma_star,
            J_at_25 = J_at_25,
            Rd_at_25 = Rd_at_25,
            Tp = Tp,
            Vcmax_at_25 = Vcmax_at_25
        )

        required_variables <-
            require_flexible_param(required_variables, flexible_param)

        check_required_variables(exdf_obj, required_variables)

        # Make sure curvature parameters lie on [0,1]
        check_zero_one <- list(
            curvature_cj = curvature_cj,
            curvature_cjp = curvature_cjp
        )

        sapply(seq_along(check_zero_one), function(i) {
            if (any(check_zero_one[[i]] < 0 | check_zero_one[[i]] > 1)) {
                stop(paste(names(check_zero_one)[i], 'must be >= 0 and <= 1'))
            }
        })

        # Make sure the Cc values are all positive
        if (any(exdf_obj[, cc_column_name] <= 0)) {
            stop('All Cc values must be positive')
        }
    }

    # Retrieve values of flexible parameters as necessary
    if (!value_set(alpha_g))     {alpha_g     <- exdf_obj[, 'alpha_g']}
    if (!value_set(Gamma_star))  {Gamma_star  <- exdf_obj[, 'Gamma_star']}
    if (!value_set(J_at_25))     {J_at_25     <- exdf_obj[, 'J_at_25']}
    if (!value_set(Rd_at_25))    {Rd_at_25    <- exdf_obj[, 'Rd_at_25']}
    if (!value_set(Tp))          {Tp          <- exdf_obj[, 'Tp']}
    if (!value_set(Vcmax_at_25)) {Vcmax_at_25 <- exdf_obj[, 'Vcmax_at_25']}

    # Extract a few columns from the exdf object to make the equations easier to
    # read, converting units as necessary
    pressure <- exdf_obj[, total_pressure_column_name] # bar

    Cc <- exdf_obj[, cc_column_name] # micromol / mol
    PCc <- Cc * pressure             # microbar

    Kc <- exdf_obj[, kc_column_name] * pressure                 # microbar
    Ko <- exdf_obj[, ko_column_name] * pressure * 1000          # microbar
    Gamma_star <- Gamma_star * pressure                         # microbar

    Vcmax_tl <- Vcmax_at_25 * exdf_obj[, vcmax_norm_column_name] # micromol / m^2 / s
    Rd_tl <- Rd_at_25 * exdf_obj[, rd_norm_column_name]          # micromol / m^2 / s
    J_tl <- J_at_25 * exdf_obj[, j_norm_column_name]             # micromol / m^2 / s

    # Make sure key inputs have reasonable values
    msg <- character()

    if (any(alpha_g < 0 | alpha_g > 1, na.rm = TRUE)) {msg <- append(msg, 'alpha_g must be >= 0 and <= 1')}
    if (any(Cc < 0, na.rm = TRUE))                    {msg <- append(msg, 'Cc must be >= 0')}
    if (any(Gamma_star < 0, na.rm = TRUE))            {msg <- append(msg, 'Gamma_star must be >= 0')}
    if (any(J_at_25 < 0, na.rm = TRUE))               {msg <- append(msg, 'J_at_25 must be >= 0')}
    if (any(Kc < 0, na.rm = TRUE))                    {msg <- append(msg, 'Kc must be >= 0')}
    if (any(Ko < 0, na.rm = TRUE))                    {msg <- append(msg, 'Ko must be >= 0')}
    if (any(pressure < 0, na.rm = TRUE))              {msg <- append(msg, 'pressure must be >= 0')}
    if (any(Rd_at_25 < 0, na.rm = TRUE))              {msg <- append(msg, 'Rd_at_25 must be >= 0')}
    if (any(Tp < 0, na.rm = TRUE))                    {msg <- append(msg, 'Tp must be >= 0')}
    if (any(Vcmax_at_25 < 0, na.rm = TRUE))           {msg <- append(msg, 'Vcmax_at_25 must be >= 0')}

    msg <- paste(msg, collapse = '. ')

    # We only bypass these checks if !perform_checks && return_exdf
    if (perform_checks || !return_exdf) {
        if (msg != '') {
            stop(msg)
        }
    }

    # Rubisco-limited carboxylation (micromol / m^2 / s)
    Wc <- PCc * Vcmax_tl / (PCc + Kc * (1.0 + POc / Ko))

    # RuBP-regeneration-limited carboxylation (micromol / m^2 / s)
    Wj <- PCc * J_tl / (atp_use * PCc + nadph_use * Gamma_star)

    # TPU-limited carboxylation (micromol / m^2 / s)
    Wp <- PCc * 3 * Tp / (PCc - Gamma_star * (1 + 3 * alpha_g))
    Wp[PCc <= Gamma_star * (1 + 3 * alpha_g)] <- Inf

    # Co-limitation between Wc and Wj
    a_cj <- curvature_cj
    b_cj <- -(Wc + Wj)
    c_cj <- Wc * Wj

    Wcj <- sapply(seq_along(b_cj), function(i) {
        quadratic_root_min(a_cj, b_cj[i], c_cj[i]) # micromol / m^2 / s
    })

    # Co-limitation between Wcj and Wp. If Wp is infinite, then we have
    # Wp >> Wcj and Wp >> curvature_cjp, so the quadratic coefficients become
    # a_cjp = 0, b_cjp = -Wp, and c_cjp = Wcj * Wp. In that case, we have
    # 0 = -Wp * Wcjp + Wcj * Wp, whose solution is simply Wcjp = Wcj.
    a_cjp <- curvature_cjp
    b_cjp <- -(Wcj + Wp)
    c_cjp <- Wcj * Wp

    Wcjp <- sapply(seq_along(b_cjp), function(i) {
        if (is.infinite(Wp[i])) {
            Wcj[i] # micromol / m^2 / s
        } else {
            quadratic_root_min(a_cjp, b_cjp[i], c_cjp[i]) # micromol / m^2 / s
        }
    })

    # Calculate corresponding net CO2 assimilations by accounting for
    # photorespiration and day respiration
    photo_resp_factor <- 1.0 - Gamma_star / PCc # dimensionless
    Ac <- photo_resp_factor * Wc - Rd_tl
    Aj <- photo_resp_factor * Wj - Rd_tl
    Ap <- photo_resp_factor * Wp - Rd_tl
    An <- photo_resp_factor * Wcjp - Rd_tl

    if (return_exdf) {
        # Make a new exdf object from the calculated variables and make sure units
        # are included
        output <- exdf(data.frame(
            alpha_g = alpha_g,
            Gamma_star = Gamma_star,
            J_tl = J_tl,
            Rd_tl = Rd_tl,
            Tp = Tp,
            Vcmax_tl = Vcmax_tl,
            Ac = Ac,
            Aj = Aj,
            Ap = Ap,
            An = An,
            Wc = Wc,
            Wj = Wj,
            Wp = Wp,
            Vc = Wcjp,
            c3_assimilation_msg = msg
        ))

        document_variables(
            output,
            c('calculate_c3_assimilation', 'alpha_g',             'dimensionless'),
            c('calculate_c3_assimilation', 'Gamma_star',          'micromol mol^(-1)'),
            c('calculate_c3_assimilation', 'J_tl',                'micromol m^(-2) s^(-1)'),
            c('calculate_c3_assimilation', 'Rd_tl',               'micromol m^(-2) s^(-1)'),
            c('calculate_c3_assimilation', 'Tp',                  'micromol m^(-2) s^(-1)'),
            c('calculate_c3_assimilation', 'Vcmax_tl',            'micromol m^(-2) s^(-1)'),
            c('calculate_c3_assimilation', 'Ac',                  'micromol m^(-2) s^(-1)'),
            c('calculate_c3_assimilation', 'Aj',                  'micromol m^(-2) s^(-1)'),
            c('calculate_c3_assimilation', 'Ap',                  'micromol m^(-2) s^(-1)'),
            c('calculate_c3_assimilation', 'An',                  'micromol m^(-2) s^(-1)'),
            c('calculate_c3_assimilation', 'Wc',                  'micromol m^(-2) s^(-1)'),
            c('calculate_c3_assimilation', 'Wj',                  'micromol m^(-2) s^(-1)'),
            c('calculate_c3_assimilation', 'Wp',                  'micromol m^(-2) s^(-1)'),
            c('calculate_c3_assimilation', 'Vc',                  'micromol m^(-2) s^(-1)'),
            c('calculate_c3_assimilation', 'c3_assimilation_msg', '')
        )
    } else {
        return(list(An = An, Ac = Ac, Aj = Aj, Ap = Ap, Wc = Wc, Wj = Wj))
    }
}
