calculate_c3_assimilation <- function(
    exdf_obj,
    J,            # micromol / m^2 / s   (at 25 degrees C; typically this value is being fitted)
    Rd,           # micromol / m^2 / s   (at 25 degrees C; typically this value is being fitted)
    TPU,          # micromol / m^2 / s   (typically this value is being fitted)
    Vcmax,        # micromol / m^2 / s   (at 25 degrees C; typically this value is being fitted)
    POc = 210000, # microbar             (typically this value is known from the experimental setup)
    atp_use = 4.0,
    nadph_use = 8.0,
    alpha = 0.0,
    curvature_cj = 1.0,
    curvature_cjp = 1.0,
    cc_column_name = 'Cc',
    gamma_star_column_name = 'Gamma_star',
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
        required_variables[[gamma_star_column_name]]     <- 'micromol mol^(-1)'
        required_variables[[j_norm_column_name]]         <- 'normalized to J at 25 degrees C'
        required_variables[[kc_column_name]]             <- 'micromol mol^(-1)'
        required_variables[[ko_column_name]]             <- 'mmol mol^(-1)'
        required_variables[[rd_norm_column_name]]        <- 'normalized to Rd at 25 degrees C'
        required_variables[[total_pressure_column_name]] <- 'bar'
        required_variables[[vcmax_norm_column_name]]     <- 'normalized to Vcmax at 25 degrees C'

        check_required_variables(exdf_obj, required_variables)

        # Make sure certain inputs lie on [0,1]
        check_zero_one <- list(
            alpha = alpha,
            curvature_cj = curvature_cj,
            curvature_cjp = curvature_cjp
        )

        sapply(seq_along(check_zero_one), function(i) {
            if (check_zero_one[[i]] < 0 || check_zero_one[[i]] > 1) {
                stop(paste(names(check_zero_one)[i], 'must be >= 0 and <= 1'))
            }
        })

        # Make sure the Cc values are all positive
        if (any(exdf_obj[, cc_column_name] <= 0)) {
            stop('All Cc values must be positive')
        }
    }

    # Extract a few columns from the exdf object to make the equations easier to
    # read, converting units as necessary
    pressure <- exdf_obj[, total_pressure_column_name] # bar

    PCc <- exdf_obj[, cc_column_name] * pressure # microbar

    Kc <- exdf_obj[, kc_column_name] * pressure                 # microbar
    Ko <- exdf_obj[, ko_column_name] * pressure * 1000          # microbar
    Gamma_star <- exdf_obj[, gamma_star_column_name] * pressure # microbar

    Vcmax_tl <- Vcmax * exdf_obj[, vcmax_norm_column_name] # micromol / m^2 / s
    Rd_tl <- Rd * exdf_obj[, rd_norm_column_name]          # micromol / m^2 / s
    J_tl <- J * exdf_obj[, j_norm_column_name]             # micromol / m^2 / s

    # Rubisco-limited carboxylation (micromol / m^2 / s)
    Wc <- PCc * Vcmax_tl / (PCc + Kc * (1.0 + POc / Ko))

    # RuBP-regeneration-limited carboxylation (micromol / m^2 / s)
    Wj <- PCc * J_tl / (atp_use * PCc + nadph_use * Gamma_star)

    # TPU-limited carboxylation (micromol / m^2 / s)
    Wp <- PCc * 3 * TPU / (PCc - Gamma_star * (1 + 3 * alpha))
    Wp[PCc <= Gamma_star * (1 + 3 * alpha)] <- Inf

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
            J_tl = J_tl,
            Rd_tl = Rd_tl,
            TPU = TPU,
            Vcmax_tl = Vcmax_tl,
            Ac = Ac,
            Aj = Aj,
            Ap = Ap,
            An = An,
            Wc = Wc,
            Wj = Wj,
            Wp = Wp,
            Vc = Wcjp
        ))

        document_variables(
            output,
            c('calculate_c3_assimilation', 'J_tl',       'micromol m^(-2) s^(-1)'),
            c('calculate_c3_assimilation', 'Rd_tl',      'micromol m^(-2) s^(-1)'),
            c('calculate_c3_assimilation', 'TPU',        'micromol m^(-2) s^(-1)'),
            c('calculate_c3_assimilation', 'Vcmax_tl',   'micromol m^(-2) s^(-1)'),
            c('calculate_c3_assimilation', 'Ac',         'micromol m^(-2) s^(-1)'),
            c('calculate_c3_assimilation', 'Aj',         'micromol m^(-2) s^(-1)'),
            c('calculate_c3_assimilation', 'Ap',         'micromol m^(-2) s^(-1)'),
            c('calculate_c3_assimilation', 'An',         'micromol m^(-2) s^(-1)'),
            c('calculate_c3_assimilation', 'Wc',         'micromol m^(-2) s^(-1)'),
            c('calculate_c3_assimilation', 'Wj',         'micromol m^(-2) s^(-1)'),
            c('calculate_c3_assimilation', 'Wp',         'micromol m^(-2) s^(-1)'),
            c('calculate_c3_assimilation', 'Vc',         'micromol m^(-2) s^(-1)')
        )
    } else {
        return(list(An = An, Ac = Ac, Aj = Aj, Ap = Ap, Wc = Wc, Wj = Wj))
    }
}
