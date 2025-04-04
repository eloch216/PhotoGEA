# Here we implement equations defined in: Grassi, G. & Magnani, F. "Stomatal,
# mesophyll conductance and biochemical limitations to photosynthesis as
# affected by drought and leaf ontogeny in ash and oak trees." Plant, Cell &
# Environment 28, 834–849 (2005).

calculate_c3_limitations_grassi <- function(
    exdf_obj,
    Wj_coef_C = 4.0,
    Wj_coef_Gamma_star = 8.0,
    cc_column_name = 'Cc',
    gamma_star_column_name = 'Gamma_star_tl',
    gmc_column_name = 'gmc_tl',
    gsc_column_name = 'gsc',
    kc_column_name = 'Kc_tl',
    ko_column_name = 'Ko_tl',
    oxygen_column_name = 'oxygen',
    total_pressure_column_name = 'total_pressure',
    vcmax_column_name = 'Vcmax_tl',
    j_column_name = NULL
)
{
    # Check inputs
    if (!is.exdf(exdf_obj)) {
        stop('calculate_c3_limitations_grassi requires an exdf object')
    }

    # Check to see if we should consider RuBP-regeneration-limited assimilation
    use_j <- !is.null(j_column_name)

    # Make sure the required variables are defined and have the correct units
    required_variables <- list()
    required_variables[[cc_column_name]]             <- unit_dictionary('Cc')
    required_variables[[gamma_star_column_name]]     <- unit_dictionary('Gamma_star_at_25')
    required_variables[[gmc_column_name]]            <- unit_dictionary('gmc_at_25')
    required_variables[[gsc_column_name]]            <- unit_dictionary('gsc')
    required_variables[[kc_column_name]]             <- unit_dictionary('Kc_at_25')
    required_variables[[ko_column_name]]             <- unit_dictionary('Ko_at_25')
    required_variables[[oxygen_column_name]]         <- unit_dictionary('oxygen')
    required_variables[[total_pressure_column_name]] <- unit_dictionary('total_pressure')
    required_variables[[vcmax_column_name]]          <- unit_dictionary('Vcmax_at_25')

    if (use_j) {
        required_variables[[j_column_name]] <- unit_dictionary('J_at_25')
    }

    # Don't throw an error if some columns are all NA
    check_required_variables(exdf_obj, required_variables, check_NA = FALSE)

    # Extract key variables to make the following equations simpler. Note that
    # we convert the units for some of these.
    Cc         <- exdf_obj[, cc_column_name] * 1e-6                                    # dimensionless from mol / mol
    Gamma_star <- exdf_obj[, gamma_star_column_name] * 1e-6                            # dimensionless from mol / mol
    gmc        <- exdf_obj[, gmc_column_name] * exdf_obj[, total_pressure_column_name] # mol / m^2 / s
    gsc        <- exdf_obj[, gsc_column_name]                                          # mol / m^2 / s
    Kc         <- exdf_obj[, kc_column_name] * 1e-6                                    # dimensionless from mol / mol
    Ko         <- exdf_obj[, ko_column_name] * 1e3                                     # dimensionless from mol / mol
    O          <- exdf_obj[, oxygen_column_name] * 1e-2                                # dimensionless from mol / mol
    Vcmax      <- exdf_obj[, vcmax_column_name] * 1e-6                                 # mol / m^2 / s

    J <- if (use_j) {
        exdf_obj[, j_column_name] * 1e-6 # mol / m^2 / s
    } else {
        NA
    }

    # Total conductance
    gtot <- 1 / (1 / gmc + 1 / gsc) # mol / m^2 / s

    # Partial derivative of A with respect to Cc, assuming Rubisco-limited
    # assimilation
    Km <- Kc * (1 + O / Ko)                                 # dimensionless
    dAdC_rubisco <- Vcmax * (Gamma_star + Km) / (Cc + Km)^2 # mol / m^2 / s

    # Include partial derivative and limitations in the exdf object
    exdf_obj[, 'dAdC_rubisco']      <- dAdC_rubisco
    exdf_obj[, 'ls_rubisco_grassi'] <- c3_limitation_grassi(gtot, dAdC_rubisco, gsc)
    exdf_obj[, 'lm_rubisco_grassi'] <- c3_limitation_grassi(gtot, dAdC_rubisco, gmc)
    exdf_obj[, 'lb_rubisco_grassi'] <- c3_limitation_grassi(gtot, dAdC_rubisco, dAdC_rubisco)

    # Document the columns that were just added
    exdf_obj <- document_variables(
        exdf_obj,
        c('calculate_c3_limitations_grassi', 'dAdC_rubisco',      'mol m^(-2) s^(-1)'),
        c('calculate_c3_limitations_grassi', 'ls_rubisco_grassi', 'dimensionless'),
        c('calculate_c3_limitations_grassi', 'lm_rubisco_grassi', 'dimensionless'),
        c('calculate_c3_limitations_grassi', 'lb_rubisco_grassi', 'dimensionless')
    )

    if (use_j) {
        # Partial derivative of A with respect to Cc, assuming
        # RuBP-regeneration-limited assimilation
        dAdC_j <- J * Gamma_star * (Wj_coef_C + Wj_coef_Gamma_star) /
            (Wj_coef_C * Cc + Wj_coef_Gamma_star * Gamma_star)^2 # mol / m^2 / s

        # Include partial derivative and limitations in the exdf object
        exdf_obj[, 'dAdC_j']      <- dAdC_j
        exdf_obj[, 'ls_j_grassi'] <- c3_limitation_grassi(gtot, dAdC_j, gsc)
        exdf_obj[, 'lm_j_grassi'] <- c3_limitation_grassi(gtot, dAdC_j, gmc)
        exdf_obj[, 'lb_j_grassi'] <- c3_limitation_grassi(gtot, dAdC_j, dAdC_j)

        # Document the columns that were just added and return the exdf
        exdf_obj <- document_variables(
            exdf_obj,
            c('calculate_c3_limitations_grassi', 'dAdC_j',      'mol m^(-2) s^(-1)'),
            c('calculate_c3_limitations_grassi', 'ls_j_grassi', 'dimensionless'),
            c('calculate_c3_limitations_grassi', 'lm_j_grassi', 'dimensionless'),
            c('calculate_c3_limitations_grassi', 'lb_j_grassi', 'dimensionless')
        )
    }

    # Return the exdf object
    return(exdf_obj)
}

# All three parts of Equation 7 can be reproduced with one functional form
c3_limitation_grassi <- function(
    g_total,       # mol / m^2 / s
    partial_deriv, # mol / m^2 / s
    g_other        # mol / m^2 / s
)
{
    (g_total / g_other) * partial_deriv / (g_total + partial_deriv) # dimensionless
}
