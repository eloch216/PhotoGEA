# Here we implement equations defined in: Grassi, G. & Magnani, F. "Stomatal,
# mesophyll conductance and biochemical limitations to photosynthesis as
# affected by drought and leaf ontogeny in ash and oak trees." Plant, Cell &
# Environment 28, 834–849 (2005).

calculate_c3_limitations <- function(
    exdf_obj,
    POc = 210000, # microbar (typically this value is known from the experimental setup)
    atp_use = 4.0,
    nadph_use = 8.0,
    cc_column_name = 'Cc',
    gamma_star_column_name = 'Gamma_star',
    gmc_column_name = 'gmc',
    gsc_column_name = 'gsc',
    kc_column_name = 'Kc',
    ko_column_name = 'Ko',
    total_pressure_column_name = 'total_pressure',
    vcmax_column_name = 'Vcmax_tl',
    j_column_name = 'J_tl'
)
{
    # Check inputs
    if (!is.exdf(exdf_obj)) {
        stop('calculate_c3_limitations requires an exdf object')
    }

    # Make sure the required variables are defined and have the correct units
    required_variables <- list()
    required_variables[[cc_column_name]]             <- 'micromol mol^(-1)'
    required_variables[[gamma_star_column_name]]     <- 'micromol mol^(-1)'
    required_variables[[gmc_column_name]]            <- 'mol m^(-2) s^(-1) bar^(-1)'
    required_variables[[gsc_column_name]]            <- 'mol m^(-2) s^(-1)'
    required_variables[[j_column_name]]              <- 'micromol m^(-2) s^(-1)'
    required_variables[[kc_column_name]]             <- 'micromol mol^(-1)'
    required_variables[[ko_column_name]]             <- 'mmol mol^(-1)'
    required_variables[[total_pressure_column_name]] <- 'bar'
    required_variables[[vcmax_column_name]]          <- 'micromol m^(-2) s^(-1)'

    check_required_variables(exdf_obj, required_variables)

    # Extract key variables to make the following equations simpler. Note that
    # we convert the units for some of these.
    Cc         <- exdf_obj[, cc_column_name] * 1e-6                                    # dimensionless from mol / mol
    Gamma_star <- exdf_obj[, gamma_star_column_name] * 1e-6                            # dimensionless from mol / mol
    gmc        <- exdf_obj[, gmc_column_name] * exdf_obj[, total_pressure_column_name] # mol / m^2 / s
    gsc        <- exdf_obj[, gsc_column_name]                                          # mol / m^2 / s
    J          <- exdf_obj[, j_column_name] * 1e-6                                     # mol / m^2 / s
    Kc         <- exdf_obj[, kc_column_name] * 1e-6                                    # dimensionless from mol / mol
    Ko         <- exdf_obj[, ko_column_name] * 1e3                                     # dimensionless from mol / mol
    O          <- POc * 1e-6 / exdf_obj[, total_pressure_column_name]                  # dimensionless from mol / mol
    Vcmax      <- exdf_obj[, vcmax_column_name] * 1e-6                                 # mol / m^2 / s

    # Total conductance
    gtot <- 1 / (1 / gmc + 1 / gsc) # mol / m^2 / s

    # Partial derivative of A with respect to Cc, assuming Rubisco-limited
    # assimilation
    Km <- Kc * (1 + O / Ko)                                 # dimensionless
    dAdC_rubisco <- Vcmax * (Gamma_star + Km) / (Cc + Km)^2 # mol / m^2 / s

    # Partial derivative of A with respect to Cc, assuming
    # RuBP-regeneration-limited assimilation
    dAdC_j <- J * Gamma_star * (atp_use + nadph_use) /
        (atp_use * Cc + nadph_use * Gamma_star)^2 # mol / m^2 / s

    # Calculate limitations and include them in the exdf object, along with
    # partial derivatives
    exdf_obj[, 'dAdC_rubisco'] <- dAdC_rubisco
    exdf_obj[, 'dAdC_j']       <- dAdC_j
    exdf_obj[, 'ls_rubisco']   <- c3_limitation(gtot, dAdC_rubisco, gsc)
    exdf_obj[, 'lm_rubisco']   <- c3_limitation(gtot, dAdC_rubisco, gmc)
    exdf_obj[, 'lb_rubisco']   <- c3_limitation(gtot, dAdC_rubisco, dAdC_rubisco)
    exdf_obj[, 'ls_j']         <- c3_limitation(gtot, dAdC_j,       gsc)
    exdf_obj[, 'lm_j']         <- c3_limitation(gtot, dAdC_j,       gmc)
    exdf_obj[, 'lb_j']         <- c3_limitation(gtot, dAdC_j,       dAdC_j)

    # Document the columns that were just added and return the exdf
    document_variables(
        exdf_obj,
        c('calculate_c3_limitations', 'dAdC_rubisco', 'mol m^(-2) s^(-1)'),
        c('calculate_c3_limitations', 'dAdC_j',       'mol m^(-2) s^(-1)'),
        c('calculate_c3_limitations', 'ls_rubisco',   'dimensionless'),
        c('calculate_c3_limitations', 'lm_rubisco',   'dimensionless'),
        c('calculate_c3_limitations', 'lb_rubisco',   'dimensionless'),
        c('calculate_c3_limitations', 'ls_j',         'dimensionless'),
        c('calculate_c3_limitations', 'lm_j',         'dimensionless'),
        c('calculate_c3_limitations', 'lb_j',         'dimensionless')
    )
}

# All three parts of Equation 7 can be reproduced with one functional form
c3_limitation <- function(
    g_total,       # mol / m^2 / s
    partial_deriv, # mol / m^2 / s
    g_other        # mol / m^2 / s
)
{
    (g_total / g_other) * partial_deriv / (g_total + partial_deriv) # dimensionless
}
