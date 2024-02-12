# Here we implement equations defined in: Warren et al. "Transfer conductance in
# second growth Douglas-fir (Pseudotsuga menziesii (Mirb.)Franco) canopies."
# Plant, Cell & Environment 26, 1215–1227 (2003).

calculate_c3_limitations_warren <- function(
    exdf_obj,
    POc = 210000,
    atp_use = 4.0,
    nadph_use = 8.0,
    curvature_cj = 1.0,
    curvature_cjp = 1.0,
    ca_column_name = 'Ca',
    cc_column_name = 'Cc',
    ci_column_name = 'Ci',
    j_norm_column_name = 'J_norm',
    kc_column_name = 'Kc',
    ko_column_name = 'Ko',
    rd_norm_column_name = 'Rd_norm',
    total_pressure_column_name = 'total_pressure',
    vcmax_norm_column_name = 'Vcmax_norm'
)
{
    # Check inputs
    if (!is.exdf(exdf_obj)) {
        stop('calculate_c3_limitations_warren requires an exdf object')
    }

    # Make sure the required variables are defined and have the correct units.
    required_variables <- list()
    required_variables[['alpha_g']]                  <- unit_dictionary[['alpha_g']]
    required_variables[['Gamma_star']]               <- unit_dictionary[['Gamma_star']]
    required_variables[['J_at_25']]                  <- unit_dictionary[['J_at_25']]
    required_variables[['Rd_at_25']]                 <- unit_dictionary[['Rd_at_25']]
    required_variables[['Tp']]                       <- unit_dictionary[['Tp']]
    required_variables[['Vcmax_at_25']]              <- unit_dictionary[['Vcmax_at_25']]
    required_variables[[ca_column_name]]             <- unit_dictionary[['Ca']]
    required_variables[[cc_column_name]]             <- unit_dictionary[['Cc']]
    required_variables[[ci_column_name]]             <- unit_dictionary[['Ci']]
    required_variables[[j_norm_column_name]]         <- unit_dictionary[['J_norm']]
    required_variables[[kc_column_name]]             <- unit_dictionary[['Kc']]
    required_variables[[ko_column_name]]             <- unit_dictionary[['Ko']]
    required_variables[[rd_norm_column_name]]        <- unit_dictionary[['Rd_norm']]
    required_variables[[total_pressure_column_name]] <- unit_dictionary[['total_pressure']]
    required_variables[[vcmax_norm_column_name]]     <- unit_dictionary[['Vcmax_norm']]

    check_required_variables(exdf_obj, required_variables)

    # Extract key variables to make the following equations simpler
    Ca <- exdf_obj[, ca_column_name]    # micromol / mol
    Cc <- exdf_obj[, cc_column_name]    # micromol / mol
    Ci <- exdf_obj[, ci_column_name]    # micromol / mol

    # If gsc is as measured and gmc is infinite, we have the measured drawdown
    # across the stomata, but no drawdown across the mesophyll:
    #  Cc_inf_gmc = Ca - (Ca - Ci) - 0 = Ci
    exdf_obj[, 'Cc_inf_gmc'] <- Ci # micromol / mol

    # If gsc is infinite and gmc is as measured, we have no drawdown across the
    # stomata, but the measured drawdown across the mesophyll:
    #  Cc_inf_gsc = Ca - 0 - (Ci - Cc) = Ca - Ci + Cc
    exdf_obj[, 'Cc_inf_gsc'] <- Ca - Ci + Cc # micromol / mol

    # Document the columns that were just added
    exdf_obj <- document_variables(
        exdf_obj,
        c('calculate_c3_limitations_warren', 'Cc_inf_gmc', 'micromol mol^(-1)'),
        c('calculate_c3_limitations_warren', 'Cc_inf_gsc', 'micromol mol^(-1)')
    )

    # Make a helping function for calculating assimilation rates for different
    # Cc column names
    an_from_cc <- function(cc_name) {
        calculate_c3_assimilation(
            exdf_obj,
            '', # alpha_g
            '', # Gamma_star
            '', # J_at_25
            '', # Rd_at_25
            '', # Tp
            '', # Vcmax_at_25
            POc = POc,
            atp_use = atp_use,
            nadph_use = nadph_use,
            curvature_cj = curvature_cj,
            curvature_cjp = curvature_cjp,
            cc_column_name = cc_name,
            j_norm_column_name = j_norm_column_name,
            kc_column_name = kc_column_name,
            ko_column_name = ko_column_name,
            rd_norm_column_name = rd_norm_column_name,
            total_pressure_column_name = total_pressure_column_name,
            vcmax_norm_column_name = vcmax_norm_column_name,
            perform_checks = FALSE,
            return_exdf = TRUE
        )[, 'An']
    }

    # Calculate the net assimilation rate assuming gmc and gsc are as measured
    An <- an_from_cc(cc_column_name) # micromol / m^2 / s

    # Calculate the net assimilation rate assuming gmc is infinite and gsc is as
    # measured
    exdf_obj[, 'An_inf_gmc'] <- an_from_cc('Cc_inf_gmc') # micromol / m^2 / s

    # Calculate the net assimilation rate assuming gmc is as measured and gsc is
    # infinite
    exdf_obj[, 'An_inf_gsc'] <- an_from_cc('Cc_inf_gsc') # micromol / m^2 / s

    # Calculate the limitations using Equations 10 and 11
    exdf_obj[, 'lm_warren'] <- c3_limitation_warren(exdf_obj[, 'An_inf_gmc'], An) # dimensionless
    exdf_obj[, 'ls_warren'] <- c3_limitation_warren(exdf_obj[, 'An_inf_gsc'], An) # dimensionless

    # Document the columns that were just added and return the exdf
    document_variables(
        exdf_obj,
        c('calculate_c3_limitations_warren', 'An_inf_gmc', 'micromol m^(-2) s^(-1)'),
        c('calculate_c3_limitations_warren', 'An_inf_gsc', 'micromol m^(-2) s^(-1)'),
        c('calculate_c3_limitations_warren', 'lm_warren',  'dimensionless'),
        c('calculate_c3_limitations_warren', 'ls_warren',  'dimensionless')
    )
}

# Compares assimilation rates assuming an infinite conductance against a base
# assimilation rate
c3_limitation_warren <- function(An_inf, An_base) {
    (An_inf - An_base) / An_inf
}
