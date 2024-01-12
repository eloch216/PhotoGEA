# Here we implement equations defined in: Warren et al. "Transfer conductance in
# second growth Douglas-fir (Pseudotsuga menziesii (Mirb.)Franco) canopies."
# Plant, Cell & Environment 26, 1215–1227 (2003).

calculate_c3_limitations_warren <- function(
    exdf_obj,
    POc = 210000,
    atp_use = 4.0,
    nadph_use = 8.0,
    alpha = 0.0,
    curvature_cj = 1.0,
    curvature_cjp = 1.0,
    ca_column_name = 'Ca',
    cc_column_name = 'Cc',
    ci_column_name = 'Ci',
    j_column_name = 'J_at_25',
    j_norm_column_name = 'J_norm',
    kc_column_name = 'Kc',
    ko_column_name = 'Ko',
    rd_column_name = 'Rd_at_25',
    rd_norm_column_name = 'Rd_norm',
    total_pressure_column_name = 'total_pressure',
    tpu_column_name = 'TPU',
    vcmax_column_name = 'Vcmax_at_25',
    vcmax_norm_column_name = 'Vcmax_norm'
)
{
    # Check inputs
    if (!is.exdf(exdf_obj)) {
        stop('calculate_c3_limitations_warren requires an exdf object')
    }

    # Make sure the required variables are defined and have the correct units.
    # We don't need to check all of them, since calculate_c3_assimilation will
    # perform checks.
    required_variables <- list()
    required_variables[['Gamma_star']]      <- unit_dictionary[['Gamma_star']]
    required_variables[[ca_column_name]]    <- 'micromol mol^(-1)'
    required_variables[[cc_column_name]]    <- 'micromol mol^(-1)'
    required_variables[[ci_column_name]]    <- 'micromol mol^(-1)'
    required_variables[[j_column_name]]     <- 'micromol m^(-2) s^(-1)'
    required_variables[[rd_column_name]]    <- 'micromol m^(-2) s^(-1)'
    required_variables[[tpu_column_name]]   <- 'micromol m^(-2) s^(-1)'
    required_variables[[vcmax_column_name]] <- 'micromol m^(-2) s^(-1)'

    check_required_variables(exdf_obj, required_variables)

    # Extract key variables to make the following equations simpler
    Ca         <- exdf_obj[, ca_column_name]    # micromol / mol
    Cc         <- exdf_obj[, cc_column_name]    # micromol / mol
    Ci         <- exdf_obj[, ci_column_name]    # micromol / mol
    Gamma_star <- exdf_obj[, 'Gamma_star']      # micromol / mol
    J          <- exdf_obj[, j_column_name]     # micromol / m^2 / s
    Rd         <- exdf_obj[, rd_column_name]    # micromol / m^2 / s
    TPU        <- exdf_obj[, tpu_column_name]   # micromol / m^2 / s
    Vcmax      <- exdf_obj[, vcmax_column_name] # micromol / m^2 / s

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
        sapply(seq_len(nrow(exdf_obj)), function(i) {
            calculate_c3_assimilation(
                exdf_obj[i, , TRUE],
                J[i],
                Gamma_star[i],
                Rd[i],
                TPU[i],
                Vcmax[i],
                POc = POc,
                atp_use = atp_use,
                nadph_use = nadph_use,
                alpha = alpha,
                curvature_cj = curvature_cj,
                curvature_cjp = curvature_cjp,
                cc_column_name = cc_name,
                j_norm_column_name = j_norm_column_name,
                kc_column_name = kc_column_name,
                ko_column_name = ko_column_name,
                rd_norm_column_name = rd_norm_column_name,
                total_pressure_column_name = total_pressure_column_name,
                vcmax_norm_column_name = vcmax_norm_column_name,
                perform_checks = TRUE,
                return_exdf = FALSE
            )[['An']]
        })
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
