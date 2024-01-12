calculate_c3_variable_j <- function(
    exdf_obj,
    Rd_at_25, # micromol / m^2 / s (at 25 degrees C; typically this value is being fitted)
    tau,      # dimensionless      (typically this value is being fitted)
    atp_use = 4.0,
    nadph_use = 8.0,
    a_column_name = 'A',
    ci_column_name = 'Ci',
    gamma_star_column_name = 'Gamma_star',
    phips2_column_name = 'PhiPS2',
    qin_column_name = 'Qin',
    rd_norm_column_name = 'Rd_norm',
    total_pressure_column_name = 'total_pressure',
    perform_checks = TRUE,
    return_exdf = TRUE
)
{
    # Define flexible parameters
    flexible_param <- list(
        Rd_at_25 = Rd_at_25,
        tau = tau
    )

    if (perform_checks) {
        if (!is.exdf(exdf_obj)) {
            stop('calculate_c3_variable_j requires an exdf object')
        }

        # Make sure the required variables are defined and have the correct units
        required_variables <- list()
        required_variables[[a_column_name]]              <- 'micromol m^(-2) s^(-1)'
        required_variables[[ci_column_name]]             <- 'micromol mol^(-1)'
        required_variables[[gamma_star_column_name]]     <- 'micromol mol^(-1)'
        required_variables[[phips2_column_name]]         <- NA
        required_variables[[qin_column_name]]            <- 'micromol m^(-2) s^(-1)'
        required_variables[[rd_norm_column_name]]        <- 'normalized to Rd at 25 degrees C'
        required_variables[[total_pressure_column_name]] <- 'bar'

        required_variables <-
            require_flexible_param(required_variables, flexible_param)

        check_required_variables(exdf_obj, required_variables)
    }

    # Retrieve values of flexible parameters as necessary
    if (!is.numeric(Rd_at_25)) {Rd_at_25 <- exdf_obj[, 'Rd_at_25']}
    if (!is.numeric(tau))      {tau      <- exdf_obj[, 'tau']}

    # Extract a few columns from the exdf object to make the equations easier to
    # read, converting units as necessary
    pressure   <- exdf_obj[, total_pressure_column_name]        # bar
    Ci         <- exdf_obj[, ci_column_name]                    # micromol / mol
    PCi        <- Ci * pressure                                 # microbar
    Gamma_star <- exdf_obj[, gamma_star_column_name] * pressure # microbar

    An     <- exdf_obj[, a_column_name]                         # micromol / m^2 / s
    PhiPS2 <- exdf_obj[, phips2_column_name]                    # dimensionless
    Qin    <- exdf_obj[, qin_column_name]                       # micromol / m^2 / s

    Rd_tl <- Rd_at_25 * exdf_obj[, rd_norm_column_name]         # micromol / m^2 / s

    # Calculate J_F (actual RuBP regeneration rate as estimated from
    # fluorescence) using Equation 5
    J_F <- tau * Qin * PhiPS2 # micromol / m^2 / s

    # Calculate gmc (mesophyll conductance to CO2) using Equation 6
    AnRd <- An + Rd_tl # micromol / m^2 / s

    gmc_top <- An * (J_F - atp_use * AnRd) # (micromol / m^2 / s)^2

    gmc_bottom <- (Ci - Gamma_star) * J_F -
        (Ci * atp_use + Gamma_star * nadph_use) * AnRd # microbar * (micromol / m^2 / s)

    gmc <- gmc_top / gmc_bottom # mol / m^2 / s / bar

    # Calculate Cc
    Cc <- Ci - An / (gmc * pressure) # micromol / mol

    if (return_exdf) {
        # Make a new exdf object from the calculated variables and make sure units
        # are included
        output <- exdf(data.frame(
            tau = tau,
            Rd_tl = Rd_tl,
            J_F = J_F,
            gmc = gmc,
            Cc = Cc
        ))

        document_variables(
            output,
            c('calculate_c3_variable_j', 'tau',   'dimensionless'),
            c('calculate_c3_variable_j', 'Rd_tl', 'micromol m^(-2) s^(-1)'),
            c('calculate_c3_variable_j', 'J_F',   'micromol m^(-2) s^(-1)'),
            c('calculate_c3_variable_j', 'gmc',   'mol m^(-2) s^(-1) bar^(-1)'),
            c('calculate_c3_variable_j', 'Cc',    'micromol mol^(-1)')
        )
    } else {
        return(list(gmc = gmc, Cc = Cc))
    }
}
