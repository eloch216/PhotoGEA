calculate_ternary_correction <- function(
    exdf_obj,
    ci_column_name = 'Ci',
    co2_s_column_name = 'CO2_s',
    csurface_column_name = 'Csurface',
    e_column_name = 'E',
    gtc_column_name = 'gtc'
)
{
    # Check the inputs
    if (!is.exdf(exdf_obj)) {
        stop('exdf_obj must be an exdf object')
    }

    # Make sure the required variables are defined and have the correct units
    required_variables <- list()
    required_variables[[ci_column_name]]       <- 'micromol mol^(-1)'
    required_variables[[co2_s_column_name]]    <- 'micromol mol^(-1)'
    required_variables[[csurface_column_name]] <- 'micromol mol^(-1)'
    required_variables[[e_column_name]]        <- 'mol m^(-2) s^(-1)'
    required_variables[[gtc_column_name]]      <- 'mol m^(-2) s^(-1)'

    check_required_variables(exdf_obj, required_variables)

    # Extract some important columns
    Ca  <- exdf_obj[, co2_s_column_name]    # micromol / mol
    Ci  <- exdf_obj[, ci_column_name]       # micromol / mol
    Cs  <- exdf_obj[, csurface_column_name] # micromol / mol
    E   <- exdf_obj[, e_column_name]        # mol / m^2 / s
    gac <- exdf_obj[, gtc_column_name]      # mol / m^2 / s

    # Define some constants
    a_b <- 2.9 # ppt; isotopic fractionation during diffusion through the laminar boundary layer
    a_s <- 4.4 # ppt; isotopic fractionation during diffusion through air

    # The weighted isotopic fractionation across the boundary layer and stomata
    a_bar <- (a_b * (Ca - Cs) + a_s * (Cs - Ci)) / (Ca - Ci) # ppt

    alpha_ac <- 1 + a_bar * 1e-3  # dimensionless; fractionation factor during diffusion through air

    # Ternary correction factor
    t <- alpha_ac * E / (2 * gac) # dimensionless

    # Store the calculated quantities in the exdf object
    exdf_obj[, 'a_bar'] <- a_bar
    exdf_obj[, 'alpha_ac'] <- alpha_ac
    exdf_obj[, 't'] <- t

    # # Document the columns that were added and return the exdf
    document_variables(
        exdf_obj,
        c('calculate_gm_busch', 'a_bar',    'ppt'),
        c('calculate_gm_busch', 'alpha_ac', 'dimensionless'),
        c('calculate_gm_busch', 't',        'dimensionless')
    )
}
