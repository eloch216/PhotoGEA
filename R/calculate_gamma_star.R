calculate_gamma_star <- function(
    exdf_obj,
    alpha_pr = 0.5,
    oxygen_column_name = 'oxygen',
    specificity_at_tleaf_column_name = 'specificity_at_tleaf',
    tleaf_column_name = 'TleafCnd'
)
{
    # Check inputs
    if (!is.exdf(exdf_obj)) {
        stop('exdf_obj must be an exdf object')
    }

    # Make sure the required variables are defined and have the correct units
    required_variables <- list()
    required_variables[[oxygen_column_name]] <- 'percent'
    required_variables[[specificity_at_tleaf_column_name]] <- 'M / M'
    required_variables[[tleaf_column_name]]  <- 'degrees C'

    check_required_variables(exdf_obj, required_variables)

    # Extract some important columns
    oxygen <- exdf_obj[, oxygen_column_name]                             # percent
    specificity_at_tleaf <- exdf_obj[, specificity_at_tleaf_column_name] # M / M
    Tleaf  <- exdf_obj[, tleaf_column_name]                              # degrees C

    # Convert temperature to Kelvin
    T <- Tleaf + 273.15

    # Equation 18 from Tromans (1998): Henry's constant for O2 in H2O, expressed
    # in mol O2 / kg H2O / atm
    H_O2 <- exp((0.046 * T^2 + 203.357 * T * log(T / 298) - (299.378 + 0.092 * T) * (T - 298) - 20.591e3) / (8.3144 * T))

    # Convert to mol O2 / kg H2O / MPa using 1 atm = 101325 Pa = 0.101325 MPa.
    # Note: this is equivalent to micromol O2 / kg H2O / Pa.
    H_O2 <- H_O2 / 0.101325

    # Equation 4 from Carroll et al. (1991): Henry's constant for CO2 in H2O,
    # expressed in mol CO2 / mol H2O / MPa.
    H_CO2 <- 1 / exp(-6.8346 + 1.2817e4 / T - 3.7668e6 / T^2 + 2.997e8 / T^3)

    # Convert to mol CO2 / kg H2O / MPa using 1 mol H2O = 0.01801528 kg.
    # Note: this is equivalent to micromol CO2 / kg H2O / Pa.
    H_CO2 <- H_CO2 / 0.01801528

    # Get the Rubisco specificity on a gas concentration basis (Pa / Pa)
    specificity_gas_basis <-
        specificity_at_tleaf * H_CO2 / H_O2

    # Convert oxygen percentage to a concentration in micromol / mol
    O2 <- (oxygen * 1e-2) * 1e6

    # Get Gamma_star in micromol / mol
    Gamma_star <- alpha_pr * O2 / specificity_gas_basis

    # Store the new variables in the exdf object
    exdf_obj[, 'Gamma_star'] <- Gamma_star
    exdf_obj[, 'H_CO2'] <- H_CO2
    exdf_obj[, 'H_O2'] <- H_O2
    exdf_obj[, 'specificity_gas_basis'] <- specificity_gas_basis

    # Document the columns that were added and return the exdf
    document_variables(
        exdf_obj,
        c('calculate_gamma_star', 'Gamma_star',            'micromol mol^(-1)'),
        c('calculate_gamma_star', 'H_CO2',                 'micromol kg^(-1) Pa^(-1)'),
        c('calculate_gamma_star', 'H_O2',                  'micromol kg^(-1) Pa^(-1)'),
        c('calculate_gamma_star', 'specificity_gas_basis', 'Pa / Pa')
    )
}
