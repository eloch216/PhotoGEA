# Define a few useful constants for the Arrhenius function
ideal_gas_constant <- 8.3145e-3  # kJ / mol / k
absolute_zero <- -273.15         # degrees C
f <- ideal_gas_constant * (25 - absolute_zero)
c_pa_to_ppm <- log(1e6 / 101325)

# Define a helping function that calculates one Arrhenius exponent
arrhenius <- function(
    scaling,           # dimensionless
    activation_energy, # kJ / mol
    temperature        # degrees C
)
{
    temperature_k <- temperature - absolute_zero # Kelvin

    return(exp(scaling - activation_energy / (ideal_gas_constant * temperature_k)))
}

calculate_arrhenius <- function(
    exdf_obj,
    arrhenius_parameters,
    tleaf_column_name = 'TleafCnd'
)
{
    if (!is.exdf(exdf_obj)) {
        stop("calculate_arrhenius requires an exdf object")
    }

    # Make sure the required variables are defined and have the correct units
    required_variables <- list()
    required_variables[[tleaf_column_name]] <- "degrees C"

    check_required_variables(exdf_obj, required_variables)

    # Get some important information to simplify the following commands
    leaf_temperature <- exdf_obj[, tleaf_column_name]
    param_names <- names(arrhenius_parameters)

    # Calculate the Arrhenius exponents and store them in the exdf object
    for (i in seq_along(arrhenius_parameters)) {
        param <- arrhenius_parameters[[i]]

        # Make sure the necessary elements are defined
        if (!all(c('c', 'Ea', 'units') %in% names(param))) {
            stop(paste0(
                'Arrhenius parameter named `', param_names[i],
                '` has the following elements: ', paste(names(param), collapse = ', '),
                '; elements named `c`, `Ea`, and `units` are required.'
            ))
        }

        # Calculate the exponents
        exdf_obj <- set_variable(
            exdf_obj,
            param_names[i],
            param$units,
            'calculate_arrhenius',
            arrhenius(param$c, param$Ea, leaf_temperature)
        )
    }

    return(exdf_obj)
}
