# Define a helping function that calculates one Johnson, Eyring, & Williams
# factor. Here we use the `arrhenius` function (defined in
# `calculate_temperature_response_arrhenius.R`) to simplify the calculations.
# The `ideal_gas_constant` variable is also defined in that file.
johnson_eyring_williams <- function(
    scaling,               # dimensionless
    activation_enthalpy,   # kJ / mol
    deactivation_enthalpy, # kJ / mol
    entropy,               # kJ / K / mol
    temperature            # degrees C
)
{
    top <- arrhenius(scaling, activation_enthalpy, temperature)

    bot <- 1.0 + arrhenius(entropy / ideal_gas_constant, deactivation_enthalpy, temperature)

    return(top / bot)
}

calculate_temperature_response_johnson <- function(
    exdf_obj,
    johnson_parameters,
    tleaf_column_name = 'TleafCnd'
)
{
    if (!is.exdf(exdf_obj)) {
        stop('calculate_temperature_response_johnson requires an exdf object')
    }

    # Make sure the required variables are defined and have the correct units
    required_variables <- list()
    required_variables[[tleaf_column_name]] <- 'degrees C'

    check_required_variables(exdf_obj, required_variables)

    # Get some important information to simplify the following commands
    leaf_temperature <- exdf_obj[, tleaf_column_name]
    param_names <- names(johnson_parameters)

    # Calculate the parameter values and store them in the exdf object
    for (i in seq_along(johnson_parameters)) {
        param <- johnson_parameters[[i]]

        # Make sure the necessary elements are defined
        if (!all(c('c', 'Ha', 'Hd', 'S', 'units') %in% names(param))) {
            stop(paste0(
                'Johnson parameter named `', param_names[i],
                '` has the following elements: ', paste(names(param), collapse = ', '),
                '; elements named `c`, `Ha`, `Hd`, `S`, and `units` are required.'
            ))
        }

        # Calculate the exponents
        exdf_obj <- set_variable(
            exdf_obj,
            param_names[i],
            param$units,
            'calculate_temperature_response_johnson',
            johnson_eyring_williams(param$c, param$Ha, param$Hd, param$S, leaf_temperature)
        )
    }

    return(exdf_obj)
}
