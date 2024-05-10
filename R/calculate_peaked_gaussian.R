# Define a few useful constants for the Arrhenius function
ideal_gas_constant <- 8.3145e-3  # kJ / mol / k
absolute_zero <- -273.15         # degrees C
f <- ideal_gas_constant * (25 - absolute_zero)
c_pa_to_ppm <- log(1e6 / 101325)

# Define a helping function that calculates one peaked Gaussian factor
peaked_gaussian <- function(
    t_opt, # degrees C
    sigma, # degrees C
    t_leaf # degrees C
)
{
    exp(-(t_leaf - t_opt)^2 / sigma^2) # dimensionless
}

calculate_peaked_gaussian <- function(
    exdf_obj,
    peaked_gaussian_parameters,
    tleaf_column_name = 'TleafCnd'
)
{
    if (!is.exdf(exdf_obj)) {
        stop("calculate_peaked_gaussian requires an exdf object")
    }

    # Make sure the required variables are defined and have the correct units
    required_variables <- list()
    required_variables[[tleaf_column_name]] <- "degrees C"

    check_required_variables(exdf_obj, required_variables)

    # Get some important information to simplify the following commands
    leaf_temperature <- exdf_obj[, tleaf_column_name]
    param_names <- names(peaked_gaussian_parameters)

    # Calculate the peaked Gaussian factors and store them in the exdf object
    for (i in seq_along(peaked_gaussian_parameters)) {
        param <- peaked_gaussian_parameters[[i]]

        # Make sure the necessary elements are defined
        if (!all(c('t_opt', 'sigma', 'units') %in% names(param))) {
            stop(paste0(
                'Peaked Gaussian parameter named `', param_names[i],
                '` has the following elements: ', paste(names(param), collapse = ', '),
                '; elements named `t_opt`, `sigma`, and `units` are required.'
            ))
        }

        # Calculate the exponents
        exdf_obj <- set_variable(
            exdf_obj,
            param_names[i],
            param$units,
            'calculate_peaked_gaussian',
            peaked_gaussian(param$t_opt, param$sigma, leaf_temperature)
        )
    }

    return(exdf_obj)
}
