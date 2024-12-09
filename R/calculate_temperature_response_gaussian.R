# Define a helping function that calculates one Gaussian factor
gaussian_response <- function(
    optimum_rate, # any units
    t_opt,        # degrees C
    sigma,        # degrees C
    t_leaf        # degrees C
)
{
    optimum_rate * exp(-(t_leaf - t_opt)^2 / sigma^2) # same units as optimum_rate
}

calculate_temperature_response_gaussian <- function(
    exdf_obj,
    gaussian_parameters,
    tleaf_column_name = 'TleafCnd'
)
{
    if (!is.exdf(exdf_obj)) {
        stop("calculate_temperature_response_gaussian requires an exdf object")
    }

    # Make sure the required variables are defined and have the correct units
    required_variables <- list()
    required_variables[[tleaf_column_name]] <- "degrees C"

    check_required_variables(exdf_obj, required_variables)

    # Get some important information to simplify the following commands
    leaf_temperature <- exdf_obj[, tleaf_column_name]
    param_names <- names(gaussian_parameters)

    # Calculate the Gaussian factors and store them in the exdf object
    for (i in seq_along(gaussian_parameters)) {
        param <- gaussian_parameters[[i]]

        # Make sure the necessary elements are defined
        if (!all(c('optimum_rate', 't_opt', 'sigma', 'units') %in% names(param))) {
            stop(paste0(
                'Gaussian parameter named `', param_names[i],
                '` has the following elements: ', paste(names(param), collapse = ', '),
                '; elements named `optimum_rate`, `t_opt`, `sigma`, and `units` are required.'
            ))
        }

        # Calculate the exponents
        exdf_obj <- set_variable(
            exdf_obj,
            param_names[i],
            param$units,
            'calculate_temperature_response_gaussian',
            gaussian_response(param$optimum_rate, param$t_opt, param$sigma, leaf_temperature)
        )
    }

    return(exdf_obj)
}
