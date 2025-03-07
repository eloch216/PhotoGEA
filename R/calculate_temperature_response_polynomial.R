# Define a helping function that evaluates one polynomial
evaluate_polynomial <- function(x, coefs) {
    sapply(x, function(val) {
        polynomial_terms <- sapply(seq_along(coefs), function(i) {
            term_order <- i - 1
            coefs[i] * val^term_order
        })
        sum(polynomial_terms)
    })
}

calculate_temperature_response_polynomial <- function(
    exdf_obj,
    polynomial_parameters,
    tleaf_column_name = 'TleafCnd'
)
{
    if (!is.exdf(exdf_obj)) {
        stop("calculate_temperature_response_polynomial requires an exdf object")
    }

    # Make sure the required variables are defined and have the correct units
    required_variables <- list()
    required_variables[[tleaf_column_name]] <- "degrees C"

    check_required_variables(exdf_obj, required_variables)

    # Get some important information to simplify the following commands
    leaf_temperature <- exdf_obj[, tleaf_column_name]
    param_names <- names(polynomial_parameters)

    # Calculate the polynomials and store them in the exdf object
    for (i in seq_along(polynomial_parameters)) {
        param <- polynomial_parameters[[i]]

        # Make sure the necessary elements are defined
        if (!all(c('coef', 'units') %in% names(param))) {
            stop(paste0(
                'Polynomial parameter named `', param_names[i],
                '` has the following elements: ', paste(names(param), collapse = ', '),
                '; elements named `coef` and `units` are required.'
            ))
        }

        # Calculate the exponents
        exdf_obj <- set_variable(
            exdf_obj,
            param_names[i],
            param$units,
            'calculate_temperature_response_polynomial',
            evaluate_polynomial(leaf_temperature, param$coef)
        )
    }

    return(exdf_obj)
}
