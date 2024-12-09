calculate_temperature_response <- function(
    exdf_obj,
    temperature_response_parameters,
    tleaf_column_name = 'TleafCnd'
)
{
    if (!is.exdf(exdf_obj)) {
        stop("calculate_temperature_response requires an exdf object")
    }

    # Get some important information to simplify the following commands
    param_names <- names(temperature_response_parameters)

    # Calculate the temperature-dependent values of each parameter
    for (i in seq_along(temperature_response_parameters)) {
        param <- temperature_response_parameters[[i]]

        # Get the type of temperature response function to use for this
        # parameter
        if (!'type' %in% names(param)) {
            stop(paste0(
                'Temperature response parameter set named `', param_names[i],
                '` does not specify a `type` value'
            ))
        }

        type <- tolower(param[['type']])

        # Make sure units are defined
        if (!'units' %in% names(param)) {
            stop(paste0(
                'Temperature response parameter set named `', param_names[i],
                '` does not specify a `units` value'
            ))
        }

        # Get the temperature response function itself
        trf <- switch(type,
            arrhenius = calculate_temperature_response_arrhenius,
            gaussian = calculate_temperature_response_gaussian,
            johnson = calculate_temperature_response_johnson,
            stop(paste0(
                'Temperature response parameter set named `', param_names[i],
                '` specifies an unsupported `type` value: `', param[['type']],
                '`. The available options are: arrhenius, gaussian, and johnson.'
            ))
        )

        # Apply the temperature response function
        param_info <- list(param)
        names(param_info) <- param_names[i]

        exdf_obj <- trf(exdf_obj, param_info, tleaf_column_name)
    }

    return(exdf_obj)
}
