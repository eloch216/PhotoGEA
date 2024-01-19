# A helping function that finds confidence intervals for all parameters
confidence_interval_all_param <- function(
    error_function,
    best_fit_parameters,
    luf,
    error_threshold_factor,
    category_name
)
{
    # Extract lower, upper, and fit_options
    lower <- luf$lower
    upper <- luf$upper
    fit_options <- luf$fit_options
    param_to_fit <- luf$param_to_fit

    # Get the best-fit error value
    best_error <- best_fit_parameters[, 'optimum_val']

    # Get the best-fit parameters as a numeric vector
    best_fit_vector <- sapply(names(fit_options), function(pn) {
        best_fit_parameters[, pn]
    })
    best_fit_vector <- as.numeric(best_fit_vector)

    # Find intervals for each fit parameter
    for (i in seq_along(param_to_fit)) {
        if (param_to_fit[i]) {
            # Find limits
            lim <- confidence_interval_one_param(
                i,
                error_function,
                param_to_fit,
                best_fit_vector,
                best_error * error_threshold_factor,
                lower[i],
                upper[i]
            )

            # Get the parameter name
            pn <- names(fit_options)[i]
            pn_lower <- paste0(pn, '_lower')
            pn_upper <- paste0(pn, '_upper')

            # Store the results
            best_fit_parameters <- set_variable(
                best_fit_parameters,
                pn_lower,
                best_fit_parameters$units[, pn],
                category_name,
                lim$conf_lower
            )

            best_fit_parameters <- set_variable(
                best_fit_parameters,
                pn_upper,
                best_fit_parameters$units[, pn],
                category_name,
                lim$conf_upper
            )
        }
    }

    best_fit_parameters
}

# A helping function that finds confidence intervals for one parameter
confidence_interval_one_param <- function(
    index,
    error_function,
    param_to_fit,
    best_fit_inputs,
    error_threshold,
    lower_lim,
    upper_lim
)
{
    # Use partial application to create a 1D error function
    erf_i <- function(x) {
        inputs <- best_fit_inputs
        inputs[index] <- x
        inputs <- inputs[param_to_fit]
        error_function(inputs) - error_threshold
    }

    # Find lower limit
    conf_lower <- tryCatch(
        {
            stats::uniroot(
                erf_i,
                c(lower_lim, best_fit_inputs[index]),
                extendInt = 'downX'
            )[['root']]
        },
        error = function(e) {-Inf}
    )

    # Find upper limit
    conf_upper <- tryCatch(
        {
            stats::uniroot(
                erf_i,
                c(best_fit_inputs[index], upper_lim),
                extendInt = 'upX'
            )[['root']]
        },
        error = function(e) {Inf}
    )

    list(
        conf_lower = conf_lower,
        conf_upper = conf_upper
    )
}
