# A helping function that finds confidence intervals for all parameters
confidence_interval_all_param <- function(
    error_function,
    best_fit_parameters,
    luf,
    relative_likelihood_threshold,
    category_name
)
{
    # Extract lower, upper, and fit_options
    lower <- luf$lower
    upper <- luf$upper
    fit_options <- luf$fit_options
    param_to_fit <- luf$param_to_fit

    # Get the best-fit parameters as a numeric vector
    best_fit_vector <- sapply(names(fit_options), function(pn) {
        best_fit_parameters[, pn]
    })
    best_fit_vector <- as.numeric(best_fit_vector)

    # Get the best-fit error value
    best_error <- error_function(best_fit_vector[param_to_fit])

    # Find intervals for each fit parameter
    for (i in seq_along(param_to_fit)) {
        # Get the parameter name
        pn <- names(fit_options)[i]
        pn_lower <- paste0(pn, '_lower')
        pn_upper <- paste0(pn, '_upper')

        # Find the limits
        lim <- if (param_to_fit[i]) {
            # This parameter was fit, so find the confidence limits
            confidence_interval_one_param(
                i,
                error_function,
                param_to_fit,
                best_fit_vector,
                best_error - log(relative_likelihood_threshold),
                lower[i],
                upper[i]
            )
        } else {
            # This parameter was not fit, so just return NA for the limits
            list(
                conf_lower = NA,
                conf_upper = NA
            )
        }

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

    # Find lower limit. If the best-fit value is NA, that means the fit failed
    # and it won't be possible to find a confidence interval. If the best-fit
    # value is at the lower limit, uniroot will complain; in that case, we just
    # use the lower limit.
    conf_lower <- if (is.na(best_fit_inputs[index])) {
        NA
    } else if (lower_lim >= best_fit_inputs[index]) {
        lower_lim
    } else {
        tryCatch(
            {
                stats::uniroot(
                    erf_i,
                    c(lower_lim, best_fit_inputs[index]),
                    extendInt = 'downX'
                )[['root']]
            },
            error = function(e) {-Inf}
        )
    }

    # Find upper limit. If the best-fit value is NA, that means the fit failed
    # and it won't be possible to find a confidence interval. If the best-fit
    # value is at the upper limit, uniroot will compain; in that case, we just
    # use the upper limit.
    conf_upper <- if (is.na(best_fit_inputs[index])) {
        NA
    } else if (best_fit_inputs[index] >= upper_lim) {
        upper_lim
    } else {
        tryCatch(
            {
                stats::uniroot(
                    erf_i,
                    c(best_fit_inputs[index], upper_lim),
                    extendInt = 'upX'
                )[['root']]
            },
            error = function(e) {Inf}
        )
    }

    # Adjust upper limit
    if (!is.na(conf_upper) && conf_upper > 10 * upper_lim) {
        conf_upper <- Inf
    }

    list(
        conf_lower = conf_lower,
        conf_upper = conf_upper
    )
}

# A helping function that calculates confidence intervals at leaf temperature
# from intervals at 25 degrees C
confidence_intervals_leaf_temperature <- function(replicate_identifiers, parameters, fname) {
    for (param in parameters) {
        param_at_25        <- paste0(param, '_at_25')
        param_at_25_lower  <- paste0(param, '_at_25_lower')
        param_at_25_upper  <- paste0(param, '_at_25_upper')
        param_tl_avg       <- paste0(param, '_tl_avg')
        param_tl_avg_lower <- paste0(param, '_tl_avg_lower')
        param_tl_avg_upper <- paste0(param, '_tl_avg_upper')
        
        param_tl_scale <- replicate_identifiers[, param_tl_avg] / replicate_identifiers[, param_at_25]
        
        replicate_identifiers[, param_tl_avg_lower] <- replicate_identifiers[, param_at_25_lower] * param_tl_scale
        replicate_identifiers[, param_tl_avg_upper] <- replicate_identifiers[, param_at_25_upper] * param_tl_scale
        
        replicate_identifiers <- document_variables(
            replicate_identifiers,
            c(fname, param_tl_avg_lower, replicate_identifiers$units[[param_at_25]]),
            c(fname, param_tl_avg_upper, replicate_identifiers$units[[param_at_25]])
        )
    }
    
    replicate_identifiers
}
