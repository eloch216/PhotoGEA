# Define trust indicator values
trust_indicators <- list(
    u0 = 'unreliable (process never limiting)',
    u1 = 'unreliable (infinite upper limit)',
    r  = 'reliable'
)

# Determine a trust indicator for a parameter estimate
trust_value <- function(
    unreliable_npts, # TRUE means that the corresponding process is never limiting at any point in the curve
    unreliable_inf   # TRUE means that the upper confidence limit for the parameter is infinity
)
{
    if (unreliable_npts) {
        trust_indicators[['u0']]
    } else if (unreliable_inf) {
        trust_indicators[['u1']]
    } else {
        trust_indicators[['r']]
    }
}

# Decide whether to remove the estimated parameter value
remove_estimate <- function(trust, param_types_to_remove) {
    if (trust %in% param_types_to_remove) {
        TRUE
    } else {
        FALSE
    }
}

# Check and/or convert the user input
convert_param_setting <- function(remove_unreliable_param) {
    if (identical(remove_unreliable_param, 0)) {
        c()
    } else if (identical(remove_unreliable_param, 1)) {
        trust_indicators[['u0']]
    } else if (identical(remove_unreliable_param, 2)) {
        c(trust_indicators[['u0']], trust_indicators[['u1']])
    } else {
        type_okay <- sapply(remove_unreliable_param, function(param_type) {
            param_type %in% as.character(trust_indicators)
        })

        if (any(!type_okay)) {
            stop(
                'If `remove_unreliable_param` is not 0, 1, or 2, its elements ',
                'must each be one of the following: ',
                paste0('`', as.character(trust_indicators), '`', collapse = ', '),
                call. = FALSE
            )
        }

        remove_unreliable_param
    }
}
