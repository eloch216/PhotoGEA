# Define trust indicator values
trust_indicators <- list(
    r       = 'reliable',
    u_inf   = 'unreliable (infinite upper limit)',
    u_never = 'unreliable (process never limiting)',
    u_over  = 'unreliable (insufficient DOF)'
)

# Determine a trust indicator for a parameter estimate
trust_value <- function(
    unreliable_npts, # TRUE means that the corresponding process is never limiting at any point in the curve
    unreliable_inf,  # TRUE means that the upper confidence limit for the parameter is infinity
    dof              # non-positive value means that the fit was overparameterized
)
{
    if (dof < 1.0) {
        trust_indicators[['u_over']]
    } else if (unreliable_npts) {
        trust_indicators[['u_never']]
    } else if (unreliable_inf) {
        trust_indicators[['u_inf']]
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
        c(
            trust_indicators[['u_never']],
            trust_indicators[['u_over']]
        )
    } else if (identical(remove_unreliable_param, 2)) {
        c(
            trust_indicators[['u_inf']],
            trust_indicators[['u_never']],
            trust_indicators[['u_over']]
        )
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
