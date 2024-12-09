# Calculate a trust value indicator for a parameter estimate using the following
# rules:
#
# - Level 0: The corresponding process is never limiting at any point in the
#   curve. This parameter estimate cannot be trusted.
#
# - Level 1: The corresponding process is limiting at more than one point in the
#   curve, but the upper confidence limit for the parameter value is infinity.
#   This parameter estimate is likely unreliable.
#
# - Level 2: The corresponding process is limiting at more than one point in the
#   curve and there is a finite upper confidence limit for the parameter value.
#   This parameter estimate is trustworthy.
trust_value <- function(
    unreliable_npts, # TRUE means that the corresponding process is never limiting at any point in the curve
    unreliable_inf   # TRUE means that the upper confidence limit for the parameter is infinity
)
{
    if (unreliable_npts) {
        0
    } else if (unreliable_inf) {
        1
    } else {
        2
    }
}

# Decide whether to remove the estimated parameter value using the following
# rules:
#
# - When remove_unreliable_param is 0, all parameter estimates should be kept.
#
# - When remove_unreliable_param is 1, severely unreliable parameter estimates
#   (where trust is 0) should be removed and other parameter estimates should be
#   kept.
#
# - When remove_unreliable_param is 2, any potentially unreliable parameter
#   estimates (where trust is 0 or 1) should be removed.
remove_estimate <- function(trust, remove_unreliable_param) {
    if (remove_unreliable_param == 0) {
        FALSE
    } else if (remove_unreliable_param == 1) {
        if (trust == 0) {
            TRUE
        } else {
            FALSE
        }
    } else {
        if (trust == 0 || trust == 1) {
            TRUE
        } else {
            FALSE
        }
    }
}

# Make sure the user input is acceptable
check_param_setting <- function(remove_unreliable_param) {
    if (!remove_unreliable_param %in% c(0, 1, 2)) {
        stop('`remove_unreliable_param` must be 0, 1, or 2')
    }
}
