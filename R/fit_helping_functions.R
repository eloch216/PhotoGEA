# A helping function for determining if a flexible parameter has been set
value_set <- function(x) {
    is.numeric(x) && !is.na(x)
}

# A helping function for checking units of "flexible" parameters (uses the unit
# dictionary)
require_flexible_param <- function(required_variables, flexible_param) {
    for (i in seq_along(flexible_param)) {
        if (!value_set(flexible_param[[i]])) {
            pn <- names(flexible_param)[i]
            required_variables[pn] <- unit_dictionary(pn)
        }
    }
    required_variables
}

# A helping function for checking units of "optional" parameters (uses the unit
# dictionary)
require_optional_param <- function(required_variables, optional_param, data_table) {
    for (i in seq_along(optional_param)) {
        if (optional_param[[i]] %in% colnames(data_table)) {
            pn <- names(optional_param)[i]
            required_variables[pn] <- unit_dictionary(pn)
        }
    }
    required_variables
}

# A helping function that combines user-supplied and default values of `lower`,
# `upper`, and `fit_options`, as used by several fitting functions in PhotoGEA.
# This is only intended to be used internally, so checks are not provided for
# some of the inputs that cannot be affected by user actions.
assemble_luf <- function(
    param_names,
    default_lower,
    default_upper,
    default_fit_options,
    user_lower,
    user_upper,
    user_fit_options
)
{
    user_inputs <- list(user_lower, user_upper, user_fit_options)

    # Make sure the user inputs are lists
    is_list <- sapply(user_inputs, is.list)

    if (any(!is_list)) {
        stop('`lower`, `upper`, and `fit_options` must be lists')
    }

    # Make sure the lists have names
    missing_names <- sapply(user_inputs, function(x) {
        if (length(x) < 1) {
            FALSE
        } else {
            is.null(names(x))
        }
    })

    if (any(missing_names)) {
        stop('`lower`, `upper`, and `fit_options` must have named elements if they are not empty lists')
    }

    # Make sure the user list names are acceptable
    all_names <- c(
        names(user_lower),
        names(user_upper),
        names(user_fit_options)
    )

    if (any(all_names == '')) {
        stop('`lower`, `upper`, and `fit_options` must not have any empty names')
    }

    if (any(!all_names %in% param_names)) {
        stop(
            '`lower`, `upper`, and `fit_options` must only provide settings for ',
            paste0('`', param_names, '`', collapse = ', ')
        )
    }

    # Make sure the fit options are acceptable
    fit_options_okay <- sapply(user_fit_options, function(x) {
        if (value_set(x)) {
            TRUE
        } else if (is.character(x)) {
            tolower(x) %in% c('fit', 'column')
        } else {
            FALSE
        }
    })

    if (any(!fit_options_okay)) {
        stop('Each element of `fit_options` must be `fit`, `column`, or a numeric value')
    }

    # Combine the default and user lists
    default_lower[names(user_lower)] <- user_lower
    default_upper[names(user_upper)] <- user_upper
    default_fit_options[names(user_fit_options)] <- user_fit_options

    # Make sure capitalization is standardized in fit options
    default_fit_options <- lapply(default_fit_options, function(x) {
        if (value_set(x)) {
            x
        } else {
            tolower(x)
        }
    })

    # Make sure at least one parameter will be fit
    if (all(default_fit_options != 'fit')) {
        stop('No entries in `fit_options` are set to `fit`')
    }

    # Make sure elements are properly ordered
    default_lower <- default_lower[param_names]
    default_upper <- default_upper[param_names]
    default_fit_options <- default_fit_options[param_names]

    # Get a logical vector indicating which parameters should be fit
    param_to_fit <- default_fit_options == 'fit'

    # Convert the bounds to numeric vectors
    default_upper <- as.numeric(default_upper)
    default_lower <- as.numeric(default_lower)

    # Make a new version of the fit options, expressed as a numeric vector.
    # Here, any entries set to `fit` or `column` will be replaced by NA.
    fit_options_vec <- sapply(default_fit_options, function(x) {
        if (value_set(x)) {
            x
        } else {
            NA
        }
    })

    fit_options_vec <- as.numeric(fit_options_vec) # make sure names are gone

    # Return
    list(
        lower = default_lower,
        upper = default_upper,
        fit_options = default_fit_options,
        fit_options_vec = fit_options_vec,
        param_to_fit = param_to_fit
    )
}
