remove_points <- function(exdf_obj, ..., method = 'remove') {
    if (!is.exdf(exdf_obj)) {
        stop('remove_points requires an exdf object')
    }

    if (!tolower(method) %in% c('remove', 'exclude')) {
        stop('`method` must be either `remove` or `exclude`')
    }

    # Create a list of the optional input arguments
    arg_list <- list(...)

    # Make sure the optional input arguments are all lists
    list_check <- sapply(arg_list, function(x) {!is.list(x)})
    if (any(list_check)) {
        stop('optional arguments to remove_points must be lists')
    }

    # Make sure the optional input arguments have names
    name_check <- sapply(arg_list, function(x) {is.null(names(x))})
    if (any(name_check)) {
        stop('optional arguments to remove_points must have names')
    }

    # Go through each set of conditions to remove the desired points
    for (point_description in arg_list) {
        # Make sure the exdf object contains the specified columns
        required_variables <-
            lapply(point_description, function(x) {return(NA)})

        check_required_variables(exdf_obj, required_variables)

        # If we are using the "exclude" mode, initialize the
        # `include_when_fitting` column if the exdf object does not already have
        # one
        if (method == 'exclude' && !'include_when_fitting' %in% colnames(exdf_obj)) {
            exdf_obj <- set_variable(
                exdf_obj,
                'include_when_fitting',
                category = 'remove_points',
                value = TRUE
            )
        }

        # Initialize the logical vector of points to keep
        points_to_keep <- rep.int(FALSE, nrow(exdf_obj))

        # Apply all the conditions, adding back any points that don't meet the
        # criteria
        for (i in seq_along(point_description)) {
            name <- names(point_description)[i]
            condition <- point_description[[i]]
            points_to_keep <- points_to_keep | !(exdf_obj[ , name] %in% condition)
        }

        # Truncate the exdf to just the points that don't meet the condition, or
        # exclude those points from any subsequent fits
        exdf_obj <- if (tolower(method) == 'remove') {
            exdf_obj[points_to_keep, , TRUE]
        } else {
            exdf_obj[!points_to_keep, 'include_when_fitting'] <- FALSE
            exdf_obj
        }
    }

    return(exdf_obj)
}
