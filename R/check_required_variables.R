# Checks whether the required variables are in the object.
check_required_variables <- function(x, required_variables, check_NA = TRUE) {
    UseMethod('check_required_variables', x)
}

# For a data frame, `required_variables` can just be a character vector of
# column names that should be included in `colnames(x)`.
check_required_variables.data.frame <- function(x, required_variables, check_NA = TRUE) {
    # If required_variables is a list, just use the names of its elements
    if (is.list(required_variables)) {
        required_variables <- names(required_variables)
    }

    missing_columns <-
        required_variables[!required_variables %in% colnames(x)]

    if (length(missing_columns) > 0) {
        msg <- paste(
            'The following required columns are not present:',
            paste(missing_columns, collapse = ', ')
        )
        stop(msg)
    }

    na_columns <- sapply(required_variables, function(rv) {
        all(is.na(x[[rv]]))
    })

    if (check_NA && any(na_columns)) {
        msg <- paste(
            'The following required columns are all NA:',
            paste(required_variables[na_columns], collapse = ', ')
        )
        stop(msg)
    }
}

# For an exdf object, `required_variables` should be a list of named elements,
# where the name of each element is the name of a column that should be included
# in `x$main_data` and the value of each element is the corresponding unit for
# that column. If the required unit is NA, no check will be performed.
check_required_variables.exdf <- function(x, required_variables, check_NA = TRUE) {
    # Use the method for data frames to check that the required variables are
    # present in `x$main_data`
    check_required_variables(x$main_data, names(required_variables), check_NA)

    # Now make sure the units are correct
    error_msg <- character(0)
    for (i in seq_along(required_variables)) {
        column_name <- names(required_variables)[i]
        required_units <- required_variables[[i]]
        actual_units <- x$units[[column_name]]
        if (!is.na(required_units) && (is.na(actual_units) || actual_units != required_units)) {
            error_msg <- append(
                error_msg,
                paste0(
                    'The `', column_name, '` column must have units of `',
                    required_units, '`, but its units are `', actual_units,
                    '`\n'
                )
            )
        }
    }

    if (length(error_msg) > 0) {
        stop(paste(
            'The following columns have incorrect units:\n',
            paste(error_msg, collapse = ', ')
        ))
    }
}
