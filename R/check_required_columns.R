# Checks whether the required columns are in the object. This function is not
# exported to the package namespace since it is intended for internal use only.
check_required_columns <- function(x, required_columns, ...) {
    UseMethod("check_required_columns", x)
}

# For a data frame, `required_columns` should just be a character vector of
# column names that should be included in `colnames(x)`.
check_required_columns.data.frame <- function(x, required_columns, ...) {
    missing_columns <-
        required_columns[!required_columns %in% colnames(x)]

    if (length(missing_columns) > 0) {
        msg <- paste(
            "The following columns are undefined:",
            paste(missing_columns, collapse = "  ")
        )
        stop(msg)
    }
}

# For an exdf object, `required_columns` should be a list of named elements,
# where the name of each element is the name of a column that should be included
# in `x$main_data` and the value of each element is the corresponding unit for
# that column. If the required unit is NA, no check will be performed.
check_required_columns.exdf <- function(x, required_columns, ...) {
    # Use the method for data frames to check that the required columns are
    # present in `x$main_data`
    check_required_columns(x$main_data, names(required_columns), ...)

    # Now make sure the units are correct
    error_msg <- character(0)
    for (i in seq_along(required_columns)) {
        column_name <- names(required_columns)[i]
        required_units <- required_columns[[i]]
        actual_units <- x$units[[column_name]]
        if (!is.na(required_units) && actual_units != required_units) {
            error_msg <- append(
                error_msg,
                paste0(
                    "The `", column_name, '` column must have units of `',
                    required_units, "`, but its units are `", actual_units,
                    "`\n"
                )
            )
        }
    }

    if (length(error_msg) > 0) {
        stop(paste(
            "The following columns have incorrect units:\n",
            paste(error_msg, collapse = "  ")
        ))
    }
}
