set_variable <- function(
    data_table,
    name,
    units = NULL,
    category = NULL,
    value = NA,
    id_column = NULL,
    value_table = NULL
)
{
    # Check for some potential issues with the inputs
    if (!is.null(value_table)) {
        if (!is.list(value_table)) {
            stop('When a `value_table` is supplied, it must be a list')
        }

        tnames <- names(value_table)

        if (is.null(tnames) || any(tnames == '')) {
            stop('When a `value_table` is supplied, all of its elements must have names')
        }
    }

    if (is.null(id_column) && !is.null(value_table)) {
        stop('When a `value_table` is supplied, an `id_column` must also be supplied')
    }

    # Set the value of the column
    data_table[, name] <- value

    # If units or a category were provided, document the column info
    if (!is.null(units) || !is.null(category)) {
        if (is.null(units)) {
            units <- 'NA'
        }

        if (is.null(category)) {
            category <- 'NA'
        }

        data_table <- document_variables(data_table, c(category, name, units))
    }

    # If an id_column was provided, set column values based on the table
    if (!is.null(id_column)) {
        required_variables <- list()
        required_variables[[id_column]] <- NA
        check_required_variables(data_table, required_variables, check_NA = FALSE)

        for (i in seq_along(value_table)) {
            data_table[as.character(data_table[, id_column]) == names(value_table)[i], name] <- value_table[[i]]
        }
    }

    return(data_table)
}
