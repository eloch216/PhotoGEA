set_variable <- function(
    exdf_obj,
    name,
    units = NULL,
    category = NULL,
    value = NA,
    id_column = NULL,
    value_table = list()
)
{
    if (!is.exdf(exdf_obj)) {
        stop("set_variable requires an exdf object")
    }

    # Set the value of the column
    exdf_obj[, name] <- value

    # If units or a category were provided, document the column info
    if (!is.null(units) || !is.null(category)) {
        if (is.null(units)) {
            units <- 'NA'
        }

        if (is.null(category)) {
            category <- 'NA'
        }

        exdf_obj <- document_variables(exdf_obj, c(category, name, units))
    }

    # If an id_column was provided, set column values based on the table
    if (!is.null(id_column)) {
        required_variables <- list()
        required_variables[[id_column]] <- NA
        check_required_variables(exdf_obj, required_variables)

        for (i in seq_along(value_table)) {
            exdf_obj[as.character(exdf_obj[, id_column]) == names(value_table)[i], name] <- value_table[[i]]
        }
    }

    return(exdf_obj)
}
