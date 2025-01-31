# A helping function for `document_variables` that sets the units and category
# for one column of an exdf
set_column_info <- function(x, name, units, category) {
    if (!is.exdf(x)) {
        stop("`x` must be a exdf")
    }

    should_be_strings <- list(
        name = name,
        units = units,
        category = category
    )

    for (i in seq_along(should_be_strings)) {
        if (length(should_be_strings[[i]]) > 1) {
            stop(
                paste0(
                    "`",
                    names(should_be_strings)[[i]],
                    "` must have length 1 (input was '",
                    should_be_strings[[i]],
                    "')"
                )
            )
        }

        if (!is.character(should_be_strings[[i]])) {
            stop(
                paste0(
                    "`",
                    names(should_be_strings)[[i]],
                    "` must be a string (input was '",
                    should_be_strings[[i]],
                    "')"
                )
            )
        }
    }

    if (!name %in% colnames(x)) {
        x[['main_data']][[name]] <- NA
    }

    x[['units']][[name]] <- units
    x[['categories']][[name]] <- category

    return(x)
}

document_variables <- function(x, ...) {
    UseMethod('document_variables', x)
}

document_variables.data.frame <- function(x, ...) {
    variable_specs <- list(...)

    for (v in variable_specs) {
        if (length(v) != 3 & !is.character(v)) {
            stop("all variable specifications must be provided as vectors of three strings")
        }

        name <- v[2]

        if (!name %in% colnames(x)) {
            x[[name]] <- NA
        }
    }

    return(x)
}

document_variables.exdf <- function(x, ...)
{
    variable_specs <- list(...)

    for (v in variable_specs) {
        if (length(v) != 3 & !is.character(v)) {
            stop("all variable specifications must be provided as vectors of three strings")
        }

        x <- set_column_info(x, v[2], v[3], v[1])
    }

    return(x)
}
