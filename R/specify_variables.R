# A helping function for `specify_variables` that sets the units and category
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

specify_variables <- function(exdf_obj, ...)
{
    if (!is.exdf(exdf_obj)) {
        stop("exdf_obj must be an exdf object")
    }

    variable_specs <- list(...)

    for (v in variable_specs) {
        if (length(v) != 3 & !is.character(v)) {
            stop("all variable specifications must be provided as vectors of three strings")
        }

        exdf_obj <- set_column_info(
            exdf_obj,
            v[2],
            v[3],
            v[1]
        )
    }

    return(exdf_obj)
}
