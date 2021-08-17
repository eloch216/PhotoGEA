# Constructor for the exdf (extended data frame) class, which is similar to a
# data.frame but is extended by specifying units and a category for each column
# in addition to names
exdf <- function(main_data, units, categories, ...) {
    # Get ready to store messages about any problems with the inputs
    errors <- character()

    # Check to make sure the main_data, units, and categories inputs are data
    # frames with the same column names
    should_be_data_frames <- list(
        main_data = main_data,
        units = units,
        categories = categories
    )

    main_data_colnames <- colnames(main_data)

    for (i in seq_along(should_be_data_frames)) {
        if (!is.data.frame(should_be_data_frames[[i]])) {
            errors <- append(
                errors,
                paste0(
                    "'",
                    names(should_be_data_frames)[i],
                    "' must be a data.frame"
                )
            )
        }
        if (!identical(colnames(should_be_data_frames[[i]]), main_data_colnames)) {
            errors <- append(
                errors,
                paste0(
                    "'",
                    names(should_be_data_frames)[i],
                    "' must have the same column names as main_data"
                )
            )
        }
    }

    # Check to make sure the units and categories each have one row
    should_have_one_row <- list(
        units = units,
        categories = categories
    )

    for (i in seq_along(should_have_one_row)) {
        if (nrow(should_have_one_row[[i]]) != 1) {
            errors <- append(
                errors,
                paste0(
                    "'",
                    names(should_have_one_row)[i],
                    "' must have exactly one row"
                )
            )
        }
    }

    # Make the exdf object
    new_exdf <- list(
        main_data = main_data,
        units = units,
        categories = categories
    )

    # Add any other properties specified by the input arguments
    extra_properties <- list(...)
    for (i in seq_along(extra_properties)) {
        new_exdf[[names(extra_properties)[i]]] <- extra_properties[[i]]
    }

    # Specify the class and return the exdf
    class(new_exdf) <- c("exdf", class(new_exdf))
    return(new_exdf)
}

# Check if an object is an exdf
is.exdf <- function(x) {
    inherits(x, "exdf")
}

# Convert an exdf to a data.frame (used by `View`, `write.csv`, and others)
as.data.frame.exdf <- function(x, ...) {
    # Make sure time columns are properly formatted for displaying
    main_data <- x[['main_data']]
    for (col in colnames(main_data)) {
        if ("POSIXlt" %in% class(main_data[[col]])) {
            main_data[[col]] <- format(main_data[[col]])
        }
    }

    # Store the categories and units as the first rows of the data.frame
    return(
        rbind(
            x[['categories']],
            x[['units']],
            main_data
        )
    )
}

# Define a helper function for making nice column names (used to improve `print`
# and `str`)
fancy_column_names <- function(x) {
    paste0(
        colnames(x[['units']]),
        " [", x[['categories']][1,], "]",
        " (", x[['units']][1,], ")"
    )
}

# Print an exdf
print.exdf <- function(x) {
    res <- x[['main_data']]
    colnames(res) <- fancy_column_names(x)
    print(res)
}

# Display the structure of an exdf
str.exdf <- function(x) {
    res <- x[['main_data']]
    colnames(res) <- fancy_column_names(x)
    str(res)
}

# Get the length of an exdf
length.exdf <- function(x) {
    length(x[['main_data']])
}

# Get the dimensions of an exdf
dim.exdf <- function(x) {
    dim(x[['main_data']])
}

# Get the dimension names of an exdf
dimnames.exdf <- function(x) {
    dimnames(x[['main_data']])
}

# Access elements of an exdf's main_data
`[.exdf` <- function(x, i, j) {
    return(x[['main_data']][i,j])
}

# Modify elements of an exdf's main_data
`[<-.exdf` <- function(x, i, j, value) {
    if (!is.character(j)) {
        stop(
            paste(
                "when modifying a exdf using [i,j], the column `j` must be",
                "specified as a name rather than an index"
            )
        )
    }

    if (!j %in% colnames(x)) {
        x[['units']][1,j] <- "NA"
        x[['categories']][1,j] <- "NA"
    }

    x[['main_data']][i,j] <- value
    return(x)
}

# Set the units and category for a column of an exdf
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

# Combine exdf objects by the rows of their main_data
rbind.exdf <- function(
    ...,
    deparse.level = 1,
    make.row.names = TRUE,
    stringsAsFactors = default.stringsAsFactors(),
    factor.exclude = TRUE
)
{
    exdf_list <- list(...)

    if (length(exdf_list) < 2) {
        stop("rbind requires two or more exdf objects")
    }

    type_check <- rep(FALSE, times = length(exdf_list))
    for (i in seq_along(exdf_list)) {
        if (is.exdf(exdf_list[[i]])) {
            type_check[i] <- TRUE
        }
    }

    if (!all(type_check)) {
        stop(
            paste(
                "exdf objects can only be combined with other exdf",
                "objects when using rbind"
            )
        )
    }

    main_data_list <- list()
    for (i in seq_along(exdf_list)) {
        main_data_list[[i]] <- exdf_list[[i]][['main_data']]
    }

    new_main_data <- do.call("rbind", main_data_list)

    return(
        exdf(
            new_main_data,
            exdf_list[[1]][['units']],
            exdf_list[[1]][['categories']]
        )
    )
}
