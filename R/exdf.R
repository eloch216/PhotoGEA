# This file defines a constructor for the exdf class and some essential S3
# methods. An exdf object is similar to a data frame but is extended by
# specifying units and a category for each column in addition to names.

# Constructor
exdf <- function(main_data, units = NULL, categories = NULL, ...) {
    # Make sure `main_data` is a data frame
    if (!is.data.frame(main_data)) {
        stop("`main_data` must be a data frame")
    }

    # If `units` or `categories` is NULL, replace it with a default data frame
    # that has the same names as `main_data` and one row where all entries are
    # NA.
    if (is.null(units)) {
        units <- main_data[1, ]
        units[1, ] <- NA
        row.names(units) <- NULL
    }

    if (is.null(categories)) {
        categories <- main_data[1, ]
        categories[1, ] <- NA
        row.names(categories) <- NULL
    }

    # Get ready to store messages about any problems with the inputs
    errors <- character()

    # Check to make sure the main_data, units, and categories inputs are data
    # frames with the same column names, unless they are NULL
    should_be_data_frames <- list(
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
is.exdf <- function(x, consistency_check = FALSE) {
    # Make sure the class of x includes `exdf`
    class_check <- inherits(x, "exdf")

    if (!class_check) {
        return(FALSE)
    }

    if (consistency_check) {
        # Make sure x has elements with the required names
        required_elements <- c('main_data', 'units', 'categories')
        element_check <- all(sapply(
            required_elements,
            function(ele) {ele %in% names(x)}
        ))

        if (!element_check) {
            warning('x must have elements named `main_data`, `units`, and `categories`')
            return(FALSE)
        }

        # Make sure the required elements are all data frames
        type_check <- all(sapply(
            required_elements,
            function(ele) {is.data.frame(x[[ele]])}
        ))

        if (!element_check) {
            warning('`x$main_data`, `x$units`, and `x$categories` must be data frames')
            return(FALSE)
        }

        # Make sure the required elements have the same column names
        columns <- colnames(x$main_data)
        column_check <- all(sapply(
            c('units', 'categories'),
            function(ele) {identical(colnames(x[[ele]]), columns)}
        ))

        if (!column_check) {
            warning('`x$main_data`, `x$units`, and `x$categories` must have the same column names')
            return(FALSE)
        }

        # Make sure the units and categories have just one row
        row_check <- all(sapply(
            c('units', 'categories'),
            function(ele) {nrow(x[[ele]] == 1)}
        ))

        if (!row_check) {
            warning('`x$units` and `x$categories` must each have exactly one row')
            return(FALSE)
        }
    }

    return(TRUE)
}

# Convert an exdf to a data.frame (used by `View`, `write.csv`, and others)
as.data.frame.exdf <- function(x, ...) {
    # Make sure time columns are properly formatted for displaying
    main_data <- x$main_data
    for (col in colnames(main_data)) {
        if ("POSIXlt" %in% class(main_data[[col]])) {
            main_data[[col]] <- format(main_data[[col]])
        }
    }

    # Store the categories and units as the first rows of the data.frame
    return(rbind(x$categories, x$units, main_data))
}

# Define a helper function for making nice column names (used to improve `print`
# and `utils::str`)
fancy_column_names <- function(x) {
    paste0(
        colnames(x$units),
        " [", x$categories[1,], "]",
        " (", x$units[1,], ")"
    )
}

# Print an exdf
print.exdf <- function(x, ...) {
    res <- x$main_data
    colnames(res) <- fancy_column_names(x)
    print(res, ...)
}

# Display the structure of an exdf
str.exdf <- function(object, ...) {
    res <- object$main_data
    colnames(res) <- fancy_column_names(object)
    utils::str(res, ...)
}

# Get the length of an exdf
length.exdf <- function(x) {
    length(x$main_data)
}

# Get the dimensions of an exdf
dim.exdf <- function(x) {
    dim(x$main_data)
}

# Get the dimension names of an exdf
dimnames.exdf <- function(x) {
    dimnames(x$main_data)
}

# Set the dimension names of an exdf; don't attempt to change the row names for
# the units and categories.
`dimnames<-.exdf` <- function(x, value) {
    dimnames(x$main_data) <- value

    unit_dimnames <- dimnames(x$units)
    unit_dimnames[[2]] <- value[[2]]
    dimnames(x$units) <- unit_dimnames

    category_dimnames <- dimnames(x$categories)
    category_dimnames[[2]] <- value[[2]]
    dimnames(x$categories) <- category_dimnames

    return(x)
}

# Access elements of an exdf
`[.exdf` <- function(x, i, j, return_exdf = FALSE) {
    if (return_exdf) {
        exdf_elements <- names(x)

        essential_elements <- c("main_data", "units", "categories")

        extra_elements <-
            exdf_elements[!exdf_elements %in% essential_elements]

        extra_element_list <- stats::setNames(
            lapply(extra_elements, function(ele) {return(x[[ele]])}),
            extra_elements
        )

        essential_element_list <- list(
            x$main_data[i, j],
            x$units[1, j],
            x$categories[1, j]
        )

        return(do.call(exdf, c(essential_element_list, extra_element_list)))
    }
    else {
        return(x$main_data[i, j])
    }
}

# Modify elements of an exdf's main_data
`[<-.exdf` <- function(x, i, j, value) {
    if (!is.character(j)) {
        stop("when modifying a exdf using [i,j], the column `j` must be specified as a name rather than an index")
    }

    if (!j %in% colnames(x)) {
        x$units[ ,j] <- "NA"
        x$categories[ ,j] <- "NA"
    }

    if (is.null(value)) {
        x$units[ , j] <- NULL
        x$categories[ , j] <- NULL
    }

    x$main_data[i, j] <- value
    return(x)
}

# Combine exdf objects by the columns of their main_data, units, and categories
cbind.exdf <- function(..., deparse.level = 1) {
    exdf_list <- list(...)

    # Make sure there is one or more exdf object
    if (length(exdf_list) < 1) {
        stop("rbind.exdf requires one or more exdf objects")
    }

    # Make sure all the objects are indeed exdf objects
    type_check <- lapply(exdf_list, function(x) {is.exdf(x)})

    if (!all(as.logical(type_check))) {
        stop("exdf objects can only be combined with other exdf objects when using cbind")
    }

    # Get the main_data, units, and categories from each exdf object
    main_data_list <- lapply(exdf_list, function(x) {x$main_data})
    units_list <- lapply(exdf_list, function(x) {x$units})
    categories_list <- lapply(exdf_list, function(x) {x$categories})

    # Make a new exdf object by combining them all with cbind
    return(exdf(
        do.call(cbind, c(main_data_list, list(deparse.level = deparse.level))),
        do.call(cbind, c(units_list, list(deparse.level = deparse.level))),
        do.call(cbind, c(categories_list, list(deparse.level = deparse.level)))
    ))
}

# Combine exdf objects by the rows of their main_data
rbind.exdf <- function(
    ...,
    deparse.level = 1,
    make.row.names = TRUE,
    stringsAsFactors = FALSE,
    factor.exclude = TRUE
)
{
    exdf_list <- list(...)

    # Make sure there is one or more exdf object
    if (length(exdf_list) < 1) {
        stop("rbind.exdf requires one or more exdf objects")
    }

    # Make sure all the objects are indeed exdf objects
    type_check <- lapply(exdf_list, function(x) {is.exdf(x)})

    if (!all(as.logical(type_check))) {
        stop("exdf objects can only be combined with other exdf objects when using rbind")
    }

    # Get the info from the first file to use as a reference
    first_exdf <- exdf_list[[1]]

    # Check to make sure all the exdf objects have the same variables, units,
    # and categories
    for (i in seq_along(exdf_list)) {
        current_exdf <- exdf_list[[i]]

        if (!identical(colnames(first_exdf$main_data), colnames(current_exdf$main_data))) {
            cat('\ncolnames from first exdf object:\n')
            print(colnames(first_exdf$main_data))
            cat('\ncolnames from current exdf object:\n')
            colnames(current_exdf$main_data)
            cat('\n')
            stop("exdf objects must all have the same column names when using rbind")
        }

        if (!identical(first_exdf$categories, current_exdf$categories)) {
            cat('\ncategories from first exdf object:\n')
            print(colnames(first_exdf$categories))
            cat('\ncategories from current exdf object:\n')
            colnames(current_exdf$categories)
            cat('\n')
            stop("exdf objects must all have the same categories when using rbind")
        }

        if (!identical(first_exdf$units, current_exdf$units)) {
            cat('\nunits from first exdf object:\n')
            print(colnames(first_exdf$units))
            cat('\nunits from current exdf object:\n')
            colnames(current_exdf$units)
            cat('\n')
            stop("exdf objects must all have the same units when using rbind")
        }
    }

    # Use rbind for data frames to combine all the main_data elements
    main_data_list <- list()
    for (i in seq_along(exdf_list)) {
        main_data_list[[i]] <- exdf_list[[i]]$main_data
    }

    new_main_data <- do.call(
        rbind,
        c(
            main_data_list,
            list(
                deparse.level = deparse.level,
                make.row.names = make.row.names,
                stringsAsFactors = stringsAsFactors,
                factor.exclude = factor.exclude
            )
        )
    )

    return(exdf(new_main_data, first_exdf$units, first_exdf$categories))
}

# Split an exdf object into a list of smaller exdf objects by the value of one
# or more factors
split.exdf <- function(x, f, drop = FALSE, lex.order = FALSE, ...)
{
    f <- interaction(f, drop = drop, lex.order = lex.order)

    lapply(
        split(x = seq_len(nrow(x)), f = f, drop = drop, ...),
        function(ind) x[ind, , return_exdf = TRUE]
    )
}

# Split an exdf object into chunks by the value of one or more factors, apply
# FUN to each of the chunks, and return the output of each call to FUN as one
# element of a list
by.exdf <- function(data, INDICES, FUN, ...)
{
    split_exdf_obj <- split(data, INDICES, drop = TRUE)
    lapply(split_exdf_obj, function(x) {FUN(x, ...)})
}
