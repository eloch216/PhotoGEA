# `apply_fit_function_across_reps` applies a fitting function to each replicate
# in an exdf pbject. This function is intended for internal use and is not
# exported to the PhotoGEA namespace.
#
# exdf_obj: an exdf object representing measurement sequences acquired from
#     multiple replicates
#
# f: a vector or list of vectors that can be coerced to factors, each with the
#    same number of rows as exdf_obj, that can be used to identify individual
#    replicates in the data.
#
# FUN: a function that takes a "replicate exdf object" as an input argument
#     (a subset of exdf_obj corresponding to the data from one replicate) and
#     applies a fitting procedure, returning a list with two elements:
#     `parameters`, which should be a list with named elements representing the
#     calculated fitting parameters and any identifying information that
#     describes the replicate; and `fits`, which should be a data frame
#     containing fitted values and any identifying information that describes
#     the replicate, which is usually accomplished by appending the fitted
#     values onto the original replicate data frame as a new column.
#
# ...: further arguments to FUN
#
# This function returns a list with two named elements: `parameters`, which
# includes the fitted parameters from each replicate, and `fits`, which includes
# the fits from each replicate.
apply_fit_function_across_reps <- function(exdf_obj, f, ..., FUN) {
    # Make sure exdf_obj is an exdf object
    if (!is.exdf(exdf_obj)) {
        stop("apply_fit_function_across_reps requires an exdf object")
    }

    # Split exdf_obj according to f, dropping any unused levels
    split_exdf_obj <- split(exdf_obj, f, drop = TRUE)

    # Apply FUN to each subgroup of exdf_obj
    full_results <- lapply(split_exdf_obj, function(x) {FUN(x, ...)})

    # Get all the `fits` elements, identify columns that are present in all of
    # them, and combine them
    list_of_fits <- lapply(full_results, function(x) {x$fits})

    common_columns_fits <-
        do.call(identify_common_columns, list_of_fits)

    list_of_fits <- lapply(list_of_fits, function(x) {
        x[ , common_columns_fits, return_exdf = TRUE]
    })

    fits <- do.call(rbind, list_of_fits)

    # Get all the `parameters` elements, identify columns that are present in
    # all of them, and combine them
    list_of_parameters <- lapply(full_results, function(x) {x$parameters})

    common_columns_parameters <-
        do.call(identify_common_columns, list_of_parameters)

    list_of_parameters <- lapply(list_of_parameters, function(x) {
        x[ , common_columns_parameters, return_exdf = TRUE]
    })

    parameters <- do.call(rbind, list_of_parameters)

    # Return the results
    return(list(fits = fits, parameters = parameters))
}

# `find_identifier_columns` gets the names and values of any columns in an exdf
# object that have a single unique value; these columns are taken to be
# "identifiers" that describe a replicate. This function is often used inside
# fitting functions that are passed to `apply_fit_function_across_reps` as its
# `FUN` input argument.
find_identifier_columns <- function(exdf_obj) {
    # Make sure exdf_obj is an exdf object
    if (!is.exdf(exdf_obj)) {
        stop("apply_fit_function_across_reps requires an exdf object")
    }

    # Find columns that have a single unique value
    id_column_indx <- sapply(
        colnames(exdf_obj),
        function(x) {length(unique(exdf_obj[ , x])) == 1}
    )

    # Make an exdf object that includes the unique value of each such column
    id_columns <- exdf_obj[1, id_column_indx, return_exdf = TRUE]
}
