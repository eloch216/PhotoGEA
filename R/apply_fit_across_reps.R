# `apply_fit_function_across_reps` applies a fitting function to each replicate
# in a data frame.
#
# main_data_frame: a data frame of measurement sequences acquired from multiple
#     replicates
#
# replicate_column_name: a string specifying the replicate column in the data
#     frame
#
# FUN: a function that takes a "replicate data frame" as an input argument
#     (a subset of the "main data frame" corresponding to the data from one
#     replicate) and applies a fitting procedure, returning a list with two
#     elements: `parameters`, which should be a list with named elements
#     representing the calculated fitting parameters and any identifying
#     information that describes the replicate; and `fits`, which should
#     be a data frame containing fitted values and any identifying information
#     that describes the replicate, which is usually accomplished by appending
#     the fitted values onto the original replicate data frame as a new column.
#
# ...: further arguments to FUN
#
# This function returns a list with two named elements: `parameters`, which
# includes the fitted parameters from each replicate, and `fits`, which includes
# the fits from each replicate.
#
apply_fit_function_across_reps <- function(
    main_data_frame,
    replicate_column_name,
    ...,
    FUN
)
{
    # Use `by` to apply the function to subsets of the main data frame, broken
    # up according to the values in the replicate column. The result should be a
    # list where each element is itself a list with two elements corresponding
    # to the output of FUN.
    full_results <- by(
        main_data_frame,
        main_data_frame[[replicate_column_name]],
        ...,
        FUN = FUN
    )

    # Use `lapply` to make a list of just the `fits` element from each entry in
    # the `full_results` list, and then use `do.call(rbind, args)` to combine
    # them into one data frame.
    fits <- do.call(
        rbind,
        lapply(
            full_results,
            function(x) {x[['fits']]}
        )
    )

    # Apply the same procedure to the `parameters` element.
    parameters <- do.call(
        rbind,
        lapply(
            full_results,
            function(x) {as.data.frame(x[['parameters']])}
        )
    )

    # Remove any row names that appeared while processing
    rownames(fits) <- NULL
    rownames(parameters) <- NULL

    return(list(fits = fits, parameters = parameters))
}

# `find_identifier_columns` gets the names and values of any columns in the data
# frame that have a single unique value; these columns are taken to be
# "identifiers" that describe a replicate. This function is often used inside
# fitting functions that are passed to `apply_fit_function_across_reps` as its
# `FUN` input argument.
find_identifier_columns <- function(data_frame) {
    id_column_indx <- sapply(
        colnames(data_frame),
        function(x) {length(unique(data_frame[[x]]))==1}
    )
    id_columns <- unique(data_frame[id_column_indx])
}
