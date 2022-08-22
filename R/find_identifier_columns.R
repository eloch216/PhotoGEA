# `find_identifier_columns` gets the names and values of any columns in an exdf
# object that have a single unique value; these columns are taken to be
# "identifiers" that describe a replicate. This function is often used inside
# fitting functions that are passed to `by.exdf` as its `FUN` input argument.
find_identifier_columns <- function(exdf_obj) {
    # Make sure exdf_obj is an exdf object
    if (!is.exdf(exdf_obj)) {
        stop("find_identifier_columns requires an exdf object")
    }

    # Find columns that have a single unique value
    id_column_indx <- sapply(
        colnames(exdf_obj),
        function(x) {length(unique(exdf_obj[ , x])) == 1}
    )

    # Make an exdf object that includes the unique value of each such column
    id_columns <- exdf_obj[1, id_column_indx, return_exdf = TRUE]
}
