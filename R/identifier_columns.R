identifier_columns <- function(x) {
    UseMethod("identifier_columns", x)
}

identifier_columns.data.frame <- function(x) {
    # Find columns that have a single unique value
    id_column_indx <-
        sapply(colnames(x), function(cn) {length(unique(x[ , cn])) == 1})

    # Restrict the data frame to the first row of just these columns
    id_columns <- x[1, id_column_indx]

    # Remove row names and return the identifiers
    rownames(id_columns) <- NULL
    return(id_columns)
}

identifier_columns.exdf <- function(x) {
    # Find identifier columns using the data frame method
    id_column_df <- identifier_columns(x$main_data)

    # Restrict the exdf object to the first row of just those columns
    id_columns <- x[1, colnames(id_column_df), return_exdf = TRUE]

    # Remove row names and return the identifiers
    rownames(id_columns$main_data) <- NULL
    return(id_columns)
}
