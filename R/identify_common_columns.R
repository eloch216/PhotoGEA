identify_common_columns <- function(...) {
    if (length(list(...)) < 1) {
        stop("Identify_common_columns requires at least one input object")
    }

    UseMethod("identify_common_columns", list(...)[[1]])
}

identify_common_columns.data.frame <- function(...) {
    df_list <- list(...)

    # Make sure all the inputs are data frames
    type_check <- sapply(df_list, is.data.frame)
    if (!all(type_check)) {
        stop("All input objects must be data frames")
    }

    # Create a list where each entry is a data frame of all the column names
    # from one element of `df_list`
    all_column_list <-
        lapply(df_list, function(x) {data.frame(n = colnames(x))})

    # Combine the column names into a data frame with one column (`n`)
    all_column_df <- do.call(rbind, all_column_list)

    # Find the number of times each column name is repeated
    all_column_lengths <- by(all_column_df, all_column_df$n, nrow)

    # Find the column names that are common to all the elements of `df_list`
    common_column_names <-
        names(all_column_lengths)[all_column_lengths == length(df_list)]

    # Put the common names back to the same order they have in the first data
    # frame
    first_columns <- colnames(df_list[[1]])

    ordered_common_column_names <-
        first_columns[first_columns %in% common_column_names]

    return(ordered_common_column_names)
}

identify_common_columns.exdf  <- function(...) {
    exdf_list <- list(...)

    # Make sure all the inputs are exdf objects
    type_check <- sapply(exdf_list, is.exdf)
    if (!all(type_check)) {
        stop("All input objects must be exdf objects")
    }

    # Create a list where each entry is a data frame of all the "fancy" column
    # names from one element of `exdf_list`
    all_column_list <-
        lapply(exdf_list, function(x) {data.frame(n = fancy_column_names(x))})

    # Combine the column names into a data frame with one column (`n`)
    all_column_df <- do.call(rbind, all_column_list)

    # Find the number of times each column name is repeated
    all_column_lengths <- by(all_column_df, all_column_df$n, nrow)

    # Find the column names that are common to all the elements of `exdf_list`
    common_column_names <-
        names(all_column_lengths)[all_column_lengths == length(exdf_list)]

    # Truncate the names (so they are no longer fancy)
    common_column_names <- gsub(" [[].+$", "", common_column_names)

    # Put the common names back to the same order they have in the first data
    # frame
    first_columns <- colnames(exdf_list[[1]])

    ordered_common_column_names <-
        first_columns[first_columns %in% common_column_names]

    return(ordered_common_column_names)
}

simple_exdf_1 <- exdf(data.frame(A = 1), data.frame(A = 'u'), data.frame(A = 'c'))
simple_exdf_2 <- exdf(data.frame(A = 2), data.frame(A = 'u'), data.frame(A = 'c'))

identify_common_columns.exdf(simple_exdf_1, simple_exdf_2)
