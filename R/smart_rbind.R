# Given a list of data frames, `smart_rbind` determines which columns are common
# to all the data frames, restricts each data frame to the common columns, and
# binds all of them together using `rbind`.
smart_rbind <- function(list_of_data_frames) {
    # Remove any NULL entries or those with zero rows
    list_of_data_frames <-
        list_of_data_frames[!sapply(list_of_data_frames, function(x) {
            is.null(x) || nrow(x) < 1
        })]

    # Create a list where each entry is a data frame of all the column names
    # from one element of `list_of_data_frames`
    all_column_list <-
        lapply(list_of_data_frames, function(x) {data.frame(n = names(x))})

    # Combine the column names into a data frame with one column (`n`)
    all_column_df <- do.call(rbind, all_column_list)

    # Find the number of times each column name is repeated
    all_column_lengths <- by(all_column_df, all_column_df$n, nrow)

    # Find the column names that are common to all the elements of
    # `list_of_data_frames`
    common_column_names <-
        names(all_column_lengths)[all_column_lengths == length(list_of_data_frames)]

    # Restrict the data frames in the list to the common columns
    list_of_data_frame_subsets <-
        lapply(list_of_data_frames, function(x) {x[common_column_names]})

    # Use rbind to combine the elements of `list_of_data_frames`
    result <- do.call(rbind, list_of_data_frame_subsets)

    # Remove any row names that pop up and return the result
    row.names(result) <- NULL
    return(result)
}
