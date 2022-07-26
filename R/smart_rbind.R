# Given a list of data frames, `smart_rbind` determines which columns are common
# to all the data frames, restricts each data frame to the common columns, and
# binds all of them together using `rbind`.
smart_rbind <- function(list_of_data_frames) {
    # Remove any NULL entries or those with zero rows
    list_of_data_frames <-
        list_of_data_frames[!sapply(list_of_data_frames, function(x) {
            is.null(x) || nrow(x) < 1
        })]

    # Find the column names that are common to all the elements of
    # `list_of_data_frames`
    common_column_names <- do.call(identify_common_columns, list_of_data_frames)

    # Restrict the data frames in the list to the common columns
    list_of_data_frame_subsets <-
        lapply(list_of_data_frames, function(x) {x[common_column_names]})

    # Use rbind to combine the elements of `list_of_data_frames`
    result <- do.call(rbind, list_of_data_frame_subsets)

    # Remove any row names that pop up and return the result
    row.names(result) <- NULL
    return(result)
}
