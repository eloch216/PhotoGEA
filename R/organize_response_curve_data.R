organize_response_curve_data <- function(
    exdf_obj,
    measurement_number_name,
    num_obs_in_seq,
    measurement_numbers_to_keep,
    independent_column_name,
    rep_column_name,
    event_column_name
)
{
    if (!is.exdf(exdf_obj)) {
        stop("organize_response_curve_data requires an exdf object")
    }

    # Add a `seq_num` column, where a value of `i` means that this row is the
    # `ith` point along a response curve
    exdf_obj[, 'seq_num'] <-
        ((exdf_obj[, measurement_number_name] - 1) %% num_obs_in_seq) + 1

    # Make a subset of the full result that only includes the desired
    # measurement points
    data_subset <-
        exdf_obj[exdf_obj[, 'seq_num'] %in% measurement_numbers_to_keep, , return_exdf = TRUE]

    # Make sure the data is ordered properly for plotting
    data_subset <- data_subset[order(data_subset[, independent_column_name]), return_exdf = TRUE]
    data_subset <- data_subset[order(data_subset[, rep_column_name]), return_exdf = TRUE]
    data_subset <- data_subset[order(data_subset[, event_column_name]), return_exdf = TRUE]

    row.names(data_subset$main_data) <- NULL

    return(data_subset)
}
