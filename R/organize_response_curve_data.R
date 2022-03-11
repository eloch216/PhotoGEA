organize_response_curve_data <- function(
    full_data_set,
    measurement_number_name,
    num_obs_in_seq,
    measurement_numbers_to_keep,
    independent_column_name,
    rep_column_name,
    event_column_name
)
{
    # Add a `seq_num` column, where a value of `i` means that this row is the
    # `ith` point along a response curve
    full_data_set$seq_num <-
        ((full_data_set[[measurement_number_name]] - 1) %% num_obs_in_seq) + 1

    # Make a subset of the full result that only includes the desired
    # measurement points
    data_subset <-
        full_data_set[full_data_set$seq_num %in% measurement_numbers_to_keep,]

    # Make sure the data is ordered properly for plotting
    data_subset <- data_subset[order(data_subset[[independent_column_name]]),]
    data_subset <- data_subset[order(data_subset[[rep_column_name]]),]
    data_subset <- data_subset[order(data_subset[[event_column_name]]),]

    row.names(data_subset) <- NULL

    return(data_subset)
}
