factorize_id_column <- function(full_data_set, id_column_name) {
    # Make sure the required columns are defined
    check_required_columns(full_data_set, id_column_name)

    if (is.factor(full_data_set[[id_column_name]])) {
        # This column is already a factor, so turn it back into a character
        # vector
        full_data_set[[id_column_name]] <-
            as.character(full_data_set[[id_column_name]])
    }

    # Get all the identifiers
    identifiers <- unique(full_data_set[[id_column_name]])

    # Split the identifiers into those that begin with "WT" and those that don't
    wt_identifiers <- grep("^WT", identifiers, value = TRUE)
    non_wt_identifiers <- identifiers[!identifiers %in% wt_identifiers]

    # Order the WT identifiers according to the value with the initial "WT"
    # removed, interpreted as a numeric value when possible
    wt_sorted <-
        wt_identifiers[order(try_as_numeric(gsub("^WT", "", wt_identifiers)))]

    # Order the non-WT identifiers according to the word before the first space,
    # interpreted as a numeric value when possible
    non_wt_sorted <-
        non_wt_identifiers[order(try_as_numeric(gsub(" .+$", "", non_wt_identifiers)))]

    # Get the full set of sorted identifiers
    sorted_identifiers <- c(wt_sorted, non_wt_sorted)

    # Convert the data frame column to a factor using the sorted names
    full_data_set[[id_column_name]] <- factor(
        full_data_set[[id_column_name]],
        levels = sorted_identifiers
    )

    return(full_data_set)
}
