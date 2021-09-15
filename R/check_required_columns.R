# Checks whether the required columns are in the exdf object (should also work
# for data frames)
check_required_columns <- function(exdf_obj, required_columns) {
    missing_columns <-
        required_columns[!required_columns %in% colnames(exdf_obj)]

    if (length(missing_columns) > 0) {
        msg <- paste(
            "The following columns are undefined:",
            paste(missing_columns, collapse = " ")
        )
        stop(msg)
    }
}
