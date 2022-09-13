organize_response_curve_data <- function(
    licor_exdf,
    identifier_columns,
    measurement_numbers_to_remove,
    column_for_ordering,
    ordering_column_tolerance = Inf
)
{
    if (!is.exdf(licor_exdf)) {
        stop('organize_response_curve_data requires an exdf object')
    }

    # Make sure the column for ordering is defined
    required_variables <- list()
    required_variables[[column_for_ordering]] <- NA

    check_required_variables(licor_exdf, required_variables)

    # Create factors from the identifier columns
    f <- lapply(identifier_columns, function(x) {licor_exdf[ , x]})

    # Add `seq_num` values for each curve
    licor_exdf <- do.call(rbind, by(licor_exdf, f, function(x) {
        x[, 'seq_num'] <- seq_len(nrow(x))
        return(x)
    }))

    # Make a subset of the full result that only includes the desired
    # measurement points
    licor_exdf <-
        licor_exdf[!licor_exdf[, 'seq_num'] %in% measurement_numbers_to_remove, , TRUE]

    # Check to make sure the identifier columns specify curves with the same
    # number of points. Don't check for any infinite values since we don't know
    # which columns should be ignored for this check.
    check_licor_data(
      licor_exdf,
      identifier_columns,
      driving_column = column_for_ordering,
      driving_column_tolerance = ordering_column_tolerance,
      col_to_ignore_for_inf = NULL
    )

    # Make sure the data is ordered properly for plotting
    for (cn in c(column_for_ordering, identifier_columns)) {
        licor_exdf <- licor_exdf[order(licor_exdf[, cn]), , TRUE]
    }

    # Remove any row names and return the exdf
    row.names(licor_exdf$main_data) <- NULL
    return(licor_exdf)
}
