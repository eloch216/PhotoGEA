organize_response_curve_data <- function(
    licor_exdf,
    identifier_columns,
    measurement_numbers_to_remove,
    column_for_ordering,
    ordering_column_tolerance = Inf,
    columns_to_average = c(),
    print_information = TRUE
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

    if (length(measurement_numbers_to_remove) > 0) {
        # Make a subset of the full result that only includes the desired
        # measurement points
        licor_exdf <-
            licor_exdf[!licor_exdf[, 'seq_num'] %in% measurement_numbers_to_remove, , TRUE]

        # Check to make sure the identifier columns specify curves with the same
        # number of points. Don't check for any infinite values since we don't
        # know which columns should be ignored for this check.
        check_response_curve_data(
          licor_exdf,
          identifier_columns,
          driving_column = column_for_ordering,
          driving_column_tolerance = ordering_column_tolerance,
          col_to_ignore_for_inf = NULL,
          print_information = print_information
        )
    }

    # Make sure the data is ordered properly for plotting
    for (cn in c(column_for_ordering, identifier_columns)) {
        licor_exdf <- licor_exdf[order(licor_exdf[, cn]), , TRUE]
    }

    # Re-create factors from the identifier columns
    f_new <- lapply(identifier_columns, function(x) {licor_exdf[ , x]})

    # Calulate average values of certain columns
    columns_to_average <-
        columns_to_average[columns_to_average %in% colnames(licor_exdf)]

    licor_exdf <- do.call(rbind, by(licor_exdf, f_new, function(x) {
        for (cn in columns_to_average) {
            avg_cn <- paste0(cn, '_avg')
            x <- set_variable(
                x,
                avg_cn,
                x[['units']][[cn]],
                'organize_response_curve_data',
                mean(x[, cn], na.rm = TRUE)
            )
        }
        x
    }))

    # Remove any row names and return the exdf
    row.names(licor_exdf$main_data) <- NULL
    return(licor_exdf)
}
