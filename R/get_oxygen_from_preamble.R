get_oxygen_from_preamble <- function(licor_exdf) {
    # Check the input
    if (!is.exdf(licor_exdf)) {
        stop('licor_exdf must be an exdf object')
    }

    # Try to get the oxygen values
    oxygen_column_name <- 'Oxygen'

    oxygen_values <- if (oxygen_column_name %in% colnames(licor_exdf)) {
        licor_exdf[, oxygen_column_name]
    } else {
        NA
    }

    # Try to get the oxygen category
    oxygen_category <- if (oxygen_column_name %in% colnames(licor_exdf)) {
        licor_exdf[['categories']][1, oxygen_column_name]
    } else {
        NULL
    }

    # Set the units, retaining any previously defined values and category
    set_variable(
        licor_exdf,
        oxygen_column_name,
        'percent',
        oxygen_category,
        oxygen_values
    )
}
