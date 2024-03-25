conversions_from_text <- function(txt) {
    utils::read.table(
        textConnection(txt),
        header = TRUE,
        na.strings = '',
        stringsAsFactors = FALSE
    )
}

gasex_column_conversions <- conversions_from_text('
    instrument_type original_name original_units original_category new_name new_units       new_category multiplier
    "Licor LI-6800" "PhiPS2"      "NA"           "FLR"             "PhiPS2" "dimensionless" "FLR"        1
    "Licor LI-6800" "PhiPs2"      "NA"           "FLR"             "PhiPS2" "dimensionless" "FLR"        1
')

standardize_gasex_columns <- function(instrument_type, exdf_obj) {
    # Get conversion definitions for this instrument type
    conversion_definitions <-
        gasex_column_conversions[gasex_column_conversions$instrument_type == instrument_type, ]

    # Apply conversions
    for (cn in colnames(exdf_obj)) {
        if (cn %in% conversion_definitions$original_name) {
            # Get the conversion definitions for a column with this name
            cd <- conversion_definitions[conversion_definitions$original_name == cn, ]

            # Check to see if the column's units and categories match the
            # conversion definition
            units_match <- exdf_obj$units[[cn]] == cd$original_units
            categories_match <- exdf_obj$categories[[cn]] == cd$original_category

            if (units_match && categories_match) {
                # Apply the conversion, and change the units and category for
                # this column
                exdf_obj[, cn]            <- cd$multiplier * exdf_obj[, cn]
                exdf_obj$units[[cn]]      <- cd$new_units
                exdf_obj$categories[[cn]] <- cd$new_category

                # Get the new column names
                cnames <- colnames(exdf_obj$main_data)
                cnames[cnames == cn] <- cd$new_name

                # Apply new column names
                colnames(exdf_obj$main_data)  <- cnames
                colnames(exdf_obj$units)      <- cnames
                colnames(exdf_obj$categories) <- cnames
            }
        }
    }

    # Return
    exdf_obj
}
