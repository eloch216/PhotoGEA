factorize_id_column <- function(x, ...) {
    UseMethod("factorize_id_column", x)
}

factorize_id_column.character <- function(x, ...) {
    # Trim any leading or trailing whitespace from the id values
    x <- trimws(x)

    # Standardize capitalization of any leading WT values
    x <- gsub('^[Ww][Tt]', 'WT', x)

    # Get all the unique ids
    identifiers <- unique(x)

    # Split each identifier into three units:
    # - An initial WT (with any capitalization), if present
    # - A trailing numeric value, if present
    # - Any other remaining content
    id_wt_part <- sapply(identifiers, function(id) {
        if (grepl('^WT', id)) {
            'WT'
        } else {
            NA
        }
    })

    id_num_part <- as.numeric(sapply(identifiers, function(id) {
        fmatch <- regexpr('[0-9.]+$', id)

        if (fmatch >= 0) {
            substring(id, fmatch)
        } else {
            NA
        }
    }))

    id_mid_part <- sapply(identifiers, function(id) {
        id <- gsub('^[Ww][Tt]', '', id)
        gsub('[0-9.]+$', '', id)
    })

    # Now get an order for the ids, determined by the three separate units
    sorted_identifiers <-
        identifiers[order(id_wt_part, id_mid_part, id_num_part)]

    # Convert the ids to a factor with the correct order and return
    factor(x, levels = sorted_identifiers)
}

factorize_id_column.data.frame <- function(x, id_column_name, ...) {
    # Make sure the required variables are defined
    check_required_variables(x, id_column_name)

    # If the column is already a factor, turn it back into a character vector
    if (is.factor(x[[id_column_name]])) {
        x[[id_column_name]] <- as.character(x[[id_column_name]])
    }

    # Convert the data frame column to a factor with properly sorted names
    x[[id_column_name]] <- factorize_id_column.character(x[[id_column_name]])

    return(x)
}

factorize_id_column.exdf <- function(x, id_column_name, ...) {
    # Convert the id column in the main_data data frame
    x[['main_data']] <-
        factorize_id_column.data.frame(x[['main_data']], id_column_name)

    return(x)
}
