control_regex_string <- function(control_group_name) {
    character_list <- sapply(seq_len(nchar(control_group_name)), function(i) {
        thischar <- substr(control_group_name, i, i)
        paste0('[', toupper(thischar), tolower(thischar), ']')
    })

    paste(c('^', character_list), collapse = '')
}

factorize_id_column <- function(x, ...) {
    UseMethod("factorize_id_column", x)
}

factorize_id_column.character <- function(
    x,
    control_group_name = 'WT',
    ...
)
{
    # Trim any leading or trailing whitespace from the id values
    x <- trimws(x)

    # Get a regex for finding it at the start of a string
    control_group_regex <- control_regex_string(control_group_name)

    # Standardize capitalization of any leading WT values
    x <- gsub(control_group_regex, control_group_name, x)

    # Get all the unique ids
    identifiers <- unique(x)

    # Split each identifier into three units:
    # - An initial WT (with any capitalization), if present
    # - A trailing numeric value, if present
    # - Any other remaining content
    id_wt_part <- sapply(identifiers, function(id) {
        if (grepl(control_group_regex, id)) {
            control_group_name
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
        id <- gsub(control_group_regex, '', id)
        gsub('[0-9.]+$', '', id)
    })

    # Now get an order for the ids, determined by the three separate units
    sorted_identifiers <-
        identifiers[order(id_wt_part, id_mid_part, id_num_part)]

    # Convert the ids to a factor with the correct order and return
    factor(x, levels = sorted_identifiers)
}

factorize_id_column.data.frame <- function(
    x,
    id_column_name,
    control_group_name = 'WT',
    ...
)
{
    # Make sure the required variables are defined
    check_required_variables(x, id_column_name)

    # If the column is already a factor, turn it back into a character vector
    if (is.factor(x[[id_column_name]])) {
        x[[id_column_name]] <- as.character(x[[id_column_name]])
    }

    # Convert the data frame column to a factor with properly sorted names
    x[[id_column_name]] <- factorize_id_column.character(
        x[[id_column_name]],
        control_group_name
    )

    return(x)
}

factorize_id_column.exdf <- function(
    x,
    id_column_name,
    control_group_name = 'WT',
    ...
)
{
    # Convert the id column in the main_data data frame
    x[['main_data']] <- factorize_id_column.data.frame(
        x[['main_data']],
        id_column_name,
        control_group_name
    )

    return(x)
}
