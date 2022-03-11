identify_common_licor_columns <- function(exdf_objs, verbose) {
    if (length(exdf_objs) < 1) {
        stop("There must be at least one Licor file")
    }

    # Get a set of all the columns from all the exdf objects, formatted to
    # include the name, category, and units (in other words, the "fancy column
    # name")
    all_column_list <- lapply(
        exdf_objs,
        function(x) {
            data.frame(
                column = fancy_column_names(x),
                filename = x[['file_name']]
            )
        }
    )

    all_column_df <- do.call(rbind, all_column_list)

    # For each unique fancy column name, get a list of all the Licor files where
    # it can be found
    column_file_info <- by(
        all_column_df,
        all_column_df$column,
        function(x) {
            x$filename
        }
    )

    # Get a vector of all the fancy column names that are found in all of the
    # files
    common_column_indx <- sapply(column_file_info, length) == length(exdf_objs)
    common_column_names <- names(common_column_indx)[common_column_indx]

    # If required, print some additional information to the user
    if (verbose) {
        # Get a vector of the fancy column names that are not found in all of the
        # files
        other_column_indx <- sapply(column_file_info, length) < length(exdf_objs)
        other_column_names <- names(other_column_indx)[other_column_indx]
        other_column_info <- column_file_info[other_column_names]

        # Create and send a message about the columns were chosen or rejected
        msg <- "The following columns were found in all of the Licor files:\n "
        msg <- append(msg, paste(common_column_names, collapse = "\n  "))


        msg <- append(msg, "\n\nThe following columns were only found in some of the Licor files:\n ")

        for (i in seq_along(other_column_info)) {
            msg <- append(msg, paste(names(other_column_info)[i], "\n   "))
            msg <- append(msg, paste(column_file_info[[i]], collapse = "\n    "))
            msg <- append(msg, "\n")
        }

        msg <- append(msg, "\n")
        cat(msg)
    }

    # Put the common names back to the same order they have in the Licor files
    columns_from_file <- fancy_column_names(exdf_objs[[1]])
    ordered_common_column_names <-
        columns_from_file[columns_from_file %in% common_column_names]

    # Truncate the names (so they are no longer fancy)
    common_columns <- gsub(" .+$", "", ordered_common_column_names)

    return(common_columns)
}
