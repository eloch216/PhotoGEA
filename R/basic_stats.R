basic_stats_chunk <- function(exdf_chunk) {
    # Find the identifier columns and their names
    id_columns <- identifier_columns(exdf_chunk)
    id_column_names <- colnames(id_columns)

    # Find the names of the non-id columns
    all_column_names <- colnames(exdf_chunk)
    other_column_names <-
        all_column_names[!all_column_names %in% id_column_names]

    # Get a subset of the chunk limited to the non-id columns so we can get
    # their units, and categories
    other_columns <- exdf_chunk[ , other_column_names, TRUE]

    # Get the mean of all the non-id columns
    mean_list <- lapply(other_column_names, function(x) {
        column <- exdf_chunk[ , x]
        if (is.numeric(column)) {
            mean(exdf_chunk[ , x])
        } else {
            NA
        }
    })
    names(mean_list) <- paste0(other_column_names, '_avg')

    # Get the standard deviation of all the non-id columns
    sd_list <- lapply(other_column_names, function(x) {
        column <- exdf_chunk[ , x]
        if (is.numeric(column)) {
            stats::sd(exdf_chunk[ , x])
        } else {
            NA
        }
    })

    # Get the standard error of the mean for all the non-id columns
    stderr_list <- lapply(sd_list, function(x) {x / sqrt(nrow(exdf_chunk))})
    names(stderr_list) <- paste0(other_column_names, '_stderr')

    # Create new exdf objects for the means and standard errors
    mean_exdf <- exdf(
        as.data.frame(mean_list),
        stats::setNames(other_columns$units, names(mean_list)),
        stats::setNames(other_columns$categories, names(mean_list))
    )

    stderr_exdf <- exdf(
        as.data.frame(stderr_list),
        stats::setNames(other_columns$units, names(stderr_list)),
        stats::setNames(other_columns$categories, names(stderr_list))
    )

    stats_exdf <- cbind(mean_exdf, stderr_exdf)
    stats_exdf <- stats_exdf[ , order(colnames(stats_exdf)), TRUE]

    return(cbind(id_columns, stats_exdf))
}

basic_stats <- function(
    exdf_obj,
    identifier_columns
)
{
    if (!is.exdf(exdf_obj)) {
        stop("basic_stats requires an exdf object")
    }

    # Make sure the identifier columns are defined
    required_variables <- list()
    for (cn in identifier_columns) {
        required_variables[[cn]] <- NA
    }

    check_required_variables(exdf_obj, required_variables)

    # Split the exdf object by the identifiers
    f <- lapply(identifier_columns, function(x) {exdf_obj[ , x]})

    split_exdf <- split(exdf_obj, f, drop = TRUE)

    # Calculate the basic stats
    stats_list <- lapply(split_exdf, basic_stats_chunk)

    # Restrict to common column names
    common_columns <- do.call(identify_common_columns, stats_list)

    stats_list <- lapply(stats_list, function(x) {x[ , common_columns, TRUE]})

    # Combine all the individual exdf objects and return them
    new_exdf_obj <- do.call(rbind, stats_list)
    row.names(new_exdf_obj$main_data) <- NULL
    return(new_exdf_obj)
}
