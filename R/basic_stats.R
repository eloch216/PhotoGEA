basic_stats_column <- function(column, name, npts) {
    if(!is.numeric(column)) {
        # This is not a numeric column, so just return an empty list
        return(list())
    } else {
        # Calculate basic stats for the column
        column_avg <- mean(column)
        column_sd <- sd(column)
        column_stderr <- column_sd / sqrt(npts)

        # Return the results as a list
        return(
            setNames(
                list(column_avg, column_sd, column_stderr),
                paste0(name, c("_avg", "_sd", "_stderr"))
            )
        )
    }
}

basic_stats_chunk <- function(chunk) {
    # Find the identifier columns and their names
    id_columns <- find_identifier_columns(chunk)
    id_column_names <- names(id_columns)

    # Find the names of the non-id columns
    all_column_names <- colnames(chunk)
    other_column_names <-
        all_column_names[!all_column_names %in% id_column_names]

    # Get the number of rows in the chunk
    npts <- nrow(chunk)

    # Calculate basic stats for each non-id column
    stats <- lapply(
        other_column_names,
        function(x) {basic_stats_column(chunk[[x]], x, npts)}
    )

    # Return everything as a data frame
    return(as.data.frame(do.call(c, c(id_columns, stats, list(npts = npts)))))
}

basic_stats <- function(
    dataframe,
    factor_column_names,
    columns_to_exclude = c('obs', 'elapsed')
)
{
    factor_columns <- lapply(factor_column_names, function(x) {dataframe[[x]]})

    subdataframe <- dataframe[,!(names(dataframe) %in% columns_to_exclude)]

    stats <- smart_rbind(by(subdataframe, factor_columns, basic_stats_chunk))

    for (i in seq_along(factor_column_names)) {
        stats <- stats[order(stats[[factor_column_names[i]]]),]
    }

    row.names(stats) <- NULL

    return(stats)
}
