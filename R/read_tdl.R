read_tdl_file <- function(
    file_name,
    rows_to_skip,
    variable_name_row,
    variable_unit_row,
    data_start_row,
    timestamp_colname
)
{
    raw_tdl <- utils::read.csv(
        file_name,
        header = FALSE,
        skip = rows_to_skip
    )

    tdl_variable_names <- raw_tdl[variable_name_row - rows_to_skip,]
    tdl_variable_units <- raw_tdl[variable_unit_row - rows_to_skip,]
    tdl_data <- raw_tdl[seq(data_start_row - rows_to_skip, nrow(raw_tdl)),]

    # Convert the data to numeric values whenever possible
    tdl_data <- as.data.frame(
        lapply(tdl_data, try_as_numeric),
        stringsAsFactors = FALSE
    )

    # Set the column names
    colnames(tdl_data) <- tdl_variable_names
    colnames(tdl_variable_units) <- tdl_variable_names

    # Remove the row names that appeared after subsetting
    row.names(tdl_data) <- NULL
    row.names(tdl_variable_units) <- NULL

    # Make sure the timestamp column is properly interpreted
    tdl_data[[timestamp_colname]] <-
        as.POSIXlt(tdl_data[[timestamp_colname]], origin = "1970-01-01")

    # Make a "categories" data frame
    tdl_variable_categories <- tdl_variable_units
    tdl_variable_categories[1,] <- "Raw TDL"

    return(
        exdf(
            tdl_data,
            tdl_variable_units,
            tdl_variable_categories,
            file_name = file_name,
            rows_to_skip = rows_to_skip,
            variable_name_row = variable_name_row,
            variable_unit_row = variable_unit_row,
            data_start_row = data_start_row,
            timestamp_colname = timestamp_colname
        )
    )
}
