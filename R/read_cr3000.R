read_cr3000 <- function(
    file_name,
    rows_to_skip = 1,
    variable_name_row = 2,
    variable_unit_row = 3,
    data_start_row = 5,
    ...
)
{
    raw_cr <- utils::read.csv(
        file_name,
        header = FALSE,
        skip = rows_to_skip,
        stringsAsFactors = FALSE,
        ...
    )

    cr_variable_names <-
        make.unique(as.character(raw_cr[variable_name_row - rows_to_skip,]))

    cr_variable_units <- raw_cr[variable_unit_row - rows_to_skip,]
    cr_data <- raw_cr[seq(data_start_row - rows_to_skip, nrow(raw_cr)),]

    # Convert the data to numeric values whenever possible
    cr_data <- as.data.frame(
        lapply(cr_data, try_as_numeric),
        stringsAsFactors = FALSE
    )

    # Set the column names
    colnames(cr_data) <- cr_variable_names
    colnames(cr_variable_units) <- cr_variable_names

    # Remove the row names that appeared after subsetting
    row.names(cr_data) <- NULL
    row.names(cr_variable_units) <- NULL

    # Make a "categories" data frame
    cr_variable_categories <- data.frame(
        matrix(nrow = 1, ncol = ncol(cr_variable_units), data = 'read_cr3000'),
        stringsAsFactors = FALSE
    )

    colnames(cr_variable_categories) <- cr_variable_names

    return(
        exdf(
            cr_data,
            cr_variable_units,
            cr_variable_categories,
            rows_to_skip = rows_to_skip,
            variable_name_row = variable_name_row,
            variable_unit_row = variable_unit_row,
            data_start_row = data_start_row
        )
    )
}
