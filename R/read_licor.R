read_licor_file <- function(
    file_name,
    preamble_data_rows,
    variable_category_row,
    variable_name_row,
    variable_unit_row,
    data_start_row,
    timestamp_colname
)
{
    # Define a helping function for reading one row (with column names) from an
    # Excel file, where the data row is specified and the column names are
    # assumed to be in the preceding row
    read_row_wc <- function(row_num)
    {
        row_data <- openxlsx::readWorkbook(
            file_name,
            colNames = TRUE,
            skipEmptyCols = FALSE,
            rows = c(row_num - 1, row_num)
        )
    }

    # Get the preamble data rows
    licor_preamble <- lapply(preamble_data_rows, read_row_wc)

    # Define a helping function for reading one row (without column names) from
    # an Excel file
    read_row_nc <- function(row_num)
    {
        row_data <- openxlsx::readWorkbook(
            file_name,
            colNames = FALSE,
            skipEmptyCols = FALSE,
            rows = row_num
        )
    }

    # Define a helping function for converting a row to a data frame
    row_to_data_frame <- function(rowdata)
    {
        data_frame <- data.frame(
            matrix(
                ncol = length(rowdata),
                nrow = 1
            ),
            stringsAsFactors = FALSE
        )
        data_frame[1,] <- rowdata
        return(data_frame)
    }

    # Define a helping function for reading a row, replacing any Unicode
    # characters, and converting the result to a data frame
    get_processed_row_data_frame <- function(row_num)
    {
        row_data <- read_row_nc(row_num)
        row_data <- replace_unicode(row_data)
        row_data <- row_to_data_frame(row_data)
        return(row_data)
    }

    # Get the column types, names and units, replacing any problematic Unicode
    # characters and converting the results to data frames when required
    licor_variable_categories <- get_processed_row_data_frame(variable_category_row)
    licor_variable_names <- get_processed_row_data_frame(variable_name_row)
    licor_variable_units <- get_processed_row_data_frame(variable_unit_row)

    # Get the main data
    licor_data <- openxlsx::readWorkbook(
        file_name,
        startRow = data_start_row,
        colNames = FALSE,
        skipEmptyCols = FALSE
    )

    # Convert the data to numeric values whenever possible
    licor_data <- as.data.frame(
        lapply(licor_data, try_as_numeric),
        stringsAsFactors = FALSE
    )

    # Apply column names
    colnames(licor_data) <- licor_variable_names[1,]
    colnames(licor_variable_categories) <- licor_variable_names[1,]
    colnames(licor_variable_units) <- licor_variable_names[1,]

    # Make sure the timestamp column is properly interpreted
    licor_data[[timestamp_colname]] <-
        as.POSIXlt(licor_data[[timestamp_colname]], origin = "1970-01-01")

    return(
        exdf(
            licor_data,
            licor_variable_units,
            licor_variable_categories,
            file_name = file_name,
            preamble = licor_preamble,
            preamble_data_rows = preamble_data_rows,
            variable_category_row = variable_category_row,
            variable_name_row = variable_name_row,
            variable_unit_row = variable_unit_row,
            data_start_row = data_start_row
        )
    )
}
