# This script includes functions for reading Licor data files.
#
# ------------------------------------------------------------------------------
#
# IMPORTANT NOTE ABOUT LICOR EXCEL FILES: by default, Licor Excel files do not
# `calculate` formula values. This causes a problem when reading them in R,
# since any data entry determined from a formula will be read as 0. To fix this
# issue for a Licor Excel file, open it in in Excel, go to the `Formulas` menu,
# and choose `Calculate Now`. (Alternatively, press F9.) Then save the file and
# close it. See https://github.com/tidyverse/readxl/issues/495 for more details.

# read_licor_file: a function for reading the data from a Licor Excel file into
# an exdf object.
#
# ------------------------------------------------------------------------------
#
# INPUTS:
#
# - file_name: a relative or absolute path to an Excel file containing Licor
#       data
#
# - preamble_data_rows: a numeric vector whose entries indicate the rows in the
#       Licor excel file that contain the preamble information
#
# - variable_category_row: the row number in the Licor Excel file containing the
#       variable category information, e.g. "GasEx", "FLR", etc
#
# - variable_name_row: the row number in the Licor Excel file containing the
#       variable names, e.g. "A", "Ci", etc
#
# - variable_unit_row: the row number in the Licor Excel file containing the
#       variable units, e.g. "Pa", "s", etc
#
# - data_start_row: the first row number of the table containing the measured
#       data
#
# - timestamp_colname: the name of the column that contains the timestamp of
#       each measurement
#
# ------------------------------------------------------------------------------
#
# OUTPUT:
#
# an exdf object with the following "extra" elements that fully includes all the
# data from the Licor Excel file:
#
# - file_name: a copy of the input argument with the same name
#
# - preamble: a list of data frames representing each line in the preamble,
#             where each data frame has one row and named columns. No Unicode
#             replacement is performed for the preamble data.
#
# - preamble_data_rows: a copy of the input argument with the same name
#
# - variable_category_row: a copy of the input argument with the same name
#
# - variable_name_row: a copy of the input argument with the same name
#
# - variable_unit_row: a copy of the input argument with the same name
#
# - data_start_row: a copy of the input argument with the same name
#
# - timestamp_colname: a copy of the input argument with the same name
#
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

# batch_read_licor_file: a function for reading the data from multiple Licor
# Excel files into an R list.
#
# ------------------------------------------------------------------------------
#
# INPUTS:
#
# - file_names: a vector of relative or absolute paths to Excel files containing
#       Licor data
#
# All other inputs are identical to those in the `read_licor_file` function.
#
# ------------------------------------------------------------------------------
#
# OUTPUT:
#
# a list with unnamed elements, each of which is an exdf object representing the
# contents of a Licor file
#
batch_read_licor_file <- function(
    file_names,
    preamble_data_rows,
    variable_category_row,
    variable_name_row,
    variable_unit_row,
    data_start_row,
    timestamp_colname
)
{
    lapply(
        file_names,
        function(filename) {
            read_licor_file(
                filename,
                preamble_data_rows,
                variable_category_row,
                variable_name_row,
                variable_unit_row,
                data_start_row,
                timestamp_colname
            )
        }
    )
}

# choose_input_licor_files: a function for interactively selecting multiple
# Licor Excel files.
#
# Important note: this function is only available on MS Windows.
#
# ------------------------------------------------------------------------------
#
# INPUTS:
#
# None
#
# ------------------------------------------------------------------------------
#
# OUTPUT:
#
# a vector of file name strings representing absolute paths to Licor Excel files
#
choose_input_licor_files <- function()
{
    if (!interactive() | .Platform$OS.type != "windows") {
        stop(
            paste(
                "The `choose_input_licor_files` function is only available in",
                "interactive R sessions running in MS Windows"
            )
        )
    }

    utils::choose.files(
        default = "",
        caption = "Select Licor Excel input files",
        multi = TRUE,
        filters = matrix(
            c(
                "Excel files (*.xlsx)", "All files (*.*)",
                "*.xlsx", "*.*"
            ),
            ncol = 2
        ),
        index = 1
    )
}
