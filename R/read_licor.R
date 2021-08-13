# This script includes functions for reading Licor data files. Since many of the
# inputs to the `read_licor_file` function below are unlikely to change, default
# values are also provided here.
#
# ------------------------------------------------------------------------------
#
# These functions require the `openxlsx` library, which can be installed using
# the following command if it is not already installed:
#
# install.packages('openxlsx')
#
# ------------------------------------------------------------------------------
#
# IMPORTANT NOTE ABOUT LICOR EXCEL FILES: by default, Licor Excel files do not
# `calculate` formula values. This causes a problem when reading them in R,
# since any data entry determined from a formula will be read as 0. To fix this
# issue for a Licor Excel file, open it in in Excel, go to the `Formulas` menu,
# and choose `Calculate Now`. (Alternatively, press F9.) Then save the file and
# close it. See https://github.com/tidyverse/readxl/issues/495 for more details.

library(openxlsx)

# Specify information about how the data in the Licor spreadsheet is arranged
PREAMBLE_DATA_ROWS <- c(3, 5, 7, 9, 11, 13)  # The rows in a Licor Excel file that contain the preamble data
VARIABLE_TYPE_ROW <- 14                      # The row in a Licor Excel file that contains the variable types
VARIABLE_NAME_ROW <- 15                      # The row in a Licor Excel file that contains the variable names
VARIABLE_UNIT_ROW <- 16                      # The row in a Licor Excel file that contains the variable units
DATA_START_ROW <- 17                         # The first main data row in a Licor Excel file

# Specify ASCII replacements for specific Unicode characters or sequences of
# Unicode characters. This list is not intended to be exhaustive, but does
# include all the Unicode characters present in a few Licor data files on
# 2021-04-26.
UNICODE_REPLACEMENTS <- data.frame(
    rbind(
        c("<ce><94>",               "Delta"),
        c("<e2><81><bb><c2><b2>",   "^(-2)"),
        c("<e2><81><bb><c2><b9>",   "^(-1)"),
        c("<c2><b5>",               "micro"),
        c("<c2><b0>",               "degrees "),
        c("<c2><b2>",               "^2")
    ),
    stringsAsFactors = FALSE
)
colnames(UNICODE_REPLACEMENTS) <- c("pattern", "replacement")

# read_licor_file: a function for reading the data from a Licor Excel file into
# an R list.
#
# In R, some problems occur for column names or units that include Unicode
# characters such as Greek letters. We use a combination of the `iconv` and
# `gsub` functions to replace them by equivalent phrases, following a strategy
# discussed on
# https://stackoverflow.com/questions/36108790/trouble-with-strings-with-u0092-unicode-characters
#
# ------------------------------------------------------------------------------
#
# INPUTS:
#
# - file_name: a relative or absolute path to an Excel file containing Licor
#              data
#
# - preamble_data_rows: a numeric vector whose entries indicate the rows in the
#                       Licor excel file that contain the preamble information
#
# - variable_type_row: the row number in the Licor Excel file containing the
#                      variable type information, e.g. "GasEx", "FLR", etc
#
# - variable_name_row: the row number in the Licor Excel file containing the
#                      variable names, e.g. "A", "Ci", etc
#
# - variable_unit_row: the row number in the Licor Excel file containing the
#                      variable units, e.g. "Pa", "s", etc
#
# - data_start_row: the first row number of the table containing the measured
#                   data in the Licor Excel file
#
# - timestamp_colname: the name of the column that contains the timestamp of
#                      each measurement
#
# ------------------------------------------------------------------------------
#
# OUTPUT:
#
# a list with the following named elements that together fully include all the
# data from the Licor Excel file:
#
# - file_name: a copy of the input argument with the same name
#
# - preamble: a list of data frames representing each line in the preamble,
#             where each data frame has one row and named columns. No Unicode
#             replacement is performed for the preamble data.
#
# - types: a data frame with one row and named columns, where the column names
#          are the variable names and the row describes the type of each
#          variable
#
# - units: a data frame with one row and named columns, where the column names
#          are the variable names and the row describeds the units of each
#          variable
#
# - main_data: a data frame with named columns, where the column names are the
#              variable names and the rows represent the measurements that were
#              stored in the original Licor Excel file
#
# - preamble_data_rows: a copy of the input argument with the same name
#
# - variable_type_row: a copy of the input argument with the same name
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
    variable_type_row,
    variable_name_row,
    variable_unit_row,
    data_start_row,
    timestamp_colname = 'time'
)
{
    # Define a helping function for reading one row (with column names) from an
    # Excel, where the data row is specified and the column names are assumed to
    # be in the preceding row
    read_row_wc <- function(row_num)
    {
        row_data <- readWorkbook(
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
        row_data <- readWorkbook(
            file_name,
            colNames = FALSE,
            skipEmptyCols = FALSE,
            rows = row_num
        )
    }

    # Define a helping function for replacing Unicode characters
    replace_unicode <- function(strings)
    {
        strings <- iconv(
            strings,
            "",
            "ASCII",
            "byte"
        )

        for (i in seq_along(UNICODE_REPLACEMENTS[['pattern']])) {
            strings <- gsub(
                x = strings,
                pattern = UNICODE_REPLACEMENTS[['pattern']][i],
                replacement = UNICODE_REPLACEMENTS[['replacement']][i]
            )
        }

        return(strings)
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
    licor_variable_types <- get_processed_row_data_frame(variable_type_row)
    licor_variable_names <- get_processed_row_data_frame(variable_name_row)
    licor_variable_units <- get_processed_row_data_frame(variable_unit_row)

    # Get the main data
    licor_data <- readWorkbook(
        file_name,
        startRow = data_start_row,
        colNames = FALSE,
        skipEmptyCols = FALSE
    )

    # Apply column names
    colnames(licor_data) <- licor_variable_names[1,]
    colnames(licor_variable_types) <- licor_variable_names[1,]
    colnames(licor_variable_units) <- licor_variable_names[1,]

    # Make sure the timestamp column is properly interpreted
    licor_data[[timestamp_colname]] <-
        as.POSIXlt(licor_data[[timestamp_colname]], origin = "1970-01-01")

    return(
        list(
            file_name = file_name,
            preamble = licor_preamble,
            types = licor_variable_types,
            units = licor_variable_units,
            main_data = licor_data,
            preamble_data_rows = preamble_data_rows,
            variable_type_row = variable_type_row,
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
#               Licor data
#
# All other inputs are identical to those in the `read_licor_file` function.
#
# ------------------------------------------------------------------------------
#
# OUTPUT:
#
# a list with unnamed elements, each of which is a list describing the contents
# of a Licor file that was created using the the `read_licor_file` function.
#
batch_read_licor_file <- function(
    file_names,
    preamble_data_rows,
    variable_type_row,
    variable_name_row,
    variable_unit_row,
    data_start_row
)
{
    lapply(
        file_names,
        function(filename) {
            read_licor_file(
                    filename,
                    preamble_data_rows,
                    variable_type_row,
                    variable_name_row,
                    variable_unit_row,
                    data_start_row
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

    choose.files(
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
