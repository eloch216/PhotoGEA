# This script includes functions for writing Licor data to Excel files, as well
# as a few default settings.
#
# ------------------------------------------------------------------------------
#
# These functions require the `openxlsx` library, which can be installed using
# the following command if it is not already installed:
#
# install.packages('openxlsx')

# add_licor_sheet: a function that adds one sheet (representing one Licor data
# file) to an Excel workbook object
#
# ------------------------------------------------------------------------------
#
# INPUTS:
#
# - wb: a workbook object
#
# - sheetnum: an integer specifying the sheet to create
#
# - sheetname: a string representing the name of the sheet. Note: the name will
#              be truncated to 31 characters, since that is the limit for the
#              length of sheet names in Excel files
#
# - licor_file: a list representing the Licor information to write, typically
#               produced using one of the functions in `read_licor.R` and
#               possibly modified by one or more of the functions in
#               `licor_data_operations.R`
#
# ------------------------------------------------------------------------------
#
# OUTPUT:
#
# a list with the following named elements that together fully include all the
# data from the Licor Excel file:
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
add_licor_sheet <- function(
    wb,
    sheetnum,
    sheetname,
    licor_file,
    variable_info_offset,
    main_data_offset
)
{
    # Make sure the sheet name complies with the Excel limitation of 31
    # characters
    truncated_sheetname <- substr(sheetname, 1, 31)

    # Add a new sheet to the workbook
    addWorksheet(
        wb = wb,
        sheetName = truncated_sheetname
    )

    # Make a helping function for writing one row to the sheet
    write_row <- function(data_to_write, start_row)
    {
        writeData(
            wb,
            sheet = sheetnum,
            startRow = start_row,
            data_to_write,
            colNames = FALSE
        )
    }

    # Write the preamble to the sheet
    preamble <<- licor_file[['preamble']]
    for (i in seq_along(licor_file[['preamble_data_rows']])) {
        write_row(
            as.list(colnames(data.frame(
                preamble[i],
                stringsAsFactors = FALSE
            ))),
            licor_file[['preamble_data_rows']][i] - 1
        )
        write_row(
            preamble[i],
            licor_file[['preamble_data_rows']][i]
        )
    }

    # Write the variable info to the sheet
    write_row(
        licor_file[['types']],
        licor_file[['variable_type_row']] + variable_info_offset
    )

    write_row(
        data.frame(
            as.list(colnames(licor_file[['types']])),
            stringsAsFactors = FALSE
        ),
        licor_file[['variable_name_row']] + variable_info_offset
    )

    write_row(
        licor_file[['units']],
        licor_file[['variable_unit_row']] + variable_info_offset
    )

    # Write the main data to the sheet
    extra_offset <- 1
    writeData(
        wb,
        sheet = sheetnum,
        startRow = licor_file[['data_start_row']] +
            variable_info_offset + main_data_offset,
        licor_file[['main_data']],
        colNames = FALSE
    )
}

# write_excel_formula: a function for writing formula values to a column in a
# Licor excel file.
#
# ------------------------------------------------------------------------------
#
# INPUTS:
#
# - licor_file: a list representing the data from a Licor file (typically
#               produced by a a call to the `read_licor_file` defined in
#               `read_licor.R`)
#
# - variable: the name of a variable in the Licor file
#
# - formula_values: a vector of formula values to include in an Excel file. The
#                   values should not begin with the "=" typically required for
#                   Excel formulas
#
# ------------------------------------------------------------------------------
#
# OUTPUT:
#
# a list describing Licor data having the same structure as the `licor_file`
# input
#
write_excel_formula <- function(
    licor_file,
    variable,
    formula_values
)
{
    licor_file[['main_data']][[variable]] <- formula_values

    class(licor_file[['main_data']][[variable]]) <- c(
        class(licor_file[['main_data']][[variable]]),
        "formula"
    )

    return(licor_file)
}

batch_write_excel_formula <- function(
    licor_files,
    variable,
    formula_values
)
{
    lapply(
        licor_files,
        function(licor_file) {
            write_excel_formula(
                    licor_file,
                    variable,
                    formula_values
            )
        }
    )
}

# sheet_name_from_file_name: a function for determining an acceptable sheet name
# from the name of a Licor Excel file.
#
# This function may be helpful when writing sheets based on individual files.
#
# ------------------------------------------------------------------------------
#
# INPUTS:
#
# - filename: the name of a Licor Excel file
#
# ------------------------------------------------------------------------------
#
# OUTPUT:
#
# a string that could be used as a sheet name in an Excel document, determined
# by removing path information and the .xlsx extension from the input filename
# and taking only the final 31 characters of the name.
#
sheet_name_from_file_name <- function(filename)
{
    # Remove any path information from the filename
    filename <- basename(filename)

    # Remove the extension
    filename <- tools::file_path_sans_ext(filename)

    # Keep only the last 31 characters
    return(substr(filename, max(nchar(filename)-30,1), nchar(filename)))
}

# choose_output_licor_file: a function for interactively selecting a single
# Licor Excel file.
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
# a string representing an absolute path to a Licor Excel file
#
choose_output_licor_file <- function()
{
    if (!interactive() | .Platform$OS.type != "windows") {
        stop(
            paste(
                "The `choose_output_licor_file` function is only available in",
                "interactive R sessions running in MS Windows"
            )
        )
    }

    choose.files(
        default = "",
        caption = "Create a Licor Excel file to store the output",
        multi = FALSE,
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
