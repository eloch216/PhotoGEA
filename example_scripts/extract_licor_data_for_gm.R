# IMPORTANT NOTE: this script has been superseded by `gm_from_tdl.R` and should
# not be used, since it continues to use Excel files for making calculations.

# This script uses the PhotoGEA and openxlsx libraries to read Licor data from
# multiple files, extract certain columns, add some new columns, fill the new
# columns with Excel formulas, and finally write the information from each file
# as a separate tab within an output Excel spreadsheet. This script has been
# tested on Windows and MacOS.
#
# ------------------------------------------------------------------------------
#
# IMPORTANT NOTE ABOUT LICOR EXCEL FILES: by default, Licor Excel files do not
# `calculate` formula values. This causes a problem when reading them in R,
# since any data entry determined from a formula will be read as 0. To fix this
# issue for a Licor Excel file, open it in in Excel, go to the `Formulas` menu,
# and choose `Calculate Now`. (Alternatively, press F9.) Then save the file and
# close it. See https://github.com/tidyverse/readxl/issues/495 for more details.
#
# ------------------------------------------------------------------------------
#
# This script is broken up into several sections to make it easier to use:
# - Components that might need to change each time this script is run
# - Components that are less likely to change each time this script is run
# - Functions used to load and process the data (shouldn't need to change)
# - The command that actually calls the functions (shouldn't need to change)
#
# Typically, it should only be necessary to specify the names of input files
# and the output file. This information is specified in the FILES_TO_PROCESS
# vector and the OUTPUT_FILENAME string. If CHOOSE_FILES_INTERACTIVELY is set to
# true, these file names can be chosen interactively via a dialog box (only
# available on MS Windows).
#
# The filenames (for both the input and output files) can be specified as
# relative or absolute paths. In the case of relative paths, they should be
# specified relative to the directory that contains this script.
#
# ------------------------------------------------------------------------------
#
# To run the script, set the R working directory to the directory that contains
# this script and type:
#
# source('extract_licor_data_for_gm.R')

library(PhotoGEA)

###                                                                   ###
### COMPONENTS THAT MIGHT NEED TO CHANGE EACH TIME THIS SCRIPT IS RUN ###
###                                                                   ###

# Specify the Licor data files to process and the name of the output file. There
# are two options for doing this: either the filenames can be defined directly
# as strings and vectors of strings, or they can be defined interactively via
# dialog boxes (only available on MS Windows).
CHOOSE_FILES_INTERACTIVELY <- TRUE

FILES_TO_PROCESS <- c() # Initialize the input file list
OUTPUT_FILENAME <- ""   # Initialize the output file name

# Specify the filenames depending on the value of the CHOOSE_FILES_INTERACTIVELY
# boolean
{
    if (CHOOSE_FILES_INTERACTIVELY) {
        FILES_TO_PROCESS <- choose_input_licor_files()
        OUTPUT_FILENAME <- choose_output_licor_file()
    }
    else {
        FILES_TO_PROCESS <- c(
            "2021-04-07-site 11 vulcan cs 36627-1-17.xlsx",
            "20210407-pluto-site13-36627-WT-3.xlsx"
        )
        OUTPUT_FILENAME <- "gm_combo_file.xlsx"
    }
}

###                                                                        ###
### COMPONENTS THAT ARE LESS LIKELY TO CHANGE EACH TIME THIS SCRIPT IS RUN ###
###                                                                        ###

# Specify the variables to extract. Note that when the file is loaded, any
# Unicode characters such as Greek letters will be converted into `ASCII`
# versions, e.g. the character Δ will be become `Delta`. The conversion rules
# are defined in the `UNICODE_REPLACEMENTS` data frame (see `read_licor.R`).
# Note that some of the variables are not present in the original Licor file.
# These variables will be added according to the `VARIABLES_TO_ADD` data frame
# defined above.
VARIABLES_TO_EXTRACT <- c(
    "obs",
    "hhmmss",
    "Oxygen",       # added as a column defined by a formula
    "O2",           # added as a column defined by a formula
    "[CO2]",        # added as a column defined by a formula
    "E",
    "A",
    "Ca",
    "Ci",
    "Pci",
    "Pca",
    "gsw",
    "gbw",
    "gtw",
    "gtc",
    "TleafCnd",
    "SVPleaf",
    "RHcham",
    "VPcham",
    "SVPcham",
    "VPDleaf",
    "Qin",
    "S",
    "K",
    "CO2_s",
    "Cs_licor",     # added as a column defined by a formula
    "CO2_r",
    "Ce_licor",     # added as a column defined by a formula
    "H2O_s",
    "H2O_r",
    "Flow",
    "Pa",
    "DeltaPcham",   # the name of this column is modified from ΔPcham
    "Tair",
    "Tleaf",
    "Flow_s",
    "Flow_r",
    "ppO2",         # added as a column defined by a formula
    "gsc",          # added as a column defined by a formula
    "gbc",          # added as a column defined by a formula
    "Csurface"      # added as a column defined by a formula
)

PREAMBLE_DATA_ROWS <- c(3, 5, 7, 9, 11, 13)
VARIABLE_CATEGORY_ROW <- 14
VARIABLE_NAME_ROW <- 15
VARIABLE_UNIT_ROW <- 16
DATA_START_ROW <- 17

VARIABLE_INFO_OFFSET <- 5
MAIN_DATA_OFFSET <- 1

###                                                         ###
### FUNCTIONS THAT WILL BE CALLED WHEN THIS SCRIPT RUNS     ###
### (THEY MAY REQUIRE MODIFICATIONS IF ANY FORMULAS CHANGE) ###
###                                                         ###

# Define a function that adds some Excel formulas to some columns in the Licor
# data
add_licor_formulas <- function(
    licor_file,
    variable_info_offset,
    main_data_offset
)
{
    # Get the row numbers to use
    rows <- 1:length(licor_file[['main_data']][,1]) - 1
    rows <- rows + licor_file[['data_start_row']] +
        variable_info_offset + main_data_offset

    # Fill in the "Oxygen" column
    licor_file <- write_excel_formula(
        licor_file,
        "Oxygen",
        "C$3"
    )

    # Fill in the "O2" column
    licor_file <- write_excel_formula(
        licor_file,
        "O2",
        paste0("(C", rows, "/100)*AF", rows)
    )

    # Fill in the "[CO2]" column
    licor_file <- write_excel_formula(
        licor_file,
        "[CO2]",
        paste0("Y", rows)
    )

    # Fill in the "Cs_licor" column
    licor_file <- write_excel_formula(
        licor_file,
        "Cs_licor",
        paste0("(1000000*Y", rows, ")/(1000000-AC", rows, "*1000)")
    )

    # Fill in the "Ce_licor" column
    licor_file <- write_excel_formula(
        licor_file,
        "Ce_licor",
        paste0("(1000000*AA", rows, ")/(1000000-AD", rows, "*1000)")
    )

    # Fill in the "ppO2" column
    licor_file <- write_excel_formula(
        licor_file,
        "ppO2",
        paste0("(C", rows, "/100)*AF", rows, "/100")
    )

    # Fill in the "gsc" column
    licor_file <- write_excel_formula(
        licor_file,
        "gsc",
        paste0("1/(1.6/L", rows, ")")
    )

    # Fill in the "gbc" column
    licor_file <- write_excel_formula(
        licor_file,
        "gbc",
        paste0("1/(1.37/M", rows, ")")
    )

    # Fill in the "Csurface" column
    licor_file <- write_excel_formula(
        licor_file,
        "Csurface",
        paste0(
            "((AN", rows, "-F", rows, "/2)*Y", rows, "-G", rows, ")/(AN", rows,
            "+F", rows, "/2)"
        )
    )

    return(licor_file)
}

# Calls `add_licor_formulas` for multiple Licor files
batch_add_licor_formulas <- function(
    licor_files,
    variable_info_offset,
    main_data_offset
)
{
    for (i in 1:length(licor_files)) {
        licor_files[[i]] <- add_licor_formulas(
            licor_files[[i]],
            variable_info_offset,
            main_data_offset
        )
    }
    return(licor_files)
}

# Define a function that loads Licor data from multiple files, processes it, and
# saves the results to an Excel file
print_all <- function(
    files_to_process,
    unicode_replacements,
    preamble_data_rows,
    variable_category_row,
    variable_name_row,
    variable_unit_row,
    data_start_row,
    variables_to_extract,
    variable_info_offset,
    main_data_offset,
    output_filename
)
{
    # Load the files
    licor_files <- batch_read_licor_file(
        files_to_process,
        preamble_data_rows,
        variable_category_row,
        variable_name_row,
        variable_unit_row,
        data_start_row,
        'time'
    )

    # Add blank columns to each file
    licor_files <- batch_specify_variables(
        licor_files,
        c("in",          "Oxygen",    "%"),
        c("in",          "O2",        "kPa"),
        c("in",          "[CO2]",     "micromol mol^(-1)"),
        c("CO2Scorr",    "Cs_licor",  "micromol mol^(-1)"),
        c("CO2Rcorr",    "Ce_licor",  "micromol mol^(-1)"),
        c("calculated",  "ppO2",      "bar"),
        c("calculated",  "gsc",       "mol m^(-2) s^(-1)"),
        c("calculated",  "gbc",       "mol m^(-2) s^(-1)"),
        c("calculated",  "Csurface",  "micromol mol^(-1)")
    )

    # Extract the desired columns from each file
    licor_files <- batch_extract_variables(
        licor_files,
        variables_to_extract
    )

    # Add formulas to some columns
    licor_files <- batch_add_licor_formulas(
        licor_files,
        variable_info_offset,
        main_data_offset
    )

    # Make a workbook object
    wb <- openxlsx::createWorkbook()

    # Add the sheets to the workbook object
    for (i in seq_along(files_to_process)) {
        add_licor_sheet(
            wb,
            i,
            sheet_name_from_file_name(files_to_process[i]),
            licor_files[[i]],
            variable_info_offset,
            main_data_offset
        )
    }

    # Create a file based on the workbook object
    openxlsx::saveWorkbook(wb, output_filename, overwrite = TRUE)
}

###                                                                   ###
### COMMAND THAT ACTUALLY CALLS THE FUNCTIONS WITH APPROPRIATE INPUTS ###
### (SHOULD NOT REQUIRE ANY MODIFICATIONS TO USE THE SCRIPT)          ###
###                                                                   ###

print_all(
    FILES_TO_PROCESS,
    UNICODE_REPLACEMENTS,
    PREAMBLE_DATA_ROWS,
    VARIABLE_CATEGORY_ROW,
    VARIABLE_NAME_ROW,
    VARIABLE_UNIT_ROW,
    DATA_START_ROW,
    VARIABLES_TO_EXTRACT,
    VARIABLE_INFO_OFFSET,
    MAIN_DATA_OFFSET,
    OUTPUT_FILENAME
)
