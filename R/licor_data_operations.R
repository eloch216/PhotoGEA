# This script includes functions for manipulating Licor data files that have
# already been loaded; for example, adding blank columns, extracting columns,
# and combining files. In general, these functions are intended to operate on
# R lists that have been created using the functions in `read_licor.R`.

# add_licor_variables: a function for adding new variables (as blank columns) to
# an R list the representing data from a Licor file.
#
# In keeping with the other variables present in the original Licor data, each
# new variable has a category, name, and unit.
#
# ------------------------------------------------------------------------------
#
# INPUTS:
#
# - licor_file: a list representing the data from a Licor file (typically
#               produced by a a call to the `read_licor_file` defined in
#               `read_licor.R`)
#
# - variables_to_add: a data frame representing a set of new variables to add.
#                     Must have three columns named "category", "name", and "units".
#
# ------------------------------------------------------------------------------
#
# OUTPUT:
#
# a list describing Licor data having the same structure as the `licor_file`
# input
#
add_licor_variables <- function(
    licor_file,
    variables_to_add
)
{
    # Make a helping function for adding a blank column to the data and the
    # variable descriptors
    add_blank <- function(variable_info)
    {
        category <- variable_info[['category']]
        name <- variable_info[['name']]
        units <- variable_info[['units']]

        licor_file[['main_data']][[name]] <<- NA
        licor_file[['categories']][[name]] <<- category
        licor_file[['units']][[name]] <<- units
    }

    # Add the required blank columns
    for (i in 1:length(variables_to_add[,1])) {
        add_blank(variables_to_add[i,])
    }

    return(licor_file)
}

# batch_add_licor_variables: a function for adding new variables to R lists
# representing multiple Licor files
#
# ------------------------------------------------------------------------------
#
# INPUTS:
#
# - licor_files: a list of unnamed lists, each representing the information from
#                a single Licor file (typically produced by a call to the
#                `batch_read_licor_file` function defined in `read_licor.R`)
#
# All other inputs are identical to those in the `add_licor_variables` function.
#
# ------------------------------------------------------------------------------
#
# OUTPUT:
#
# a list with unnamed elements, each of which is a list describing the contents
# of a Licor file (as if it had been created by the `read_licor` function)
#
batch_add_licor_variables <- function(
    licor_files,
    variables_to_add
)
{
    lapply(
        licor_files,
        function(licor_file) {
            add_licor_variables(
                    licor_file,
                    variables_to_add
            )
        }
    )
}

# extract_licor_variables: a function for choosing a subset of Licor variables
# from an R list representing a Licor file
#
# In keeping with the other variables present in the original Licor data, each
# new variable has a category, name, and unit.
#
# ------------------------------------------------------------------------------
#
# INPUTS:
#
# - licor_file: a list representing the data from a Licor file (typically
#               produced by a a call to the `read_licor_file` defined in
#               `read_licor.R`)
#
# - variables_to_extract: a vector of variable names that should be extracted
#                         from the Licor data
#
# ------------------------------------------------------------------------------
#
# OUTPUT:
#
# a list describing Licor data having the same structure as the `licor_file`
# input
#
extract_licor_variables <- function(
    licor_file,
    variables_to_extract
)
{
    # Only attempt to extract variables that are included in the Licor file
    variables_to_extract <-
        variables_to_extract[variables_to_extract %in% colnames(licor_file[['main_data']])]

    licor_file[['categories']] <-
        licor_file[['categories']][,variables_to_extract]

    licor_file[['units']] <-
        licor_file[['units']][,variables_to_extract]

    licor_file[['main_data']] <-
        licor_file[['main_data']][,variables_to_extract]

    return(licor_file)
}

# batch_extract_licor_variables: a function for extracting variables from R
# lists representing multiple Licor files
#
# ------------------------------------------------------------------------------
#
# INPUTS:
#
# - licor_files: a list of unnamed lists, each representing the information from
#                a single Licor file (typically produced by a call to the
#                `batch_read_licor_file` function defined in `read_licor.R`)
#
# All other inputs are identical to those in the `extract_licor_variables`
# function.
#
# ------------------------------------------------------------------------------
#
# OUTPUT:
#
# a list with unnamed elements, each of which is a list describing the contents
# of a Licor file (as if it had been created by the `read_licor` function)
#
batch_extract_licor_variables <- function(
    licor_files,
    variables_to_extract
)
{
    lapply(
        licor_files,
        function(licor_file) {
            extract_licor_variables(
                    licor_file,
                    variables_to_extract
            )
        }
    )
}

# combine_licor_files: a function for combining the information from multiple
# Licor files into a single list. Here, only the filenames, categories, units,
# and data are retained (i.e., any parameters specified when reading the files
# will be lost).
#
# ------------------------------------------------------------------------------
#
# INPUTS:
#
# - licor_files: a list of unnamed lists, each representing the information from
#                a single Licor file (typically produced by a call to the
#                `batch_read_licor_file` function defined in `read_licor.R`)
#
# ------------------------------------------------------------------------------
#
# OUTPUT:
#
# a list describing Licor data having the same structure as the output of a call
# to the `read_licor_file` function.
#
combine_licor_files <- function(
    licor_files
)
{
    # Make sure at least one file is defined
    if (length(licor_files) < 1) {
        stop("The input list is required to contain at least one Licor file")
    }

    # Get the info from the first file to use as a reference
    first_file <- licor_files[[1]]

    # Set up the combo list
    initial_data_frame <- data.frame(
        matrix(
            ncol = ncol(first_file[['main_data']]),
            nrow = 0
        )
    )

    colnames(initial_data_frame) <- colnames(first_file[['main_data']])

    combo_info <- list(
        file_name = character(0),
        categories = first_file[['categories']],
        units = first_file[['units']],
        main_data = initial_data_frame
    )

    # Get the remaining info
    for (i in seq_along(licor_files)) {
        # Get the current file
        current_file <- licor_files[[i]]

        # Check for possible errors
        if (!identical(
            colnames(first_file[['main_data']]),
            colnames(current_file[['main_data']])
        ))
        {
            msg <- paste0(
                "The column names specified in Licor file '",
                current_file[['file_name']],
                "' do not agree with the column names specified in '",
                first_file[['file_name']],
                "', so the two files cannot be combined"
            )
            stop(msg)
        }

        if (!identical(first_file[['categories']], current_file[['categories']])) {
            msg <- paste0(
                "The categories specified in Licor file '",
                current_file[['file_name']],
                "' do not agree with the categories specified in '",
                first_file[['file_name']],
                "', so the two files cannot be combined"
            )
            stop(msg)
        }

        if (!identical(first_file[['units']], current_file[['units']])) {
            msg <- paste0(
                "The units specified in Licor file '",
                current_file[['file_name']],
                "' do not agree with the units specified in '",
                first_file[['file_name']],
                "', so the two files cannot be combined"
            )
            stop(msg)
        }

        # Add this file to the combined info
        combo_info[['main_data']] <- rbind(
            combo_info[['main_data']],
            current_file[['main_data']]
        )

        combo_info[['file_name']] <- c(
            combo_info[['file_name']],
            current_file[['file_name']]
        )
    }

    return(combo_info)
}
