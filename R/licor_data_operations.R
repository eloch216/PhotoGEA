# This script includes functions for manipulating Licor data files that have
# already been loaded; for example, adding blank columns, extracting columns,
# and combining files. In general, these functions are intended to operate on
# R lists that have been created using the functions in `read_licor.R`.

# add_licor_variables: a function for adding new variables (as blank columns) to
# an R list the representing data from a Licor file.
#
# In keeping with the other variables present in the original Licor data, each
# new variable has a type, name, and unit.
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
#                     Must have three columns named "type", "name", and "units".
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
        type <- variable_info[['type']]
        name <- variable_info[['name']]
        units <- variable_info[['units']]

        licor_file[['main_data']][[name]] <<- NA
        licor_file[['types']][[name]] <<- type
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
# new variable has a type, name, and unit.
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
    licor_file[['types']] <-
        licor_file[['types']][,variables_to_extract]

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
# Licor files into a single list.
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
    # Get the info from the first file
    combo_info <- licor_files[[1]]

    # Get the remaining info
    for (i in 2:length(licor_files)) {
        combo_info[['main_data']] <- rbind(
            combo_info[['main_data']],
            licor_files[[i]][['main_data']]
        )
    }

    return(combo_info)
}
