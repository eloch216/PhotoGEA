# This script includes functions for manipulating Licor data files that have
# already been loaded; for example, adding blank columns, extracting columns,
# and combining files. In general, these functions are intended to operate on
# R lists that have been created using the functions in `read_licor.R`.

# specify_variables: a function for specifying the units and categories of
# columns of an exdf object; any new columns will be initialized to NA.
#
# ------------------------------------------------------------------------------
#
# INPUTS:
#
# - exdf_obj: an exdf object
#
# any additional inputs should be vectors of strings that describe the variable
# specifications, where the first element is the category, second is the name,
# and third is the units
#
# ------------------------------------------------------------------------------
#
# OUTPUT:
#
# an exdf object with new and/or updated columns
#
specify_variables <- function(exdf_obj, ...)
{
    if (!is.exdf(exdf_obj)) {
        stop("exdf_obj must be an exdf object")
    }

    variable_specs <- list(...)

    for (v in variable_specs) {
        if (length(v) != 3 & !is.character(v)) {
            stop("all variable specifications must be provided as vectors of three strings")
        }

        exdf_obj <- set_column_info(
            exdf_obj,
            v[2],
            v[3],
            v[1]
        )
    }

    return(exdf_obj)
}

# batch_specify_variables: a function for specifying the units and categories of
# columns of multiple exdf objects; any new columns will be initialized to NA.
#
# ------------------------------------------------------------------------------
#
# INPUTS:
#
# - exdf_objs: a list of exdf objects
#
# any additional inputs should be vectors of strings that describe the variable
# specifications, where the first element is the category, second is the name,
# and third is the units
#
# ------------------------------------------------------------------------------
#
# OUTPUT:
#
# a list of exdf objects with new and/or updated columns
#
batch_specify_variables <- function(exdf_objs, ...)
{
    lapply(exdf_objs, function(x) {specify_variables(x, ...)})
}

# extract_variables: a function for choosing a subset of variables from an exdf
# object
#
# ------------------------------------------------------------------------------
#
# INPUTS:
#
# - exdf_obj: an exdf object
#
# - variables_to_extract: a vector of variable names that should be extracted
#
# ------------------------------------------------------------------------------
#
# OUTPUT:
#
# an exdf object with a subset of the original columns
#
extract_variables <- function(
    exdf_obj,
    variables_to_extract
)
{
    # Check for variables that are not included in the Licor file
    not_in_file <-
        variables_to_extract[!variables_to_extract %in% colnames(exdf_obj)]

    if (length(not_in_file) > 0) {
        msg <- paste0(
            "the following variables are not included in the exdf object: ",
            paste(not_in_file, collapse = ", ")
        )
        stop(msg)
    }

    exdf_obj[['main_data']] <- exdf_obj[['main_data']][,variables_to_extract]
    exdf_obj[['units']] <- exdf_obj[['units']][,variables_to_extract]
    exdf_obj[['categories']] <- exdf_obj[['categories']][,variables_to_extract]

    return(exdf_obj)
}

# batch_extract_variables: a function for choosing a subset of variables from
# multiple exdf objects
#
# ------------------------------------------------------------------------------
#
# INPUTS:
#
# - exdf_objs: a list of exdf objects
#
# - variables_to_extract: a vector of variable names that should be extracted
#
# ------------------------------------------------------------------------------
#
# OUTPUT:
#
# a list of exdf objects, each having a subset of its orignal columns
#
batch_extract_variables <- function(exdf_objs, variables_to_extract) {
    lapply(exdf_objs, function(x) {extract_variables(x, variables_to_extract)})
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

    # Check to make sure all the files have the same variables, units, and
    # categories
    for (i in seq_along(licor_files)) {
        current_file <- licor_files[[i]]

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
    }

    return(do.call("rbind.exdf", licor_files))
}
