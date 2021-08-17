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
batch_specify_variables <- function(exdf_objs, ...) {
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

# combine_exdf: a function for combining the information from multiple exdf
# objects into a single list. Here, any "extra information" specified when by
# the objects will be lost (i.e., only the main data, units, and categories will
# be retained).
#
# ------------------------------------------------------------------------------
#
# INPUTS:
#
# - exdf_objs: a list of exdf objects
#
# ------------------------------------------------------------------------------
#
# OUTPUT:
#
# a single exdf object with all the data from the original objects
#
combine_exdf <- function(exdf_objs)
{
    # Make sure at least one file is defined
    if (length(exdf_objs) < 1) {
        stop("The input list is required to contain at least one exdf object")
    }

    # Get the info from the first file to use as a reference
    first_exdf <- exdf_objs[[1]]

    # Check to make sure all the files have the same variables, units, and
    # categories
    for (i in seq_along(exdf_objs)) {
        current_exdf <- exdf_objs[[i]]

        if (!identical(
            colnames(first_exdf[['main_data']]),
            colnames(current_exdf[['main_data']])
        ))
        {
            msg <- paste0(
                "The column names specified in exdf object created from '",
                current_exdf[['file_name']],
                "' do not agree with the column names specified ",
                "in the exdf object created from '",
                first_exdf[['file_name']],
                "', so the two files cannot be combined"
            )
            stop(msg)
        }

        if (!identical(first_exdf[['categories']], current_exdf[['categories']])) {
            msg <- paste0(
                "The categories specified in exdf object created from '",
                current_exdf[['file_name']],
                "' do not agree with the categories specified ",
                "in the exdf object created from '",
                first_exdf[['file_name']],
                "', so the two files cannot be combined"
            )
            stop(msg)
        }

        if (!identical(first_exdf[['units']], current_exdf[['units']])) {
            msg <- paste0(
                "The units specified in exdf object created from '",
                current_exdf[['file_name']],
                "' do not agree with the units specified ",
                "in the exdf object created from '",
                first_exdf[['file_name']],
                "', so the two files cannot be combined"
            )
            stop(msg)
        }
    }

    return(do.call("rbind.exdf", exdf_objs))
}
