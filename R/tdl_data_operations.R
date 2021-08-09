
# combine_tdl_files: a function for combining the information from multiple
# TDL files into a single list. Here, only the filenames, units, and data are
# retained (i.e., any parameters specified when reading the files will be lost).
#
# ------------------------------------------------------------------------------
#
# INPUTS:
#
# - tdl_files: a list of unnamed lists, each representing the information from
#       a single TDL file (typically produced by a call to the
#       `batch_read_tdl_file` function defined in `read_tdl.R`)
#
# ------------------------------------------------------------------------------
#
# OUTPUT:
#
# a list describing TDL data having the same structure as the output of a call
# to the `read_tdl_file` function.
#
combine_tdl_files <- function(
    tdl_files
)
{
    # Make sure at least one file is defined
    if (length(tdl_files) < 1) {
        stop("The input list is required to contain at least one TDL file")
    }

    # Get the info from the first file to use as a reference
    first_file <- tdl_files[[1]]

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
        units = first_file[['units']],
        main_data = initial_data_frame
    )

    # Get the remaining info
    for (i in seq_along(tdl_files)) {
        # Get the current file
        current_file <- tdl_files[[i]]

        # Check for possible errors
        if (!identical(
            colnames(first_file[['main_data']]),
            colnames(current_file[['main_data']])
        ))
        {
            msg <- paste0(
                "The column names specified in TDL file '",
                current_file[['file_name']],
                "' do not agree with the column names specified in '",
                first_file[['file_name']],
                "', so the two files cannot be combined"
            )
            stop(msg)
        }

        if (!identical(first_file[['units']], current_file[['units']])) {
            msg <- paste0(
                "The units specified in TDL file '",
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
