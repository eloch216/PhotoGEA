read_licor_file <- function(
    file_name,
    timestamp_colname = NA,
    instrument_type = 'Licor LI-6800',
    file_type = 'AUTO'
)
{
    # Try to determine the file type from its name, if necessary
    file_type <- if (file_type == 'AUTO') {
        extension <- tools::file_ext(file_name)
        if (extension == '') {
            'plaintext'
        } else if (extension == 'xlsx') {
            'Excel'
        } else {
            stop(paste('Cannot automatically determine file_type for file:', file_name))
        }
    } else {
        file_type
    }

    # Try to load the file using the appropriate method
    licor_exdf <- if (instrument_type == 'Licor LI-6800' && file_type == 'plaintext') {
        read_licor_6800_plaintext(file_name)
    } else if (instrument_type == 'Licor LI-6800' && file_type == 'Excel') {
        read_licor_6800_Excel(file_name)
    } else {
        stop(paste('Unsupported (instrument_type file_type) option:', instrument_type, file_type))
    }

    # Make sure the timestamp column is properly interpreted
    if (!is.na(timestamp_colname)) {
        licor_exdf$main_data[[timestamp_colname]] <- as.POSIXlt(
            licor_exdf$main_data[[timestamp_colname]],
            origin = "1970-01-01"
        )
    }
    licor_exdf$timestamp_colname <- timestamp_colname

    return(licor_exdf)
}
