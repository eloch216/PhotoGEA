get_oxygen_info_from_preamble <- function(licor_exdf) {
    # Add a new column to the Licor file in preparation for adding the
    # oxygen information
    licor_exdf <- document_variables(
        licor_exdf,
        c("in", "Oxygen", "%")
    )

    # Try to get the oxygen information from the Licor file's preamble
    preamble <- licor_exdf[['preamble']]

    oxygen <- if ('Oxygen' %in% colnames(preamble)) {
        try_as_numeric(preamble[['Oxygen']])
    } else if ('SysConst:Oxygen' %in% colnames(preamble)) {
        try_as_numeric(preamble[['SysConst:Oxygen']])
    } else {
        msg <- paste0(
            "Could not automatically get oxygen information from Licor file:\n'",
            licor_exdf[['file_name']],
            "'\nConsider adding oxygen values with the `set_variable` function ",
            "rather than using `get_oxygen_info_from_preamble`"
        )
        stop(msg)
    }

    # Store it in the Licor file and return the updated file
    licor_exdf[,'Oxygen'] <- oxygen

    return(licor_exdf)
}
