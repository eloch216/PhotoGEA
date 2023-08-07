get_sample_valve_from_filename <- function(
    exdf_obj,
    reference_table = NULL
)
{
    # Make sure the filename is included in the exdf object
    if (!'file_name' %in% names(exdf_obj)) {
        stop('exdf_obj$file_name must be defined; see read_gasex_file for more info')
    }

    # We need to find the "site number" of the Licor data, which is encoded in
    # the file name. This will tell us which valves to use for the sample CO2
    # data from the TDL file.
    trimmed_name <-
        tools::file_path_sans_ext(basename(exdf_obj[['file_name']]))

    # Search for the following pattern: one or more spaces, followed by 'site',
    # possibly followed by one or more space, followed by one or more digits,
    # followed by one or more spaces.
    site_phrase <-
        regmatches(trimmed_name, regexpr('site *[0-9]+', trimmed_name))

    # Make sure we found something
    if (length(site_phrase) == 0) {
        msg <- paste0(
            "Could not extract TDL site number from Licor file:\n'",
            exdf_obj[['file_name']],
            "'\nThe filename must include the phrase `site NN` or `siteNN`, ",
            "where `NN` is the TDL valve number that should be used for the ",
            "sample line"
        )
        stop(msg)
    }

    # Remove any whitespace and the word "site"
    site_phrase <- sub(' ', '', site_phrase)
    site_phrase <- sub('site', '', site_phrase)

    # Make sure the result is treated as a numeric value
    sample_valve <- as.numeric(site_phrase)

    # Add it as a new column in the exdf
    exdf_obj <- set_variable(
        exdf_obj,
        'valve_number_s',
        '',
        'get_sample_valve_from_filename',
        sample_valve
    )

    # Add the corresponding reference valve if necessary and return the exdf
    if (is.null(reference_table)) {
        exdf_obj
    } else {
        set_variable(
            exdf_obj,
            'valve_number_r',
            '',
            'get_sample_valve_from_filename',
            reference_table[[as.character(sample_valve)]]
        )
    }
}
