get_sample_valve_from_licor_filename <- function(licor_exdf, reference_table) {
    # We need to find the "site number" of the Licor data, which is encoded in
    # the file name. This will tell us which valves to use for the sample CO2
    # data from the TDL file.
    trimmed_name <-
        tools::file_path_sans_ext(basename(licor_exdf[['file_name']]))

    # Search for the following pattern: one or more spaces, followed by 'site',
    # possibly followed by one or more space, followed by one or more digits,
    # followed by one or more spaces.
    site_phrase <-
        regmatches(trimmed_name, regexpr("site *[0-9]+", trimmed_name))

    # Make sure we found something
    if (length(site_phrase) == 0) {
        msg <- paste0(
            "Could not extract TDL site number from Licor file:\n'",
            licor_exdf[['file_name']],
            "'\nThe filename must include the phrase `siteNN`, where `NN` ",
            "is the TDL valve number that should be used for the sample data"
        )
        stop(msg)
    }

    # Remove any whitespace and the word "site"
    site_phrase <- sub(" ", "", site_phrase)
    site_phrase <- sub("site", "", site_phrase)

    # Make sure the result is treated as a numeric value
    sample_valve <- as.numeric(site_phrase)

    # Add it as a new column in the exdf
    licor_exdf <- set_variable(
        licor_exdf,
        'valve_number_s',
        "",
        "calibrated TDL (sample)",
        sample_valve
    )

    # Now add the corresponding reference valve and return the exdf
    set_variable(
        licor_exdf,
        'valve_number_r',
        '',
        'calibrated TDL (reference)',
        reference_table[[as.character(sample_valve)]]
    )
}
