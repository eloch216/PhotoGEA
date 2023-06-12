get_genotype_info_from_licor_filename <- function(licor_exdf) {
    # Add some new columns to the Licor file in preparation for adding the plant
    # information
    licor_exdf <- document_variables(
        licor_exdf,
        c("plant specification", "genotype",                 "NA"),
        c("plant specification", "event",                    "NA"),
        c("plant specification", "replicate",                "NA"),
        c("plant specification", "genotype_event",           "NA"),
        c("plant specification", "event_replicate",          "NA"),
        c("plant specification", "genotype_event_replicate", "NA"),
        c("plant specification", "original_file",            "NA")
    )

    # Get the filename without the path
    name <- basename(licor_exdf[['file_name']])

    # There are two naming conventions we need to check

    # Search for the following pattern: one space, followed by one or
    # more alphanumeric characters, followed by a dash, followed by one or more
    # alphanumeric characters, followed by a dash, followed by one or more
    # alphanumeric characters, followed by a period. Essentially, we expect the
    # filename to end with ' XXX-YYY-ZZZ.xlsx', where `XXX` is the genotype,
    # `YYY` is the event, and `ZZZ` is the replicate.
    plant_specification1 <-
        regmatches(name, regexpr(" [[:alnum:]]+-[[:alnum:]]+-[[:alnum:]]+\\.xlsx", name))

    # Search for the following pattern: one space, followed by one or
    # more alphanumeric characters, followed by a space, followed by one or more
    # alphanumeric characters, followed by a space, followed by one or more
    # alphanumeric characters, followed by a period. Essentially, we expect the
    # filename to end with ' XXX YYY ZZZ.xlsx', where `XXX` is the genotype,
    # `YYY` is the event, and `ZZZ` is the replicate.
    plant_specification2 <-
        regmatches(name, regexpr(" [[:alnum:]]+ [[:alnum:]]+ [[:alnum:]]+\\.xlsx", name))

    # Make sure we found something
    if (length(plant_specification1) == 0 & length(plant_specification2) == 0) {
        msg <- paste0(
            "Could not extract plant specification information from Licor file:\n'",
            licor_exdf[['file_name']],
            "'\nThe filename must end with ` GGG-EEE-RRR.xlsx` or ` GGG EEE RRR.xlsx`, ",
            "where `GGG`, `EEE`, and `RRR` are alphanumeric specifiers for ",
            "the genotype, event, and replicate represented by the file"

        )
        stop(msg)
    }

    # Extract the info
    g <- character(0)
    e <- character(0)
    r <- character(0)
    if (length(plant_specification1) > 0) {
        # Remove the period, the extension, and the whitespace
        plant_specification1 <- sub(" ", "", plant_specification1)
        plant_specification1 <- sub("\\.xlsx", "", plant_specification1)

        # Split the specification by the dashes
        plant_specification1 <- strsplit(plant_specification1, "-")[[1]]

        g <- plant_specification1[1]
        e <- fix_wt(plant_specification1[2])
        r <- plant_specification1[3]
    } else {
        # Remove the period and the extension
        plant_specification2 <- sub("\\.xlsx", "", plant_specification2)

        # Split the specification by the spaces
        plant_specification2 <- strsplit(plant_specification2, " ")[[1]]

        g <- plant_specification2[2]
        e <- fix_wt(plant_specification2[3])
        r <- plant_specification2[4]
    }

    # Store the info in the file and return it
    licor_exdf[,'genotype'] <- g
    licor_exdf[,'event'] <- e
    licor_exdf[,'replicate'] <- r
    licor_exdf[,'genotype_event'] <- paste0(g, "-", e)
    licor_exdf[,'event_replicate'] <- paste0(e, "-", r)
    licor_exdf[,'genotype_event_replicate'] <- paste0(g, "-", e, "-", r)
    licor_exdf[,'original_file'] <- licor_exdf[['file_name']]
    return(licor_exdf)
}

batch_get_genotype_info_from_licor_filename <- function(licor_exdfs) {
    lapply(licor_exdfs, function(x) {get_genotype_info_from_licor_filename(x)})
}
