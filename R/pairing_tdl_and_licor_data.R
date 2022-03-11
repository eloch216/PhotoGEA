# Specify the respiration rate (in units of micromol / m^2 / s) to store in the
# Licor file. (The rate can be supplied as a positive or negative number, and
# its absolute value will be stored.)
specify_respiration <- function(licor_exdf, respiration) {
    # Add a new column to the Licor file in preparation for adding the
    # respiration information
    licor_exdf <- specify_variables(
        licor_exdf,
        c("in", "respiration", "micromol m^(-2) s^(-1)")
    )

    # Store it in the Licor file and return the updated file
    licor_exdf[,'respiration'] <- abs(respiration)

    return(licor_exdf)
}

batch_specify_respiration <- function(licor_exdfs, respiration) {
    lapply(licor_exdfs, function(x) {specify_respiration(x, respiration)})
}

# Specify the oxygen level (as a percentage, typically either 2 or 21) to store
# in the Licor file
specify_oxygen <- function(licor_exdf, oxygen) {
    # Add a new column to the Licor file in preparation for adding the
    # oxygen information
    licor_exdf <- specify_variables(
        licor_exdf,
        c("in", "Oxygen", "%")
    )

    # Store it in the Licor file and return the updated file
    licor_exdf[,'Oxygen'] <- oxygen

    return(licor_exdf)
}

batch_specify_oxygen <- function(licor_exdfs, oxygen) {
    lapply(licor_exdfs, function(x) {specify_oxygen(x, oxygen)})
}

get_oxygen_info_from_preamble <- function(licor_exdf) {
    # Add a new column to the Licor file in preparation for adding the
    # oxygen information
    licor_exdf <- specify_variables(
        licor_exdf,
        c("in", "Oxygen", "%")
    )

    # Try to get the oxygen information from the Licor file's preamble
    oxygen <- c()
    preamble <- licor_exdf[['preamble']]
    for (pe in preamble) {
        if ("Oxygen" %in% colnames(pe)) {
            # Convert the value to numeric form, if possible, and stop going
            # through the preamble
            oxygen <- try_as_numeric(pe[['Oxygen']])
            break
        }
    }

    # Make sure we actually got the info
    if (length(oxygen) == 0) {
        msg <- paste0(
            "Could not automatically get oxygen information from Licor file:\n'",
            licor_exdf[['file_name']],
            "'\nConsider adding it with the `specify_oxygen` function rather ",
            "than using `get_oxygen_info_from_preamble`"
        )
        stop(msg)
    }

    # Store it in the Licor file and return the updated file
    licor_exdf[,'Oxygen'] <- oxygen

    return(licor_exdf)
}

batch_get_oxygen_info_from_preamble <- function(licor_exdfs) {
    lapply(licor_exdfs, function(x) {get_oxygen_info_from_preamble(x)})
}

get_genotype_info_from_licor_filename <- function(licor_exdf) {
    # Add some new columns to the Licor file in preparation for adding the plant
    # information
    licor_exdf <- specify_variables(
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
        e <- plant_specification2[3]
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

pair_licor_and_tdl <- function(
    licor_file,
    tdl_data,
    licor_timestamp_column_name,
    tdl_timestamp_column_name,
    max_allowed_time_difference  # minutes
)
{
    # Add some new columns to the Licor file in preparation for adding the TDL
    # data
    licor_file <- specify_variables(
        licor_file,
        c("calibrated TDL",              "cycle_num",              ""),
        c("calibrated TDL (sample)",     "valve_number_s",         ""),
        c("calibrated TDL (sample)",     "tdl_time_s",             ""),
        c("calibrated TDL (sample)",     "calibrated_12c_s",       "ppm"),
        c("calibrated TDL (sample)",     "calibrated_13c_s",       "ppm"),
        c("calibrated TDL (sample)",     "total_mixing_ratio_s",   "ppm"),
        c("calibrated TDL (sample)",     "total_isotope_ratio_s",  "ppt"),
        c("calibrated TDL (reference)",  "valve_number_r",         ""),
        c("calibrated TDL (reference)",  "tdl_time_r",             ""),
        c("calibrated TDL (reference)",  "calibrated_12c_r",       "ppm"),
        c("calibrated TDL (reference)",  "calibrated_13c_r",       "ppm"),
        c("calibrated TDL (reference)",  "total_mixing_ratio_r",   "ppm"),
        c("calibrated TDL (reference)",  "total_isotope_ratio_r",  "ppt")
    )

    # Make sure the time columns have the correct class, or we won't be able to
    # store timestamps in them properly
    licor_file[['main_data']][['tdl_time_s']] <-
        as.POSIXlt(licor_file[['main_data']][['tdl_time_s']])

    licor_file[['main_data']][['tdl_time_r']] <-
        as.POSIXlt(licor_file[['main_data']][['tdl_time_r']])

    # Next, we need to find the "site number" of the Licor data, which is
    # encoded in the file name. This will tell us which valves to use for the
    # sample and reference CO2 data from the TDL file.
    trimmed_name <-
        tools::file_path_sans_ext(basename(licor_file[['file_name']]))

    # Search for the following pattern: one or more spaces, followed by 'site',
    # possibly followed by one or more space, followed by one or more digits,
    # followed by one or more spaces.
    site_phrase <-
        regmatches(trimmed_name, regexpr(" +site *[0-9]+ +", trimmed_name))

    # Make sure we found something
    if (length(site_phrase) == 0) {
        msg <- paste0(
            "Could not extract TDL site number from Licor file:\n'",
            licor_file[['file_name']],
            "'\nThe filename must include the phrase ` siteNN `, where `NN` ",
            "is the TDL valve number that should be used for the sample data"
        )
        stop(msg)
    }

    # Remove the whitespace and the word "site"
    site_phrase <- sub(" ", "", site_phrase)
    site_phrase <- sub("site", "", site_phrase)

    # Make sure the result is treated as a numeric value
    sample_valve <- as.numeric(site_phrase)

    # The reference valve is one less than the sample valve
    reference_valve <- sample_valve - 1

    # Store this info in the Licor file
    licor_file[['main_data']][['valve_number_s']] <- sample_valve
    licor_file[['main_data']][['valve_number_r']] <- reference_valve

    # Now extract TDL data for each time point
    for (i in seq_len(nrow(licor_file[['main_data']]))) {
        # Find the TDL cycle that contains the closest time to the Licor data
        # point time
        licor_time <- licor_file[,licor_timestamp_column_name][i]
        time_differences <- abs(difftime(
            tdl_data[[tdl_timestamp_column_name]],
            licor_time,
            units = 'min'
        ))
        min_time_difference <- min(time_differences)
        indx_of_closest_tdl_pnt <- match(min_time_difference, time_differences)
        cycle_of_closest_tdl_pnt <- tdl_data[['cycle_num']][indx_of_closest_tdl_pnt]

        if (min_time_difference > max_allowed_time_difference) {
            msg <- paste0(
                "Could not find a time point in the TDL data corresponding to ",
                format(licor_time),
                ", a point in the Licor data.\nThe nearest time in the TDL data was ",
                format(tdl_data[[tdl_timestamp_column_name]][indx_of_closest_tdl_pnt]),
                ",\nwhich is ", min_time_difference, " minutes away, exceeding ",
                "the maximum allowed value of ", max_allowed_time_difference,
                " minutes"
            )
            stop(msg)
        }

        # Find the TDL times for the sample and reference measurements
        tdl_time_sample <-
            tdl_data[tdl_data[['cycle_num']] == cycle_of_closest_tdl_pnt &
                tdl_data[['valve_number']] == sample_valve, tdl_timestamp_column_name]

        tdl_time_reference <-
            tdl_data[tdl_data[['cycle_num']] == cycle_of_closest_tdl_pnt &
                tdl_data[['valve_number']] == reference_valve, tdl_timestamp_column_name]

        # Find the sample CO2 data
        calibrated_12c_s <-
            tdl_data[tdl_data[['cycle_num']] == cycle_of_closest_tdl_pnt &
                tdl_data[['valve_number']] == sample_valve, 'calibrated_12c']

        calibrated_13c_s <-
            tdl_data[tdl_data[['cycle_num']] == cycle_of_closest_tdl_pnt &
                tdl_data[['valve_number']] == sample_valve, 'calibrated_13c']

        total_mixing_ratio_s <-
            tdl_data[tdl_data[['cycle_num']] == cycle_of_closest_tdl_pnt &
                tdl_data[['valve_number']] == sample_valve, 'total_mixing_ratio']

        total_isotope_ratio_s <-
            tdl_data[tdl_data[['cycle_num']] == cycle_of_closest_tdl_pnt &
                tdl_data[['valve_number']] == sample_valve, 'total_isotope_ratio']

        # Find the reference CO2 data
        calibrated_12c_r <-
            tdl_data[tdl_data[['cycle_num']] == cycle_of_closest_tdl_pnt &
                tdl_data[['valve_number']] == reference_valve, 'calibrated_12c']

        calibrated_13c_r <-
            tdl_data[tdl_data[['cycle_num']] == cycle_of_closest_tdl_pnt &
                tdl_data[['valve_number']] == reference_valve, 'calibrated_13c']

        total_mixing_ratio_r <-
            tdl_data[tdl_data[['cycle_num']] == cycle_of_closest_tdl_pnt &
                tdl_data[['valve_number']] == reference_valve, 'total_mixing_ratio']

        total_isotope_ratio_r <-
            tdl_data[tdl_data[['cycle_num']] == cycle_of_closest_tdl_pnt &
                tdl_data[['valve_number']] == reference_valve, 'total_isotope_ratio']

        # Store results in the Licor file
        licor_file[['main_data']][['cycle_num']][i] <- cycle_of_closest_tdl_pnt
        licor_file[['main_data']][['tdl_time_s']][i] <- tdl_time_sample
        licor_file[['main_data']][['tdl_time_r']][i] <- tdl_time_reference

        licor_file[['main_data']][['calibrated_12c_s']][i] <- calibrated_12c_s
        licor_file[['main_data']][['calibrated_13c_s']][i] <- calibrated_13c_s
        licor_file[['main_data']][['total_mixing_ratio_s']][i] <- total_mixing_ratio_s
        licor_file[['main_data']][['total_isotope_ratio_s']][i] <- total_isotope_ratio_s

        licor_file[['main_data']][['calibrated_12c_r']][i] <- calibrated_12c_r
        licor_file[['main_data']][['calibrated_13c_r']][i] <- calibrated_13c_r
        licor_file[['main_data']][['total_mixing_ratio_r']][i] <- total_mixing_ratio_r
        licor_file[['main_data']][['total_isotope_ratio_r']][i] <- total_isotope_ratio_r
    }

    return(licor_file)
}

batch_pair_licor_and_tdl <- function(
    licor_files,
    tdl_data,
    licor_timestamp_column_name,
    tdl_timestamp_column_name,
    max_allowed_time_difference
)
{
    lapply(
        licor_files,
        function(licor_file) {
            pair_licor_and_tdl(
                licor_file,
                tdl_data,
                licor_timestamp_column_name,
                tdl_timestamp_column_name,
                max_allowed_time_difference
            )
        }
    )
}
