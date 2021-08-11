# Specify the respiration rate (in units of micromol / m^2 / s) to store in the
# Licor file. (The rate can be supplied as a positive or negative number, and
# its absolute value will be stored.)
specify_respiration <- function(licor_file, respiration) {
    # Add some new columns to the Licor file in preparation for adding the
    # respiration information
    variables_to_add <- data.frame(
        rbind(
            c("in", "respiration", "micromol m^(-2) s^(-1)")
        ),
        stringsAsFactors = FALSE
    )
    colnames(variables_to_add) <- c("type", "name", "units")

    licor_file <- add_licor_variables(licor_file, variables_to_add)

    # Store it in the Licor file and return the updated file
    licor_file[['main_data']][['respiration']] <- abs(respiration)

    return(licor_file)
}

batch_specify_respiration <- function(licor_files, respiration) {
    lapply(
        licor_files,
        function(licor_file) {
            specify_respiration(licor_file, respiration)
        }
    )
}

# Specify the oxygen level (as a percentage, typically either 2 or 21) to store
# in the Licor file
specify_oxygen <- function(licor_file, oxygen) {
    # Add some new columns to the Licor file in preparation for adding the
    # oxygen information
    variables_to_add <- data.frame(
        rbind(
            c("in", "Oxygen", "%")
        ),
        stringsAsFactors = FALSE
    )
    colnames(variables_to_add) <- c("type", "name", "units")

    licor_file <- add_licor_variables(licor_file, variables_to_add)

    # Store it in the Licor file and return the updated file
    licor_file[['main_data']][['Oxygen']] <- oxygen

    return(licor_file)
}

batch_specify_oxygen <- function(licor_files, oxygen) {
    lapply(
        licor_files,
        function(licor_file) {
            specify_oxygen(licor_file, oxygen)
        }
    )
}

get_oxygen_info_from_preamble <- function(licor_file) {
    # Add some new columns to the Licor file in preparation for adding the
    # oxygen information
    variables_to_add <- data.frame(
        rbind(
            c("in", "Oxygen", "%")
        ),
        stringsAsFactors = FALSE
    )
    colnames(variables_to_add) <- c("type", "name", "units")

    licor_file <- add_licor_variables(licor_file, variables_to_add)

    # Try to get the oxygen information from the Licor file's preamble
    oxygen <- c()
    preamble <- licor_file[['preamble']]
    for (pe in preamble) {
        if ("Oxygen" %in% colnames(pe)) {
            oxygen <- pe[['Oxygen']]
        }
    }

    # Make sure we actually got the info
    if (length(oxygen) == 0) {
        msg <- paste0(
            "Could not automatically get oxygen information from Licor file:\n'",
            licor_file[['file_name']],
            "'\nConsider adding it with the `specify_oxygen` function rather ",
            "than using `get_oxygen_info_from_preamble`"
        )
        stop(msg)
    }

    # Store it in the Licor file and return the updated file
    licor_file[['main_data']][['Oxygen']] <- oxygen

    return(licor_file)
}

batch_get_oxygen_info_from_preamble <- function(licor_files) {
    lapply(
        licor_files,
        function(licor_file) {
            get_oxygen_info_from_preamble(licor_file)
        }
    )
}

get_genotype_info_from_licor_filename <- function(licor_file) {
    # Add some new columns to the Licor file in preparation for adding the plant
    # information
    variables_to_add <- data.frame(
        rbind(
            c("plant specification",  "genotype",       ""),
            c("plant specification",  "event",          ""),
            c("plant specification",  "replicate",      ""),
            c("plant specification",  "event_rep",      ""),
            c("plant specification",  "original_file",  "")
        ),
        stringsAsFactors = FALSE
    )
    colnames(variables_to_add) <- c("type", "name", "units")

    licor_file <- add_licor_variables(licor_file, variables_to_add)

    # Get the filename without the path
    name <- basename(licor_file[['file_name']])

    # Search for the following pattern: one or more spaces, followed by one or
    # more alphanumeric characters, followed by a dash, followed by one or more
    # alphanumeric characters, followed by a dash, followed by one or more
    # alphanumeric characters, followed by a period. Essentially, we expect the
    # filename to end with ' XXX-YYY-ZZZ.xlsx', where `XXX` is the genotype,
    # `YYY` is the event, and `ZZZ` is the replicate.
    plant_specification <-
        regmatches(name, regexpr(" +[[:alnum:]]+-[[:alnum:]]+-[[:alnum:]]+\\.", name))

    # Make sure we found something
    if (length(plant_specification) == 0) {
        msg <- paste0(
            "Could not extract plant specification information from Licor file:\n'",
            licor_file[['file_name']],
            "'\nThe filename must end with `GGG-EEE-RRR`, where `GGG`, `EEE`, ",
            "and `RRR` are alphanumeric specifiers for the genotype, event, ",
            "and replicate represented by the file"

        )
        stop(msg)
    }

    # Remove the whitespace and the period
    plant_specification <- sub(" ", "", plant_specification)
    plant_specification <- sub("\\.", "", plant_specification)

    # Split the specification by the dashes
    plant_specification <- strsplit(plant_specification, "-")[[1]]

    g <- plant_specification[1]
    e <- plant_specification[2]
    r <- plant_specification[3]

    # Store the info in the file and return it
    licor_file[['main_data']][['genotype']] <- g
    licor_file[['main_data']][['event']] <- e
    licor_file[['main_data']][['replicate']] <- r
    licor_file[['main_data']][['event_rep']] <- paste0(e, "-", r)
    licor_file[['main_data']][['original_file']] <- licor_file[['file_name']]
    return(licor_file)
}

batch_get_genotype_info_from_licor_filename <- function(licor_files) {
    lapply(
        licor_files,
        function(licor_file) {
            get_genotype_info_from_licor_filename(licor_file)
        }
    )
}

pair_licor_and_tdl <- function(
    licor_file,
    tdl_data
)
{
    # Add some new columns to the Licor file in preparation for adding the TDL
    # data
    variables_to_add <- data.frame(
        rbind(
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
        ),
        stringsAsFactors = FALSE
    )
    colnames(variables_to_add) <- c("type", "name", "units")

    licor_file <- add_licor_variables(licor_file, variables_to_add)

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
    # followed by one or more digits, followed by one or more spaces.
    site_phrase <-
        regmatches(trimmed_name, regexpr(" +site[0-9]+ +", trimmed_name))

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
        licor_time <- licor_file[['main_data']][['time']][i]
        time_differences <- abs(difftime(
            tdl_data[['TIMESTAMP']],
            licor_time,
            units = 'sec'
        ))
        min_time_difference <- min(time_differences)
        indx_of_closest_tdl_pnt <- match(min_time_difference, time_differences)
        cycle_of_closest_tdl_pnt <- tdl_data[['cycle_num']][indx_of_closest_tdl_pnt]

        # Find the TDL times for the sample and reference measurements
        tdl_time_sample <-
            tdl_data[which(tdl_data[['cycle_num']] == cycle_of_closest_tdl_pnt &
                tdl_data[['valve_number']] == sample_valve), 'TIMESTAMP']

        tdl_time_reference <-
            tdl_data[which(tdl_data[['cycle_num']] == cycle_of_closest_tdl_pnt &
                tdl_data[['valve_number']] == reference_valve), 'TIMESTAMP']

        # Find the sample CO2 data
        calibrated_12c_s <-
            tdl_data[which(tdl_data[['cycle_num']] == cycle_of_closest_tdl_pnt &
                tdl_data[['valve_number']] == sample_valve), 'calibrated_12c']

        calibrated_13c_s <-
            tdl_data[which(tdl_data[['cycle_num']] == cycle_of_closest_tdl_pnt &
                tdl_data[['valve_number']] == sample_valve), 'calibrated_13c']

        total_mixing_ratio_s <-
            tdl_data[which(tdl_data[['cycle_num']] == cycle_of_closest_tdl_pnt &
                tdl_data[['valve_number']] == sample_valve), 'total_mixing_ratio']

        total_isotope_ratio_s <-
            tdl_data[which(tdl_data[['cycle_num']] == cycle_of_closest_tdl_pnt &
                tdl_data[['valve_number']] == sample_valve), 'total_isotope_ratio']

        # Find the reference CO2 data
        calibrated_12c_r <-
            tdl_data[which(tdl_data[['cycle_num']] == cycle_of_closest_tdl_pnt &
                tdl_data[['valve_number']] == reference_valve), 'calibrated_12c']

        calibrated_13c_r <-
            tdl_data[which(tdl_data[['cycle_num']] == cycle_of_closest_tdl_pnt &
                tdl_data[['valve_number']] == reference_valve), 'calibrated_13c']

        total_mixing_ratio_r <-
            tdl_data[which(tdl_data[['cycle_num']] == cycle_of_closest_tdl_pnt &
                tdl_data[['valve_number']] == reference_valve), 'total_mixing_ratio']

        total_isotope_ratio_r <-
            tdl_data[which(tdl_data[['cycle_num']] == cycle_of_closest_tdl_pnt &
                tdl_data[['valve_number']] == reference_valve), 'total_isotope_ratio']

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
    tdl_data
)
{
    lapply(
        licor_files,
        function(licor_file) {
            pair_licor_and_tdl(
                licor_file,
                tdl_data
            )
        }
    )
}
