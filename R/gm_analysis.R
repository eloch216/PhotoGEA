source("read_tdl.R")
source("tdl_data_operations.R")
source("tdl_calculations.R")
source("read_licor.R")
source("pairing_tdl_and_licor_data.R")
source("gm_calculations.R")
source("calculate_cc.R")
source("basic_stats.R")
source("save_file.R")

# Define constants that will determine the behavior of some functions in this
# script

PERFORM_CALCULATIONS <- TRUE

SAVE_RESULTS <- TRUE

RESPIRATION <- -0.710568448235977

MIN_GM <- 0.1
MAX_GM <- 3.0

MIN_CC <- 0.0

# Names of important columns in the TDL data
TDL_TIMESTAMP_COLUMN_NAME <- 'TIMESTAMP'
TDL_VALVE_COLUMN_NAME <- 'valve_number'
TDL_RAW_12C_COLUMN_NAME <- 'Conc12C_Avg'
TDL_RAW_13C_COLUMN_NAME <- 'Conc13C_Avg'

# Specify the variables to extract from the TDL data files. Note that when the
# files are loaded, any Unicode characters such as Greek letters will be
# converted into `ASCII` versions, e.g. the character Δ will be become `Delta`.
# The conversion rules are defined in the `UNICODE_REPLACEMENTS` data frame
# (see `read_licor.R`).
TDL_COLUMNS_TO_EXTRACT <- c(
    TDL_TIMESTAMP_COLUMN_NAME,
    TDL_VALVE_COLUMN_NAME,
    TDL_RAW_12C_COLUMN_NAME,
    TDL_RAW_13C_COLUMN_NAME
)

# Names of important columns in the Licor data
GM_COLUMN_NAME <- "gmc"
CI_COLUMN_NAME <- "Ci"
CC_COLUMN_NAME <- "Cc"
A_COLUMN_NAME <- "A"

# Specify the variables to extract from the Licor data files. Note that when the
# files are loaded, any Unicode characters such as Greek letters will be
# converted into `ASCII` versions, e.g. the character Δ will be become `Delta`.
# The conversion rules are defined in the `UNICODE_REPLACEMENTS` data frame
# (see `read_licor.R`).
LICOR_VARIABLES_TO_EXTRACT <- c(
    "obs",
    "time",
    "E",
    "A",
    "Ca",
    "Ci",
    "Pci",
    "Pca",
    "gsw",
    "gbw",
    "gtw",
    "gtc",
    "TleafCnd",
    "SVPleaf",
    "RHcham",
    "VPcham",
    "SVPcham",
    "VPDleaf",
    "Qin",
    "S",
    "K",
    "CO2_s",
    "CO2_r",
    "H2O_s",
    "H2O_r",
    "Flow",
    "Pa",
    "DeltaPcham",   # the name of this column is modified from ΔPcham
    "Tair",
    "Tleaf",
    "Flow_s",
    "Flow_r"
)

# Specify variables to analyze, i.e., variables where the average, standard
# deviation, and standard error will be determined for each event across the
# different reps. Note that when the file is loaded, any Unicode characters
# such as Greek letters will be converted into `ASCII` versions, e.g. the
# character Δ will be become `Delta`. The conversion rules are defined in the
# `UNICODE_REPLACEMENTS` data frame (see `read_licor.R`).
VARIABLES_TO_ANALYZE <- c(
    A_COLUMN_NAME,
    CI_COLUMN_NAME,
    CC_COLUMN_NAME,
    GM_COLUMN_NAME,
    "gsw",
    "Ci-Cc"
)

if (PERFORM_CALCULATIONS) {
    # Get all the TDL information and process it

    tdl_files <- batch_read_tdl_file(
        choose_input_tdl_files(),
        rows_to_skip = 1,
        variable_name_row = 2,
        variable_unit_row = 3,
        data_start_row = 5,
        timestamp_colname = TDL_TIMESTAMP_COLUMN_NAME
    )

    tdl_files <- batch_extract_licor_variables(
        tdl_files,
        c(
            TDL_TIMESTAMP_COLUMN_NAME,
            TDL_VALVE_COLUMN_NAME,
            'Conc12C_Avg',
            'Conc13C_Avg'
        )
    )

    tdl_files <- combine_tdl_files(tdl_files)

    tdl_files <- identify_tdl_cycles(
        tdl_files,
        valve_column_name = TDL_VALVE_COLUMN_NAME,
        cycle_start_valve = 20,
        expected_cycle_length_minutes = 2.7,
        expected_cycle_num_pts = 9
    )

    processed_tdl_data <- process_tdl_cycles(
        tdl_files[['main_data']],
        valve_column_name = TDL_VALVE_COLUMN_NAME,
        noaa_valve = 2,
        calibration_0_valve = 20,
        calibration_1_valve = 21,
        calibration_2_valve = 23,
        calibration_3_valve = 26,
        raw_12c_colname = TDL_RAW_12C_COLUMN_NAME,
        raw_13c_colname = TDL_RAW_13C_COLUMN_NAME,
        noaa_cylinder_co2_concentration = 294.996,  # ppm
        noaa_cylinder_isotope_ratio = -8.40,        # ppt
        calibration_isotope_ratio = -11.505,        # ppt
        f_other = 0.00474,                          # fraction of CO2 that is not 13C16O16O or 12C16O16O
        R_VPDB = 0.0111797
    )

    # Get all the Licor information and process it

    licor_files <- batch_read_licor_file(
        choose_input_licor_files(),
        PREAMBLE_DATA_ROWS,
        VARIABLE_TYPE_ROW,
        VARIABLE_NAME_ROW,
        VARIABLE_UNIT_ROW,
        DATA_START_ROW
    )

    licor_files <- batch_extract_licor_variables(
        licor_files,
        LICOR_VARIABLES_TO_EXTRACT
    )

    licor_files <- batch_get_genotype_info_from_licor_filename(licor_files)

    licor_files <- batch_get_oxygen_info_from_preamble(licor_files)

    licor_files <- batch_specify_respiration(licor_files, RESPIRATION)

    licor_files <- batch_pair_licor_and_tdl(
        licor_files,
        processed_tdl_data[['tdl_data']]
    )

    licor_files <- combine_licor_files(licor_files)

    licor_files <- calculate_gm(licor_files)

    licor_files <- calculate_cc(
        licor_files,
        CC_COLUMN_NAME,
        CI_COLUMN_NAME,
        A_COLUMN_NAME,
        GM_COLUMN_NAME
    )

    # Make a copy of the data where we have removed extreme gm and Cc values
    licor_files_no_outliers <- licor_files
    licor_files_no_outliers[['main_data']] <-
        licor_files_no_outliers[['main_data']][which(
            licor_files_no_outliers[['main_data']][['gmc']] > MIN_GM &
            licor_files_no_outliers[['main_data']][['gmc']] < MAX_GM &
            licor_files_no_outliers[['main_data']][['Cc']] > MIN_CC),]

    # Get stats for each event by averaging over all corresponding reps
    event_stats <- basic_stats(
        licor_files_no_outliers[['main_data']],
        'genotype',
        'event',
        VARIABLES_TO_ANALYZE,
        'sa'
    )

    # Get stats for each rep by averaging over all corresponding observations
    rep_stats <- basic_stats(
        licor_files_no_outliers[['main_data']],
        'genotype',
        'event_replicate',
        VARIABLES_TO_ANALYZE,
        'sa'
    )
}

if (SAVE_RESULTS) {
    base_dir <- getwd()
    if (interactive() & .Platform$OS.type == "windows") {
        base_dir <- choose.dir(caption="Select folder for output files")
    }

    write.csv(processed_tdl_data[['tdl_data']], file.path(base_dir, "tdl_data_processed.csv"), row.names=FALSE)
    write.csv(processed_tdl_data[['calibration_zero']], file.path(base_dir, "tdl_calibration_zero.csv"), row.names=FALSE)
    write.csv(processed_tdl_data[['calibration_12CO2']], file.path(base_dir, "tdl_calibration_12CO2.csv"), row.names=FALSE)
    write.csv(processed_tdl_data[['calibration_13CO2_data']], file.path(base_dir, "tdl_calibration_13CO2_data.csv"), row.names=FALSE)
    write.csv(processed_tdl_data[['calibration_13CO2_fit']], file.path(base_dir, "tdl_calibration_13CO2_fit.csv"), row.names=FALSE)

    save_licor_file(licor_files, file.path(base_dir, "gm_calculations_outliers_included.csv"))
    save_licor_file(licor_files_no_outliers, file.path(base_dir, "gm_calculations_outliers_excluded.csv"))
    write.csv(event_stats, file.path(base_dir, "gm_stats_by_event_outliers_excluded.csv"), row.names=FALSE)
    write.csv(rep_stats, file.path(base_dir, "gm_stats_by_rep_outliers_excluded.csv"), row.names=FALSE)
}
