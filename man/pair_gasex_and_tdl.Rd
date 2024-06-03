\name{pair_gasex_and_tdl}

\alias{pair_gasex_and_tdl}

\title{Pair gas exchange and TDL data}

\description{
  Identifies the closest TDL cycle corresponding to each entry in the gas
  exchange data and adds the TDL data to the gas exchange data.
}

\usage{
  pair_gasex_and_tdl(
    gasex_exdf,
    tdl_exdf,
    max_allowed_time_difference = 1,
    gasex_timestamp_column_name = 'time',
    tdl_timestamp_column_name = 'TIMESTAMP'
  )
}

\arguments{
  \item{gasex_exdf}{
    An \code{exdf} object representing data from a photosynthetic gas exchange
    measurement system.
  }

  \item{tdl_exdf}{
    An \code{exdf} object representing calibrated data from a tunable diode
    laser absorption spectroscopy system. Typically this is the output from
    applying \code{\link{process_tdl_cycle_erml}} or
    \code{\link{process_tdl_cycle_polynomial}} to a set of uncalibrated TDL data.
  }

  \item{max_allowed_time_difference}{
    The maximum time difference (in minutes) to allow between gas exchange and
    TDL timestamp values.
  }

  \item{gasex_timestamp_column_name}{
    The name of the column in \code{gasex_exdf} that contains the timestamp
    values.
  }

  \item{tdl_timestamp_column_name}{
    The name of the column in \code{tdl_exdf} that contains the timestamp
    values.
  }
}

\details{
  When making combined gas exchange and isotope discrimination measurements
  using a portable photosynthetic gas exchange system (such as a Licor LI-6800)
  coupled with a tunable diode laser (TDL) absorption spectroscopy system, the
  TDL's gas handling system cycles through several gas lines (or sites) by
  opening and closing valves. When analyzing such data, a key step is to
  combine TDL and gas exchange data that were measured at the same times.

  The \code{pair_gasex_and_tdl} function performs this operation by locating the
  TDL cycle whose timestamp is closest to each Licor file entry. Then, the
  12C, 13C, total CO2, and delta_13C values measured by the TDL from the Licor's
  sample and reference lines during that cycle are added to the gas exchange
  data as new columns.
}

\value{
  An \code{exdf} object based on \code{gasex_exdf} that includes TDL values
  measured at the same times as the original gas exchange logs. Several new
  columns are added: \code{'cycle_num'}, \code{'tdl_time_s'},
  \code{'calibrated_12c_s'}, \code{'calibrated_13c_s'}, \code{'total_CO2_s'},
  \code{'delta_C13_s'}, \code{'tdl_time_r'}, \code{'calibrated_12c_r'},
  \code{'calibrated_13c_r'}, \code{'total_CO2_r'}, and \code{'delta_C13_r'}.
  Variables with \code{'_s'} in the name refer to TDL measurements from the
  Licor sample line, and \code{'_r'} indicates the reference line. The category
  of each new column is \code{pair_gasex_and_tdl} to indicate that it was
  created using this function.
}

\examples{
## In this example we load gas exchange and TDL data files, calibrate the TDL
## data, and pair the data tables together

# Read the TDL data file, making sure to interpret the time zone as US Central
# time
tdl_data <- read_gasex_file(
  PhotoGEA_example_file_path('tdl_for_gm.dat'),
  'TIMESTAMP',
  list(tz = 'US/Central')
)

# Identify cycles within the TDL data
tdl_data <- identify_tdl_cycles(
  tdl_data,
  valve_column_name = 'valve_number',
  cycle_start_valve = 20,
  expected_cycle_length_minutes = 2.7,
  expected_cycle_num_valves = 9,
  timestamp_colname = 'TIMESTAMP'
)

# Use reference tanks to calibrate the TDL data
processed_tdl <- consolidate(by(
  tdl_data,
  tdl_data[, 'cycle_num'],
  process_tdl_cycle_erml,
  noaa_valve = 2,
  calibration_0_valve = 20,
  calibration_1_valve = 21,
  calibration_2_valve = 23,
  calibration_3_valve = 26,
  noaa_cylinder_co2_concentration = 294.996,
  noaa_cylinder_isotope_ratio = -8.40,
  calibration_isotope_ratio = -11.505
))

# Read the gas exchange data, making sure to interpret the time stamp in the US
# Central time zone
licor_data <- read_gasex_file(
  PhotoGEA_example_file_path('licor_for_gm_site11.xlsx'),
  'time',
  list(tz = 'US/Central')
)

# Get TDL valve information from Licor file name; for this TDL system, the
# reference valve is 12 when the sample valve is 11
licor_data <- get_sample_valve_from_filename(licor_data, list('11' = 12))

# Pair the Licor and TDL data by locating the TDL cycle corresponding to each
# Licor measurement
licor_data <- pair_gasex_and_tdl(licor_data, processed_tdl$tdl_data)

# View some of the results
licor_data[, c('A', 'delta_C13_r', 'delta_C13_s', 'total_CO2_r', 'total_CO2_s')]
}


\concept{exdf}
