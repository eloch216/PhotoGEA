\name{calculate_isotope_discrimination}

\alias{calculate_isotope_discrimination}

\title{Calculate photosynthetic isotope discrimination}

\description{
  Calculates photosynthetic carbon isotope discrimination from combined gas
  exchange and tunable diode laser absorption spectroscopy measurements.
}

\usage{
  calculate_isotope_discrimination(
    exdf_obj,
    co2_r_column_name = 'CO2_r',
    co2_s_column_name = 'CO2_s',
    delta_C13_r_column_name = 'delta_C13_r',
    delta_C13_s_column_name = 'delta_C13_s',
    h2o_r_column_name = 'H2O_r',
    h2o_s_column_name = 'H2O_s',
    tdl_12C_r_column_name = 'calibrated_12c_r',
    tdl_12C_s_column_name = 'calibrated_12c_s'
  )
}

\arguments{
  \item{exdf_obj}{
    An \code{exdf} object representing combined data from a gas exchange +
    isotope discrimination measurement system. Typically \code{exdf_obj} is
    produced by calling \code{\link{pair_gasex_and_tdl}}.
  }

  \item{co2_r_column_name}{
    The name of the column in \code{exdf_obj} that contains the CO2
    concentration in the gas exchange reference line (incoming air) as measured
    by the gas exchange system in \code{micromol mol^(-1)}.
  }

  \item{co2_s_column_name}{
    The name of the column in \code{exdf_obj} that contains the CO2
    concentration in the gas exchange sample line (outgoing air) in
    \code{micromol mol^(-1)}.
  }

  \item{delta_C13_r_column_name}{
    The name of the column in \code{exdf_obj} that contains the CO2 isotope
    ratio in the gas exchange reference line (incoming air) in \code{ppt}.
  }

  \item{delta_C13_s_column_name}{
    The name of the column in \code{exdf_obj} that contains the CO2 isotope
    ratio in the gas exchange sample line (outgoing air) in \code{ppt}.
  }

  \item{h2o_r_column_name}{
    The name of the column in \code{exdf_obj} that contains the H2O
    concentration in the gas exchange reference line (incoming air) as measured
    by the gas exchange system in \code{mmol mol^(-1)}.
  }

  \item{h2o_s_column_name}{
    The name of the column in \code{exdf_obj} that contains the H2O
    concentration in the gas exchange sample line (outgoing air) as measured
    by the gas exchange system in \code{mmol mol^(-1)}.
  }

  \item{tdl_12C_r_column_name}{
    The name of the column in \code{exdf_obj} that contains the 12CO2
    concentration in the gas exchange reference line (incoming air) as measured
    by the TDL in \code{ppm}.
  }

  \item{tdl_12C_s_column_name}{
    The name of the column in \code{exdf_obj} that contains the 12CO2
    concentration in the gas exchange sample line (outgoing air) as measured
    by the TDL in \code{ppm}.
  }
}

\details{
  As described in Ubierna et al. (2018), photosynthetic 13C discrimination can
  be determined from combined gas exchange and tunable diode laser (TDL)
  absorption spectroscopy measurements according to:

  \code{Delta_obs = xsi * (delta_out - delta_in) / (1 + delta_out - xsi * (delta_out - delta_in))},

  where \code{Delta_obs} is the observed discrimination, \code{delta_in} and
  \code{delta_out} are the carbon isotope ratios in dry air flowing in and out
  of the leaf chamber. \code{xsi} is given by

  \code{xsi = C_in / (C_in - C_out)},

  where \code{C_in} and \code{C_out} are the mole fractions of 12CO2 in dry air
  flowing in and out of the leaf chamber. (See equations 5 and 6 in Ubierna et
  al. (2018)).

  In practice, there are multiple options for calculating \code{Delta_obs} and
  \code{xsi} because CO2 concentrations are measured by both the gas exchange
  system and the TDL. For example, we can alternately calculate \code{xsi} as
  \code{xsi_tdl = C_in_tdl / (C_in_tdl - C_out_tdl)} or
  \code{xsi_gasex = C_in_gasex / (C_in_gasex - C_out_gasex)}. Likewise, we can
  also calculate \code{Delta_obs_tdl} using \code{xsi_tdl} or
  \code{Delta_obs_gasex} using \code{xsi_gasex}. The TDL values are typically
  preferred in subsequent calculations, but it can be useful to compare the two
  different versions as a consistency check; the TDL and gas exchange values
  should be similar to each other.

  There are two subtelties associated with \code{xsi_gasex}. One is that the gas
  exchange system generally measures the total CO2 concentration, not just the
  12CO2 concentration. Typically there is much less 13CO2 than 12CO2 so this is
  usually not a large source of error.

  The other issue is that the gas exchange system generally measures CO2
  concentrations in wet air. Thus, it is important to use "corrected" values of
  CO2 concentrations that account for the "dilution effect" due to water vapor
  in the air. This effect is described in the Licor LI-6400 manual: "This is a
  correction we don’t do, at least when computing CO2 concentration in the
  LI-6400. The dilution effect is simply this: as you add molecules of a gas
  (water vapor, for example) to a mixture, the fraction of that mixture that is
  made up of something else (mole fraction of CO2, for instance) has to
  decrease, since the total number of molecules in the mixture has increased.
  Now for an airsteam flowing though a chamber containing a transpiring leaf (or
  in a chamber sitting on moist soil), there very definitely is dilution.
  However, we ignore that effect when computing CO2 concentration, but account
  for it when computing photosynthetic rate (or soil CO2 efflux). Thus, the
  LI-6400 IRGA is always indicating the actual CO2 concentration, not what the
  CO2 concentration would be if there were no water vapor in it."

  To account for the dilution effect, we define a "corrected" CO2 concentration
  as \code{CO2_corrected = CO2 / (1 - H2O)}, where \code{H2O} is the water vapor
  concentration in the air. Note: the TDL always measures concentrations in dry
  air, so no correction is required.

  References:

  Ubierna, N., Holloway-Phillips, M.-M. and Farquhar, G. D. "Using Stable Carbon
  Isotopes to Study C3 and C4 Photosynthesis: Models and Calculations." in
  Photosynthesis: Methods and Protocols (ed. Covshoff, S.) 155–196 (Springer,
  2018) [\doi{10.1007/978-1-4939-7786-4_10}].
}

\value{
  An \code{exdf} object based on \code{exdf_obj} that includes several new
  columns: \code{CO2_r_corrected}, \code{CO2_s_corrected},
  \code{Delta_obs_gasex}, \code{Delta_obs_tdl}, \code{xsi_gasex}, and
  \code{xsi_tdl}.
}

\examples{
## In this example we load gas exchange and TDL data files, calibrate the TDL
## data, pair the data tables together, and then calculate isotope
## discrimination

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

# Calculate isotope discrimination
licor_data <- calculate_isotope_discrimination(licor_data)

# View some of the results
licor_data[, c('A', 'xsi_gasex', 'xsi_tdl', 'Delta_obs_gasex', 'Delta_obs_tdl')]
}

\concept{exdf}
