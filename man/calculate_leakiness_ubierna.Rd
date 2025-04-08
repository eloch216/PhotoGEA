\name{calculate_leakiness_ubierna}

\alias{calculate_leakiness_ubierna}

\title{Calculate leakiness}

\description{
  Calculates leakiness (\code{phi}) from combined gas exchange and isotope
  discrimination measurements as described in Ubierna et al. (2013). This
  function can accomodate alternative colum names for the variables taken from
  \code{exdf_obj}; it also checks the units of each required column and will
  produce an error if any units are incorrect.
}

\usage{
  calculate_leakiness_ubierna(
    exdf_obj,
    e = -3,
    a_bar_column_name = 'a_bar',
    a_column_name = 'A',
    ci_column_name = 'Ci',
    co2_s_column_name = 'CO2_s',
    csurface_column_name = 'Csurface',
    delta_c13_r_column_name = 'delta_C13_r',
    delta_obs_tdl_column_name = 'Delta_obs_tdl',
    rl_column_name = 'RL',
    t_column_name = 't'
  )
}

\arguments{
  \item{exdf_obj}{An \code{exdf} object.}

  \item{e}{
    The isotopic fractionation during day respiration in \code{ppt}.
  }

  \item{a_bar_column_name}{
    The name of the column in \code{exdf_obj} that contains the weighted
    isotopic fractionation across the boundary layer and stomata in \code{ppt}.
    Values of \code{a_bar} are typically calculated using
    \code{\link{calculate_ternary_correction}}.
  }

  \item{a_column_name}{
    The name of the column in \code{exdf_obj} that contains the net CO2
    assimilation rate in \code{micromol m^(-2) s^(-1)}.
  }

  \item{ci_column_name}{
    The name of the column in \code{exdf_obj} that contains the intercellular
    CO2 concentration in \code{micromol mol^(-1)}.
  }

  \item{co2_s_column_name}{
    The name of the column in \code{exdf_obj} that contains the CO2
    concentration in the sample line (outgoing air) in \code{micromol mol^(-1)}.
  }

  \item{csurface_column_name}{
    The name of the column in \code{exdf_obj} that contains the CO2
    concentration at the leaf surface in \code{micromol mol^(-1)}. Values of
    \code{Csurface} are typically calculated using
    \code{\link{calculate_gas_properties}}.
  }

  \item{delta_c13_r_column_name}{
    The name of the column in \code{exdf_obj} that contains the CO2 isotope
    ratio in the reference line (incoming air) in \code{ppt}.
  }

  \item{delta_obs_tdl_column_name}{
    The name of the column in \code{exdf_obj} that contains the observed isotope
    discrimination values in \code{ppt}.
  }

  \item{rl_column_name}{
    The name of the column in \code{exdf_obj} that contains the rate of day
    respiration in \code{micromol m^(-2) s^(-1)}.
  }

  \item{t_column_name}{
    The name of the column in \code{exdf_obj} that contains the ternary
    correction factor (\code{dimensionless}). Values of \code{t} are typically
    calculated using \code{\link{calculate_ternary_correction}}
  }
}

\details{
  This function uses the model for photosynthetic discrimination against 13C in
  C4 plants to determine leakiness values, as described in Ubierna et al.
  (2013). In particular, the following equations from that source are
  implemented in the code:

  \itemize{
    \item Isotopic fractionation due to day respiration (\code{e_prime}) is
          calculated using Equation 21.
    \item Leakiness including respiratory and photorespiratory fractionations
          under high light (\code{phi_i}) is calculated using Equation 16.
    \item Leakiness including respiratory and photorespiratory fractionations
          and Cs under high light (\code{phi_is}) is calculated using Equation
          15.
    \item Leakiness ignoring respiratory and photorespiratory fractionations and
          Cs (\code{phi_sim}) is calculated using Equation 17.
  }

  References:

  Ubierna, N., Sun, W., Kramer, D. M. and Cousins, A. B. "The efficiency of C4
  photosynthesis under low light conditions in Zea mays, Miscanthus x giganteus
  and Flaveria bidentis." Plant, Cell & Environment 36, 365–381 (2013)
  [\doi{10.1111/j.1365-3040.2012.02579.x}].
}

\value{
  An \code{exdf} object based on \code{exdf_obj} that includes the following
  additional columns, calculated as described above: \code{e_prime},
  \code{phi_i}, \code{phi_is}, and \code{phi_sim}. The category for each of
  these new columns is \code{calculate_leakiness_ubierna} to indicate that they
  were created using this function.
}

\examples{
## In this example we load gas exchange and TDL data files, calibrate the TDL
## data, pair the data tables together, and then calculate leakiness. The
## results from this example are not meaningful because these measurements
## were not collected from C4 plants.

# Read the TDL data file, interpreting the times as if they were specified in
# the local time zone. (The measurements were actually made in the US Central
# time zone, but `tz = 'US/Central'` is not recognized by some operating
# systems.)
tdl_data <- read_gasex_file(
  PhotoGEA_example_file_path('tdl_for_gm.dat'),
  'TIMESTAMP',
  list(tz = Sys.timezone())
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

# Read the gas exchange data file, interpreting the times as if they were
# specified in the local time zone. (The measurements were actually made in the
# US Central time zone, but `tz = 'US/Central'` is not recognized by some
# operating systems.)
licor_data <- read_gasex_file(
  PhotoGEA_example_file_path('licor_for_gm_site11.xlsx'),
  'time',
  list(tz = Sys.timezone())
)

# Get TDL valve information from Licor file name; for this TDL system, the
# reference valve is 12 when the sample valve is 11
licor_data <- get_sample_valve_from_filename(licor_data, list('11' = 12))

# Get oxygen info from the Licor file preamble (needed for calculate_gamma_star)
licor_data <- get_oxygen_from_preamble(licor_data)

# Pair the Licor and TDL data by locating the TDL cycle corresponding to each
# Licor measurement
licor_data <- pair_gasex_and_tdl(licor_data, processed_tdl$tdl_data)

# Calculate total pressure (needed for calculate_gas_properties)
licor_data <- calculate_total_pressure(licor_data)

# Calculate Csurface (needed for calculate_ternary_correction)
licor_data <- calculate_gas_properties(licor_data)

# Calculate ternary correction
licor_data <- calculate_ternary_correction(licor_data)

# Calculate isotope discrimination (needed for calculate_leakiness_ubierna)
licor_data <- calculate_isotope_discrimination(licor_data)

# Set respiration (needed for calculate_leakiness_ubierna)
licor_data <- set_variable(
  licor_data,
  'RL',
  'micromol m^(-2) s^(-1)',
  value = 1.2
)

# Calculate leakiness
licor_data <- calculate_leakiness_ubierna(licor_data)

# View some of the results
licor_data[, c('replicate', 'CO2_s', 'Delta_obs_tdl', 'phi_i', 'phi_sim')]
}

\concept{exdf}
