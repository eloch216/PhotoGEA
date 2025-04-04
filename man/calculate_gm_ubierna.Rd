\name{calculate_gm_ubierna}

\alias{calculate_gm_ubierna}

\title{Calculate mesophyll conductance to CO2 diffusion}

\description{
  Calculates mesophyll conductance to CO2 diffusion (\code{gmc}) from combined
  gas exchange and isotope discrimination measurements as described in Ubierna
  et al. (2018). This function can accomodate alternative colum names for the
  variables taken from \code{exdf_obj}; it also checks the units of each
  required column and will produce an error if any units are incorrect.
}

\usage{
  calculate_gm_ubierna(
    exdf_obj,
    e = -3,
    f = 11,
    a_bar_column_name = 'a_bar',
    a_column_name = 'A',
    ci_column_name = 'Ci',
    co2_s_column_name = 'CO2_s',
    csurface_column_name = 'Csurface',
    delta_c13_r_column_name = 'delta_C13_r',
    delta_obs_tdl_column_name = 'Delta_obs_tdl',
    gamma_star_column_name = 'Gamma_star_tl',
    rl_column_name = 'RL',
    total_pressure_column_name = 'total_pressure',
    t_column_name = 't'
  )
}

\arguments{
  \item{exdf_obj}{An \code{exdf} object.}

  \item{e}{
    The isotopic fractionation during day respiration in \code{ppt}.
  }

  \item{f}{
    The isotopic fractionation during photorespiration in \code{ppt}.
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

  \item{gamma_star_column_name}{
    The name of the column in \code{exdf_obj} that contains the chloroplastic
    CO2 concentration at which CO2 gains from Rubisco carboxylation are exactly
    balanced by CO2 losses from Rubisco oxygenation, at leaf temperature,
    expressed in \code{micromol mol^(-1)}.  Values of \code{Gamma_star} at leaf
    temperature are typically calculated using
    \code{\link{calculate_gamma_star}} or
    \code{\link{calculate_temperature_response}}.
  }

  \item{rl_column_name}{
    The name of the column in \code{exdf_obj} that contains the rate of
    non-photorespiratory CO2 release in the light, in
    \code{micromol m^(-2) s^(-1)}.
  }

  \item{total_pressure_column_name}{
    The name of the column in \code{exdf_obj} that contains the total pressure
    in \code{bar}.
  }

  \item{t_column_name}{
    The name of the column in \code{exdf_obj} that contains the ternary
    correction factor (\code{dimensionless}). Values of \code{t} are typically
    calculated using \code{\link{calculate_ternary_correction}}
  }
}

\details{
  This function uses the comprehensive model for photosynthetic discrimination
  against 13C in C3 plants to calculate mesophyll conductance, as described in
  Ubierna et al. (2018). In particular, the following equations from that source
  are implemented in the code:

  \itemize{
    \item Isotopic fractionation due to day respiration (\code{e_prime}) is
          calculated using Equations 28 and 30.
    \item Isotopic discrimination due to photorespiration (\code{Delta_f}), due
          to day respiration (\code{Delta_e}), and that would occur if Ci = Cc
          in the absence of any respiratory fractionation (\code{Delta_i}) are
          calculated using Equations 34, 33, and 31, respectively.
    \item Mesophyll conductance to CO2 diffusion (\code{gmc}) is calculated
          using Equation 44. This equation is broken up into two factors called
          \code{Delta_difference} and \code{equation_top} which are separately
          returned in the output from \code{calculate_gm_ubierna}.
  }

  For an alternative method for calculating \code{gmc}, see
  \code{\link{calculate_gm_busch}}.

  References:

  Ubierna, N., Holloway-Phillips, M.-M. and Farquhar, G. D. "Using Stable Carbon
  Isotopes to Study C3 and C4 Photosynthesis: Models and Calculations." in
  Photosynthesis: Methods and Protocols (ed. Covshoff, S.) 155–196 (Springer,
  2018) [\doi{10.1007/978-1-4939-7786-4_10}].
}

\value{
  An \code{exdf} object based on \code{exdf_obj} that includes the following
  additional columns, calculated as described above: \code{e_prime},
  \code{Delta_i}, \code{Delta_e}, \code{Delta_f}, \code{Delta_difference},
  \code{equation_top}, and \code{gmc}. The category for each of these new
  columns is \code{calculate_gm_ubierna} to indicate that they were created
  using this function.
}

\examples{
## In this example we load gas exchange and TDL data files, calibrate the TDL
## data, pair the data tables together, and then calculate mesophyll conductance

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

# Set Rubisco specificity (needed for calculate_gamma_star)
licor_data <- set_variable(
    licor_data,
    'rubisco_specificity_tl',
    'M / M',
    value = 90
)

# Calculate Gamma_star (needed for calculate_gm_ubierna)
licor_data <- calculate_gamma_star(licor_data)

# Calculate isotope discrimination (needed for calculate_gm_ubierna)
licor_data <- calculate_isotope_discrimination(licor_data)

# Set respiration (needed for calculate_gm_ubierna)
licor_data <- set_variable(
  licor_data,
  'RL',
  'micromol m^(-2) s^(-1)',
  value = 1.2
)

# Calculate mesophyll conductance
licor_data <- calculate_gm_ubierna(licor_data)

# Calculate Cc using the new values of mesophyll conductance
licor_data <- calculate_temperature_response(
  licor_data,
  c3_temperature_param_flat['gmc_norm']
)

licor_data <- set_variable(
  licor_data,
  'gmc_at_25',
  units = licor_data$units$gmc,
  value = licor_data[, 'gmc']
)

licor_data <- apply_gm(licor_data)

# View some of the results
licor_data[, c('replicate', 'CO2_s', 'Delta_obs_tdl', 'gmc', 'Ci', 'Cc')]
}

\concept{exdf}
