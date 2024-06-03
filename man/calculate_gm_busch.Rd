\name{calculate_gm_busch}

\alias{calculate_gm_busch}

\title{Calculate mesophyll conductance to CO2 diffusion}

\description{
  Calculates mesophyll conductance to CO2 diffusion (\code{gmc}) from combined
  gas exchange and isotope discrimination measurements as described in Busch
  et al. (2020). This function can accomodate alternative colum names for the
  variables taken from \code{exdf_obj}; it also checks the units of each
  required column and will produce an error if any units are incorrect.
}

\usage{
  calculate_gm_busch(
    exdf_obj,
    e = -3,
    f = 11,
    e_star_equation = 20,
    gm_type = 'dis',
    a_bar_column_name = 'a_bar',
    a_column_name = 'A',
    ci_column_name = 'Ci',
    co2_s_column_name = 'CO2_s',
    csurface_column_name = 'Csurface',
    delta_c13_r_column_name = 'delta_C13_r',
    delta_obs_growth_column_name = 'Delta_obs_growth',
    delta_obs_tdl_column_name = 'Delta_obs_tdl',
    gamma_star_column_name = 'Gamma_star',
    rd_column_name = 'Rd',
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

  \item{e_star_equation}{
    The equation from Busch et al. (2020) to use for calculating \code{e_star};
    must be \code{19} or \code{20}.
  }

  \item{gm_type}{
    Determines whether day respiration is assumed to be isotopically connected
    to the CBB cycle (\code{gm_type = 'con'}) or isotopically disconnected from
    the CBB cycle (\code{gm_type = 'dis'}). This choice will determine which
    equations are used to calculate mesophyll conductance; when \code{gm_type}
    is \code{'con'}, Equations 2 and 21 will be used; otherwise, Equations 13
    and 22 will be used.
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

  \item{delta_obs_growth_column_name}{
    The name of the column in \code{exdf_obj} that contains the observed
    discrimination under the typical CO2 concentration in the plant's
    environment during its growth (in \code{ppt}). This is only required when
    using Equation 20 for \code{e_star} (see \code{e_star_equation}).
  }

  \item{delta_obs_tdl_column_name}{
    The name of the column in \code{exdf_obj} that contains the observed isotope
    discrimination values in \code{ppt}.
  }

  \item{gamma_star_column_name}{
    The name of the column in \code{exdf_obj} that contains the CO2 compensation
    point in the absence of day respiration in \code{micromol mol^(-1)}.  Values
    of \code{Gamma_star} are typically calculated using
    \code{\link{calculate_gamma_star}} or \code{\link{calculate_arrhenius}}.
  }

  \item{rd_column_name}{
    The name of the column in \code{exdf_obj} that contains the rate of day
    respiration in \code{micromol m^(-2) s^(-1)}.
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
  This function uses a model for photosynthetic discrimination against 13C in
  C3 plants to determine mesophyll conductance values as described in Busch et
  al. (2020). That paper provides two alternate ways to calculate \code{e_star},
  and two alternate ways to calculate mesophyll conductance \code{gmc}; this
  function allows the user to choose between them. In more detail:

  \itemize{
    \item Isotopic fractionation due to day respiration
          (\code{e_prime = e + e_star}) is calculated with \code{e_star} given
          by either Equation 19 or 20 depending on the value of
          \code{e_star_equation}.
    \item Isotopic discrimination assuming infinite mesophyll conductance
          (\code{Delta_i}) is calculated by setting \code{Cc = Ci} in either
          Equation 2 or 13, depending on the value of \code{gm_type}.
    \item Mesophyll conductance to CO2 (\code{gmc}) is calculated using either
          Equation 21 or 22, depending on the value of \code{gm_type}.
  }

  Note 1: Setting \code{e_star_equation = 19} and \code{gm_type = 'con'} should
  produce identical or similar results to \code{\link{calculate_gm_ubierna}}.

  Note 2: Using \code{e_star_equation = 20} and \code{gm_type = 'dis'} is
  expected to be more accurate, as discussed in Busch et al. (2020); however, be
  aware that this method requires a value for \code{Delta_obs_growth}, which may
  not always be available unless it is intentionally measured.

  References:

  Busch, F. A., Holloway-Phillips, M., Stuart-Williams, H. and Farquhar, G. D.
  "Revisiting carbon isotope discrimination in C3 plants shows respiration rules
  when photosynthesis is low." Nat. Plants 6, 245–258 (2020)
  [\doi{10.1038/s41477-020-0606-6}].
}

\value{
  An \code{exdf} object based on \code{exdf_obj} that includes the following
  additional columns, calculated as described above: \code{e_prime},
  \code{e_star}, \code{Delta_i}, and \code{gmc}, as well as the values of a few
  intermediate calculations such as \code{Delta_i_term_1} and
  \code{Delta_i_term_2}. The category for each of these new columns is
  \code{calculate_gm_busch} to indicate that they were created using this
  function.
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
    'specificity_at_tleaf',
    'M / M',
    value = 90
)

# Calculate Gamma_star (needed for calculate_gm_busch)
licor_data <- calculate_gamma_star(licor_data)

# Calculate isotope discrimination (needed for calculate_gm_busch)
licor_data <- calculate_isotope_discrimination(licor_data)

# Set Delta_obs_growth to the average of Delta_obs_tdl over the first 6 points,
# where the ambient CO2 concentration was set to the atmospheric value (420 ppm)
# (needed for calculate_gm_busch).
licor_data <- set_variable(
  licor_data,
  'Delta_obs_growth',
  'ppt',
  value = mean(licor_data[1:6, 'Delta_obs_tdl'])
)

# Set respiration (needed for calculate_gm_busch)
licor_data <- set_variable(
  licor_data,
  'Rd',
  'micromol m^(-2) s^(-1)',
  value = 1.2
)

# Calculate mesophyll conductance
licor_data <- calculate_gm_busch(licor_data)

# Calculate Cc using the new values of mesophyll conductance
licor_data <- apply_gm(licor_data)

# View some of the results
licor_data[, c('replicate', 'CO2_s', 'Delta_obs_tdl', 'e_prime', 'gmc', 'Ci', 'Cc')]
}

\concept{exdf}
