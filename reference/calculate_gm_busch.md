# Calculate mesophyll conductance to CO2 diffusion

Calculates mesophyll conductance to CO2 diffusion (`gmc`) from combined
gas exchange and isotope discrimination measurements as described in
Busch et al. (2020). This function can accomodate alternative colum
names for the variables taken from `exdf_obj`; it also checks the units
of each required column and will produce an error if any units are
incorrect.

## Usage

``` r
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
    gamma_star_column_name = 'Gamma_star_tl',
    rl_column_name = 'RL',
    total_pressure_column_name = 'total_pressure',
    t_column_name = 't'
  )
```

## Arguments

- exdf_obj:

  An `exdf` object.

- e:

  The isotopic fractionation during day respiration in `ppt`.

- f:

  The isotopic fractionation during photorespiration in `ppt`.

- e_star_equation:

  The equation from Busch et al. (2020) to use for calculating `e_star`;
  must be `19` or `20`.

- gm_type:

  Determines whether day respiration is assumed to be isotopically
  connected to the CBB cycle (`gm_type = 'con'`) or isotopically
  disconnected from the CBB cycle (`gm_type = 'dis'`). This choice will
  determine which equations are used to calculate mesophyll conductance;
  when `gm_type` is `'con'`, Equations 2 and 21 will be used; otherwise,
  Equations 13 and 22 will be used.

- a_bar_column_name:

  The name of the column in `exdf_obj` that contains the weighted
  isotopic fractionation across the boundary layer and stomata in `ppt`.
  Values of `a_bar` are typically calculated using
  [`calculate_ternary_correction`](https://eloch216.github.io/PhotoGEA/reference/calculate_ternary_correction.md).

- a_column_name:

  The name of the column in `exdf_obj` that contains the net CO2
  assimilation rate in `micromol m^(-2) s^(-1)`.

- ci_column_name:

  The name of the column in `exdf_obj` that contains the intercellular
  CO2 concentration in `micromol mol^(-1)`.

- co2_s_column_name:

  The name of the column in `exdf_obj` that contains the CO2
  concentration in the sample line (outgoing air) in
  `micromol mol^(-1)`.

- csurface_column_name:

  The name of the column in `exdf_obj` that contains the CO2
  concentration at the leaf surface in `micromol mol^(-1)`. Values of
  `Csurface` are typically calculated using
  [`calculate_gas_properties`](https://eloch216.github.io/PhotoGEA/reference/calculate_gas_properties.md).

- delta_c13_r_column_name:

  The name of the column in `exdf_obj` that contains the CO2 isotope
  ratio in the reference line (incoming air) in `ppt`.

- delta_obs_growth_column_name:

  The name of the column in `exdf_obj` that contains the observed
  discrimination under the typical CO2 concentration in the plant's
  environment during its growth (in `ppt`). This is only required when
  using Equation 20 for `e_star` (see `e_star_equation`).

- delta_obs_tdl_column_name:

  The name of the column in `exdf_obj` that contains the observed
  isotope discrimination values in `ppt`.

- gamma_star_column_name:

  The name of the column in `exdf_obj` that contains the chloroplastic
  CO2 concentration at which CO2 gains from Rubisco carboxylation are
  exactly balanced by CO2 losses from Rubisco oxygenation, at leaf
  temperature, expressed in `micromol mol^(-1)`. Values of `Gamma_star`
  at leaf temperature are typically calculated using
  [`calculate_gamma_star`](https://eloch216.github.io/PhotoGEA/reference/calculate_gamma_star.md)
  or
  [`calculate_temperature_response`](https://eloch216.github.io/PhotoGEA/reference/calculate_temperature_response.md).

- rl_column_name:

  The name of the column in `exdf_obj` that contains the rate of
  non-photorespiratory CO2 release in the light, in
  `micromol m^(-2) s^(-1)`.

- total_pressure_column_name:

  The name of the column in `exdf_obj` that contains the total pressure
  in `bar`.

- t_column_name:

  The name of the column in `exdf_obj` that contains the ternary
  correction factor (`dimensionless`). Values of `t` are typically
  calculated using
  [`calculate_ternary_correction`](https://eloch216.github.io/PhotoGEA/reference/calculate_ternary_correction.md)

## Details

This function uses a model for photosynthetic discrimination against 13C
in C3 plants to determine mesophyll conductance values as described in
Busch et al. (2020). That paper provides two alternate ways to calculate
`e_star`, and two alternate ways to calculate mesophyll conductance
`gmc`; this function allows the user to choose between them. In more
detail:

- Isotopic fractionation due to day respiration (`e_prime = e + e_star`)
  is calculated with `e_star` given by either Equation 19 or 20
  depending on the value of `e_star_equation`.

- Isotopic discrimination assuming infinite mesophyll conductance
  (`Delta_i`) is calculated by setting `Cc = Ci` in either Equation 2 or
  13, depending on the value of `gm_type`.

- Mesophyll conductance to CO2 (`gmc`) is calculated using either
  Equation 21 or 22, depending on the value of `gm_type`.

Note 1: Setting `e_star_equation = 19` and `gm_type = 'con'` should
produce identical or similar results to
[`calculate_gm_ubierna`](https://eloch216.github.io/PhotoGEA/reference/calculate_gm_ubierna.md).

Note 2: Using `e_star_equation = 20` and `gm_type = 'dis'` is expected
to be more accurate, as discussed in Busch et al. (2020); however, be
aware that this method requires a value for `Delta_obs_growth`, which
may not always be available unless it is intentionally measured.

References:

Busch, F. A., Holloway-Phillips, M., Stuart-Williams, H. and Farquhar,
G. D. "Revisiting carbon isotope discrimination in C3 plants shows
respiration rules when photosynthesis is low." Nat. Plants 6, 245–258
(2020)
\[[doi:10.1038/s41477-020-0606-6](https://doi.org/10.1038/s41477-020-0606-6)
\].

## Value

An `exdf` object based on `exdf_obj` that includes the following
additional columns, calculated as described above: `e_prime`, `e_star`,
`Delta_i`, and `gmc`, as well as the values of a few intermediate
calculations such as `Delta_i_term_1` and `Delta_i_term_2`. The category
for each of these new columns is `calculate_gm_busch` to indicate that
they were created using this function.

## Examples

``` r
## In this example we load gas exchange and TDL data files, calibrate the TDL
## data, pair the data tables together, and then calculate mesophyll conductance

# Read the TDL data file, making sure to interpret the time zone as US Central
# time
tdl_data <- read_gasex_file(
  PhotoGEA_example_file_path('tdl_for_gm.dat'),
  'TIMESTAMP',
  list(tz = 'America/Chicago')
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
  list(tz = 'America/Chicago')
)

# Get TDL valve information from Licor file name; for this TDL system, the
# reference valve is 12 when the sample valve is 11
licor_data <- get_sample_valve_from_filename(licor_data, list('11' = 12))

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
  'RL',
  'micromol m^(-2) s^(-1)',
  value = 1.2
)

# Calculate mesophyll conductance
licor_data <- calculate_gm_busch(licor_data)

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
licor_data[, c('replicate', 'CO2_s', 'Delta_obs_tdl', 'e_prime', 'gmc', 'Ci', 'Cc')]
#>    replicate   CO2_s Delta_obs_tdl     e_prime       gmc       Ci        Cc
#> 1          1 417.363      8.039825 -0.18823755 0.1741137 286.7663 101.17189
#> 2          1 420.552      8.137268 -0.19234146 0.1842765 276.9970 102.18116
#> 3          1 418.796      7.893969  0.22659156 0.1673504 291.8317  99.70574
#> 4          1 419.493      8.029239  0.04858619 0.1692784 291.1437 101.74081
#> 5          1 420.102      8.964915 -0.79223714 0.1802174 293.7559 116.20312
#> 6          1 421.133      8.673921 -0.50309810 0.1734677 295.9070 112.18043
#> 7          1 262.873      6.434975  2.07663004 0.1596937 182.7755  55.63881
#> 8          1 262.720      6.747186  1.59617713 0.1656192 181.1386  58.38280
#> 9          1 262.633      6.326945  2.46892812 0.1623944 179.6180  54.31518
#> 10         1 262.271      6.358834  2.17658584 0.1497082 191.7498  55.81379
#> 11         1 262.112      7.450206  0.92544892 0.1536710 199.2717  66.86048
#> 12         1 262.176      6.843263  1.42169819 0.1425367 204.1826  61.70023
```
