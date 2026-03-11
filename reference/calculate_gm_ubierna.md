# Calculate mesophyll conductance to CO2 diffusion

Calculates mesophyll conductance to CO2 diffusion (`gmc`) from combined
gas exchange and isotope discrimination measurements as described in
Ubierna et al. (2018). This function can accomodate alternative colum
names for the variables taken from `exdf_obj`; it also checks the units
of each required column and will produce an error if any units are
incorrect.

## Usage

``` r
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
```

## Arguments

- exdf_obj:

  An `exdf` object.

- e:

  The isotopic fractionation during day respiration in `ppt`.

- f:

  The isotopic fractionation during photorespiration in `ppt`.

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

This function uses the comprehensive model for photosynthetic
discrimination against 13C in C3 plants to calculate mesophyll
conductance, as described in Ubierna et al. (2018). In particular, the
following equations from that source are implemented in the code:

- Isotopic fractionation due to day respiration (`e_prime`) is
  calculated using Equations 28 and 30.

- Isotopic discrimination due to photorespiration (`Delta_f`), due to
  day respiration (`Delta_e`), and that would occur if Ci = Cc in the
  absence of any respiratory fractionation (`Delta_i`) are calculated
  using Equations 34, 33, and 31, respectively.

- Mesophyll conductance to CO2 diffusion (`gmc`) is calculated using
  Equation 44. This equation is broken up into two factors called
  `Delta_difference` and `equation_top` which are separately returned in
  the output from `calculate_gm_ubierna`.

For an alternative method for calculating `gmc`, see
[`calculate_gm_busch`](https://eloch216.github.io/PhotoGEA/reference/calculate_gm_busch.md).

References:

Ubierna, N., Holloway-Phillips, M.-M. and Farquhar, G. D. "Using Stable
Carbon Isotopes to Study C3 and C4 Photosynthesis: Models and
Calculations." in Photosynthesis: Methods and Protocols (ed. Covshoff,
S.) 155–196 (Springer, 2018)
\[[doi:10.1007/978-1-4939-7786-4_10](https://doi.org/10.1007/978-1-4939-7786-4_10)
\].

## Value

An `exdf` object based on `exdf_obj` that includes the following
additional columns, calculated as described above: `e_prime`, `Delta_i`,
`Delta_e`, `Delta_f`, `Delta_difference`, `equation_top`, and `gmc`. The
category for each of these new columns is `calculate_gm_ubierna` to
indicate that they were created using this function.

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
#>    replicate   CO2_s Delta_obs_tdl       gmc       Ci        Cc
#> 1          1 417.363      8.039825 0.1741069 286.7663 101.16460
#> 2          1 420.552      8.137268 0.1842769 276.9970 102.18154
#> 3          1 418.796      7.893969 0.1673113 291.8317  99.66079
#> 4          1 419.493      8.029239 0.1692570 291.1437 101.71697
#> 5          1 420.102      8.964915 0.1803390 293.7559 116.32286
#> 6          1 421.133      8.673921 0.1735334 295.9070 112.25005
#> 7          1 262.873      6.434975 0.1593770 182.7755  55.38615
#> 8          1 262.720      6.747186 0.1653412 181.1386  58.17641
#> 9          1 262.633      6.326945 0.1620221 179.6180  54.02724
#> 10         1 262.271      6.358834 0.1494173 191.7498  55.54918
#> 11         1 262.112      7.450206 0.1535156 199.2717  66.72640
#> 12         1 262.176      6.843263 0.1423399 204.1826  61.50326
```
