# Calculate leakiness

Calculates leakiness (`phi`) from combined gas exchange and isotope
discrimination measurements as described in Ubierna et al. (2013). This
function can accomodate alternative colum names for the variables taken
from `exdf_obj`; it also checks the units of each required column and
will produce an error if any units are incorrect.

## Usage

``` r
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
```

## Arguments

- exdf_obj:

  An `exdf` object.

- e:

  The isotopic fractionation during day respiration in `ppt`.

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

- rl_column_name:

  The name of the column in `exdf_obj` that contains the rate of day
  respiration in `micromol m^(-2) s^(-1)`.

- t_column_name:

  The name of the column in `exdf_obj` that contains the ternary
  correction factor (`dimensionless`). Values of `t` are typically
  calculated using
  [`calculate_ternary_correction`](https://eloch216.github.io/PhotoGEA/reference/calculate_ternary_correction.md)

## Details

This function uses the model for photosynthetic discrimination against
13C in C4 plants to determine leakiness values, as described in Ubierna
et al. (2013). In particular, the following equations from that source
are implemented in the code:

- Isotopic fractionation due to day respiration (`e_prime`) is
  calculated using Equation 21.

- Leakiness including respiratory and photorespiratory fractionations
  under high light (`phi_i`) is calculated using Equation 16.

- Leakiness including respiratory and photorespiratory fractionations
  and Cs under high light (`phi_is`) is calculated using Equation 15.

- Leakiness ignoring respiratory and photorespiratory fractionations and
  Cs (`phi_sim`) is calculated using Equation 17.

References:

Ubierna, N., Sun, W., Kramer, D. M. and Cousins, A. B. "The efficiency
of C4 photosynthesis under low light conditions in Zea mays, Miscanthus
x giganteus and Flaveria bidentis." Plant, Cell & Environment 36,
365–381 (2013)
\[[doi:10.1111/j.1365-3040.2012.02579.x](https://doi.org/10.1111/j.1365-3040.2012.02579.x)
\].

## Value

An `exdf` object based on `exdf_obj` that includes the following
additional columns, calculated as described above: `e_prime`, `phi_i`,
`phi_is`, and `phi_sim`. The category for each of these new columns is
`calculate_leakiness_ubierna` to indicate that they were created using
this function.

## Examples

``` r
## In this example we load gas exchange and TDL data files, calibrate the TDL
## data, pair the data tables together, and then calculate leakiness. The
## results from this example are not meaningful because these measurements
## were not collected from C4 plants.

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
#>    replicate   CO2_s Delta_obs_tdl     phi_i   phi_sim
#> 1          1 417.363      8.039825 0.5647984 0.5652561
#> 2          1 420.552      8.137268 0.5787469 0.5791113
#> 3          1 418.796      7.893969 0.5548495 0.5550261
#> 4          1 419.493      8.029239 0.5625286 0.5627514
#> 5          1 420.102      8.964915 0.6097656 0.6098926
#> 6          1 421.133      8.673921 0.5935674 0.5936957
#> 7          1 262.873      6.434975 0.4796541 0.4793111
#> 8          1 262.720      6.747186 0.4966989 0.4966153
#> 9          1 262.633      6.326945 0.4760775 0.4752982
#> 10         1 262.271      6.358834 0.4709619 0.4705850
#> 11         1 262.112      7.450206 0.5193723 0.5192365
#> 12         1 262.176      6.843263 0.4874634 0.4875021
```
