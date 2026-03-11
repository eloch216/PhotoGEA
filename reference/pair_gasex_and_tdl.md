# Pair gas exchange and TDL data

Identifies the closest TDL cycle corresponding to each entry in the gas
exchange data and adds the TDL data to the gas exchange data.

## Usage

``` r
pair_gasex_and_tdl(
    gasex_exdf,
    tdl_exdf,
    max_allowed_time_difference = 1,
    gasex_timestamp_column_name = 'time',
    tdl_timestamp_column_name = 'TIMESTAMP'
  )
```

## Arguments

- gasex_exdf:

  An `exdf` object representing data from a photosynthetic gas exchange
  measurement system.

- tdl_exdf:

  An `exdf` object representing calibrated data from a tunable diode
  laser absorption spectroscopy system. Typically this is the output
  from applying
  [`process_tdl_cycle_erml`](https://eloch216.github.io/PhotoGEA/reference/process_tdl_cycle_erml.md)
  or
  [`process_tdl_cycle_polynomial`](https://eloch216.github.io/PhotoGEA/reference/process_tdl_cycle_polynomial.md)
  to a set of uncalibrated TDL data.

- max_allowed_time_difference:

  The maximum time difference (in minutes) to allow between gas exchange
  and TDL timestamp values.

- gasex_timestamp_column_name:

  The name of the column in `gasex_exdf` that contains the timestamp
  values.

- tdl_timestamp_column_name:

  The name of the column in `tdl_exdf` that contains the timestamp
  values.

## Details

When making combined gas exchange and isotope discrimination
measurements using a portable photosynthetic gas exchange system (such
as a Licor LI-6800) coupled with a tunable diode laser (TDL) absorption
spectroscopy system, the TDL's gas handling system cycles through
several gas lines (or sites) by opening and closing valves. When
analyzing such data, a key step is to combine TDL and gas exchange data
that were measured at the same times.

The `pair_gasex_and_tdl` function performs this operation by locating
the TDL cycle whose timestamp is closest to each Licor file entry. Then,
the 12C, 13C, total CO2, and delta_13C values measured by the TDL from
the Licor's sample and reference lines during that cycle are added to
the gas exchange data as new columns.

## Value

An `exdf` object based on `gasex_exdf` that includes TDL values measured
at the same times as the original gas exchange logs. Several new columns
are added: `'cycle_num'`, `'tdl_time_s'`, `'calibrated_12c_s'`,
`'calibrated_13c_s'`, `'total_CO2_s'`, `'delta_C13_s'`, `'tdl_time_r'`,
`'calibrated_12c_r'`, `'calibrated_13c_r'`, `'total_CO2_r'`, and
`'delta_C13_r'`. Variables with `'_s'` in the name refer to TDL
measurements from the Licor sample line, and `'_r'` indicates the
reference line. The category of each new column is `pair_gasex_and_tdl`
to indicate that it was created using this function.

## Examples

``` r
## In this example we load gas exchange and TDL data files, calibrate the TDL
## data, and pair the data tables together

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

# View some of the results
licor_data[, c('A', 'delta_C13_r', 'delta_C13_s', 'total_CO2_r', 'total_CO2_s')]
#>           A delta_C13_r delta_C13_s total_CO2_r total_CO2_s
#> 1  32.01867   -5.438269   -2.306566    702.8716    426.2619
#> 2  31.91890   -5.344929   -2.178155    708.0634    429.6733
#> 3  31.85562   -5.169295   -2.109784    702.7235    427.6489
#> 4  31.76382   -5.212031   -2.105925    702.7264    428.1308
#> 5  31.69923   -5.117179   -1.660248    702.6373    428.8285
#> 6  31.57078   -5.119033   -1.789007    702.7397    430.1372
#> 7  20.11287   -4.778251   -2.248652    445.1621    268.6527
#> 8  20.13905   -4.946493   -2.257818    448.4566    268.1469
#> 9  20.15698   -4.493983   -1.997756    444.8754    267.8980
#> 10 20.16227   -4.754437   -2.246718    444.6889    267.8129
#> 11 20.16253   -4.914201   -1.980538    444.6759    267.9302
#> 12 20.12078   -5.024895   -2.335218    444.8834    268.4257
```
