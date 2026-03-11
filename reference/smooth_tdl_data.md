# Smoothing data from one TDL valve

Tool for applying a smoothing function to the time series corresponding
to measurements from a single valve in a tunable diode laser (TDL) data
set.

## Usage

``` r
smooth_tdl_data(
    tdl_exdf,
    column_to_be_smoothed,
    valve_column_name,
    valve_number,
    smoothing_function
  )
```

## Arguments

- tdl_exdf:

  An `exdf` object representing data from a TDL data logger.

- column_to_be_smoothed:

  The name of the column in `tdl_exdf` that contains the data to be
  smoothed; typically, this is `'Conc12C_Avg'` or `'Conc12C_Avg'`.

- valve_column_name:

  The name of the column in `tdl_exdf` that contains the valve number;
  typically, this is `'valve_number'`.

- valve_number:

  The value of the `valve_column_name` column that indicates the valve
  to be smoothed.

- smoothing_function:

  A function that accepts two vectors `Y` and `X` (in that order) and
  returns a smoothed version of `Y(X)`; typically, `smoothing_function`
  is based on
  [`smooth.spline`](https://rdrr.io/r/stats/smooth.spline.html) or a
  filter from the `signal` package.

## Details

The output from a TDL is highly sensitive to electronic and atmospheric
noise, and it is often helpful to smooth the data from one or more
valves before attempting to apply calibration corrections or determine
the content of an unknown gas mixture. `smooth_tdl_data` is a
convenience function that extracts a time series corresponding to data
from one valve, applies a smoothing operation, and replaces the original
data in `tdl_exdf` with the smoothed version. The smoothing function is
user-supplied to allow more flexbility.

In addition to the `column_to_be_smoothed` and `valve_column_name`
columns, the `tdl_exdf` must also contain an `'elapsed_time'` column,
which is typically created by a call to
[`identify_tdl_cycles`](https://eloch216.github.io/PhotoGEA/reference/identify_tdl_cycles.md).

## Value

An `exdf` object based on `tdl_exdf`, where the time series of
`column_to_be_smoothed` vs. `'elapsed_time'` has been replaced by a
smoothed version obtained by applying the `smoothing_function`.

## Examples

``` r
# Example: Smoothing the 12C signal from one TDL valve using a spline fit
tdl_file <- read_gasex_file(
  PhotoGEA_example_file_path('tdl_sampling_1.dat'),
  'TIMESTAMP'
)

tdl_file <- identify_tdl_cycles(
  tdl_file,
  valve_column_name = 'valve_number',
  cycle_start_valve = 20,
  expected_cycle_length_minutes = 2.7,
  expected_cycle_num_valves = 9,
  timestamp_colname = 'TIMESTAMP'
)

spline_smoothing_function <- function(Y, X) {
    ss <- smooth.spline(X, Y)
    return(ss$y)
}

spline_smoothed_tdl_file <- smooth_tdl_data(
  tdl_file, 'Conc12C_Avg', 'valve_number', 20, spline_smoothing_function
)
```
