# Calculate the total pressure in bar

Calculates the total pressure in `bar`. Licor gas exchange measurement
systems report both the abient air pressure (`Pa`) and the chamber
overpressure (`DeltaPcham`) in `kPa`; the total pressure in the chamber
is therefore given by the sum of these two columns. This function can
accomodate alternative column names for the variables taken from Licor
log files in case they change at some point in the future. This function
also checks the units of each required column and will produce an error
if any units are incorrect.

## Usage

``` r
calculate_total_pressure(
    exdf_obj,
    pa_column_name = 'Pa',
    deltapcham_column_name = 'DeltaPcham'
  )
```

## Arguments

- exdf_obj:

  An `exdf` object that contains pressure measurements.

- pa_column_name:

  The name of the column in `exdf_obj` that contains the ambient air
  pressure in `kPa`.

- deltapcham_column_name:

  The name of the column in `exdf_obj` that contains the chamber
  overpressure in `kPa`.

## Details

If `deltapcham_column_name` is `NA`, this function will simply convert
the values in the `pa_column_name` to units of `bar`. Otherwise, the
values from the `pa_column_name` and `deltapcham_column_name` columns
will be added together and converted to `bar`.

## Value

An `exdf` object based on `exdf_obj` that includes the total pressure
values in a new column called `total_pressure`. The category of this new
column is `calculate_total_pressure` to indicate that it was created
using this function.

## Examples

``` r
# Read an example Licor file included in the PhotoGEA package and calculate the
# total pressure.
licor_file <- read_gasex_file(
  PhotoGEA_example_file_path('ball_berry_1.xlsx')
)

licor_file <- calculate_total_pressure(licor_file)

licor_file$units$total_pressure      # View the units of the new `total_pressure` column
#> [1] "bar"
licor_file$categories$total_pressure # View the category of the new `total_pressure` column
#> [1] "calculate_total_pressure"
licor_file[, 'total_pressure']       # View the values of the new `total_pressure` column
#>  [1] 0.9925101 0.9925749 0.9926076 0.9926272 0.9924095 0.9924228 0.9924457
#>  [8] 0.9925986 0.9926574 0.9928725 0.9929653 0.9929991 0.9928978 0.9929064
#> [15] 0.9929203 0.9929226 0.9929827 0.9929333 0.9929934 0.9926874 0.9926688
#> [22] 0.9923497 0.9924241 0.9923444 0.9923882 0.9923838 0.9923927 0.9921658
```
