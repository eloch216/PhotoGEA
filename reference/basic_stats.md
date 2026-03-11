# Calculate basic stats (mean and standard error)

Calculates basic stats (mean and standard error) for each applicable
column in an `exdf` object split up according to the values of one or
more identifier columns.

## Usage

``` r
basic_stats(
    exdf_obj,
    identifier_columns,
    na.rm = TRUE
  )
```

## Arguments

- exdf_obj:

  An `exdf` object.

- identifier_columns:

  The name(s) of one or more columns in a vector or list that can be
  used to split `exdf_obj` into chunks.

- na.rm:

  A logical value indicating whether or not to remove NA values before
  calculating means and standard errors.

## Details

This function first splits up `exdf_obj` into chunks according to the
values of the `identifier_columns`. For each chunk, columns that have a
single unique value are identified and excluded from the statistical
calculations. For the remaining numeric columns, the mean and standard
error are calculated.

## Value

An `exdf` object including the mean and standard error for each
applicable column, where each row represents one value of the
`identifier_columns`. The column names are determined by appending
`'_avg'` and `'_stderr'` to the original names.

## Examples

``` r
# Read an example Licor file included in the PhotoGEA package
licor_file <- read_gasex_file(
  PhotoGEA_example_file_path('ball_berry_1.xlsx')
)

# Calculate the average assimilation and stomatal conductance values for each
# species. (Note: this is not a meaningful calculation!)
basic_stats(
  licor_file[ , c('species', 'K', 'A', 'gsw'), TRUE],
  'species'
)
#> 
#> Converting an `exdf` object to a `data.frame` before printing
#> 
#>   species [UserDefCon] (NA) K [Const] (NA)
#> 1                   soybean            0.5
#> 2                   tobacco            0.5
#>   A_avg [GasEx] (micromol m^(-2) s^(-1))
#> 1                               23.82305
#> 2                               17.95311
#>   A_stderr [GasEx] (micromol m^(-2) s^(-1)) gsw_avg [GasEx] (mol m^(-2) s^(-1))
#> 1                                  2.470260                           0.4043737
#> 2                                  3.026751                           0.2192409
#>   gsw_stderr [GasEx] (mol m^(-2) s^(-1))
#> 1                             0.04320554
#> 2                             0.02502419
```
