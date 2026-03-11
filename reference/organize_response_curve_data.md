# Reorganize response curve data for analysis and plotting

Prepares a set of response curves for future processing and analysis by
numbering and reordering the points, (optionally) removing recovery
points, and (optionally) calculating average values of key variables
across each curve.

## Usage

``` r
organize_response_curve_data(
    licor_exdf,
    identifier_columns,
    measurement_numbers_to_remove,
    column_for_ordering,
    ordering_column_tolerance = Inf,
    columns_to_average = c(),
    print_information = TRUE
  )
```

## Arguments

- licor_exdf:

  An `exdf` object representing response curve data from a Licor gas
  exchange measurement system.

- identifier_columns:

  A vector or list of strings representing the names of columns in
  `licor_exdf` that, taken together, uniquely identify each curve. This
  often includes names like `plot`, `event`, `replicate`, etc.

- measurement_numbers_to_remove:

  A vector of integers specifying which points to remove from each
  curve; for example, if each curve has 16 points and the 10^th^ and
  11^th^ points along the sequence should not be included in subsequent
  analysis, `measurement_numbers_to_remove` could be specified as
  `c(10, 11)`. If `measurement_numbers_to_remove` is set to
  [`c()`](https://rdrr.io/r/base/c.html), no points will be removed.

- column_for_ordering:

  The name of a column that is systematically varied to produce each
  curve; for example, in a light response curve, this would typically by
  `Qin`.

- ordering_column_tolerance:

  To be passed to
  [`check_response_curve_data`](https://eloch216.github.io/PhotoGEA/reference/check_response_curve_data.md)
  as the `driving_column_tolerance` input argument.

- columns_to_average:

  A list of columns whose average values should be calculated; see below
  for details.

- print_information:

  To be passed to
  [`check_response_curve_data`](https://eloch216.github.io/PhotoGEA/reference/check_response_curve_data.md).

## Details

For an `exdf` object consisting of multiple response curves that can be
identified using the values of its `identifier_columns`, this function
performs the following actions:

- Assigns a sequential number to each measurement in each curve,
  beginning with 1. In other words, the first point in the curve is
  given number 1, the second is given number 2, etc. These numbers are
  stored as a new column called `seq_num`.

- (Optionally) extracts a subset of the data. If
  `measurement_numbers_to_remove` is
  [`c()`](https://rdrr.io/r/base/c.html), then this step will be
  skipped; otherwise, values of `seq_num` specified by
  `measurement_numbers_to_remove` will be removed, and then
  `check_response_curve_data` will be called to make sure the remaining
  points all follow the same sequence of setpoint values (within the
  tolerance set by `ordering_column_tolerance`), treating the
  `column_for_ordering` as the `driving_column`.

- Reorders the data according to ascending values of the
  `column_for_ordering`.

- (Optionally) calculates average values of important columns. If
  `columns_to_average` is [`c()`](https://rdrr.io/r/base/c.html), then
  this step will be skipped; otherwise, for each curve, the mean value
  of each column specified in `columns_to_average` will be stored in a
  new column whose name is based on the original column name, but with
  `'_avg'` added at the end. For example, the average value of the `Qin`
  column would be stored in `Qin_avg`.

Removing certain points is often helpful for A-Ci curves, where the
CO~2~ concentration begins at the ambient value, is decreased to a low
value, is reset to atmospheric for several measurements to allow the
plant to reacclimate, and then is increased to higher values. In this
case, only the first measurement at ambient CO~2~ is used for plotting
or additional analysis, and the "recovery" points should be removed.

Reordering the points is often helpful for plotting. For example, the
points in an A-Ci curve would not be ordered according to their Ci
values in a curve measured using a sequence as described above. This can
cause issues when making line plots, so it may be convenient to reorder
them according to their Ci values.

Calculating average values of certain columns is especially useful for
estimating `Jmax` values using
[`calculate_jmax`](https://eloch216.github.io/PhotoGEA/reference/calculate_jmax.md),
since this operation requires average values of leaf temperature and
incident photon flux across each curve.

## Value

An `exdf` object based on `licor_exdf` but processed as described above.

## Examples

``` r
# Read an example Licor file included in the PhotoGEA package and organize it.
# This file includes several 7-point light-response curves that can be uniquely
# identified by the values of its 'species' and 'plot' columns. Since these are
# light-response curves, each one follows a pre-set sequence of `Qin` values.
licor_file <- read_gasex_file(
  PhotoGEA_example_file_path('ball_berry_1.xlsx')
)

# Split the data into individual curves, keep all seven measurement points in
# each curve, and order them by their incident light values (since these are
# light response curves). The curves were measured from high to low values of
# `Qin`, so after organizing the curves, their order will be reversed from the
# original version. Also add the average value of TleafCnd and Qin for each
# curve.
licor_file <- organize_response_curve_data(
  licor_file,
  c('species', 'plot'),
  c(),
  'Qin',
  columns_to_average = c('TleafCnd', 'Qin')
)

# View a subset of the data, including the new `seq_num` column
print(licor_file[, c('species', 'plot', 'seq_num', 'Qin', 'A', 'Qin_avg'), TRUE])
#> 
#> Converting an `exdf` object to a `data.frame` before printing
#> 
#>    species [UserDefCon] (NA) plot [UserDefCon] (NA) seq_num [NA] (NA)
#> 1                    soybean                     1a                 7
#> 2                    soybean                     1a                 6
#> 3                    soybean                     1a                 5
#> 4                    soybean                     1a                 4
#> 5                    soybean                     1a                 3
#> 6                    soybean                     1a                 2
#> 7                    soybean                     1a                 1
#> 8                    soybean                     1b                 7
#> 9                    soybean                     1b                 6
#> 10                   soybean                     1b                 5
#> 11                   soybean                     1b                 4
#> 12                   soybean                     1b                 3
#> 13                   soybean                     1b                 2
#> 14                   soybean                     1b                 1
#> 15                   tobacco                      2                 7
#> 16                   tobacco                      2                 6
#> 17                   tobacco                      2                 5
#> 18                   tobacco                      2                 4
#> 19                   tobacco                      2                 3
#> 20                   tobacco                      2                 2
#> 21                   tobacco                      2                 1
#> 22                   soybean                      5                 7
#> 23                   soybean                      5                 6
#> 24                   soybean                      5                 5
#> 25                   soybean                      5                 4
#> 26                   soybean                      5                 3
#> 27                   soybean                      5                 2
#> 28                   soybean                      5                 1
#>    Qin [LeafQ] (micromol m^(-2) s^(-1)) A [GasEx] (micromol m^(-2) s^(-1))
#> 1                               200.141                           5.679174
#> 2                               350.106                          12.240074
#> 3                               500.135                          16.617929
#> 4                               649.951                          23.290614
#> 5                               799.983                          25.430106
#> 6                              1099.800                          22.996635
#> 7                              2000.030                          35.400473
#> 8                               200.070                           9.276031
#> 9                               350.060                          14.666531
#> 10                              499.982                          17.009070
#> 11                              650.070                          26.383035
#> 12                              799.873                          25.472460
#> 13                             1099.890                          30.250627
#> 14                             1999.850                          47.950460
#> 15                              200.236                           7.849870
#> 16                              350.073                          12.044692
#> 17                              500.036                          12.106347
#> 18                              650.001                          16.345888
#> 19                              800.046                          22.218165
#> 20                             1099.980                          25.492278
#> 21                             2000.020                          29.614535
#> 22                              200.029                          10.768250
#> 23                              349.991                          17.091593
#> 24                              499.921                          22.068194
#> 25                              650.108                          22.924484
#> 26                              799.910                          32.769396
#> 27                             1099.960                          36.218245
#> 28                             1999.980                          45.780746
#>    Qin_avg [organize_response_curve_data] (micromol m^(-2) s^(-1))
#> 1                                                         800.0209
#> 2                                                         800.0209
#> 3                                                         800.0209
#> 4                                                         800.0209
#> 5                                                         800.0209
#> 6                                                         800.0209
#> 7                                                         800.0209
#> 8                                                         799.9707
#> 9                                                         799.9707
#> 10                                                        799.9707
#> 11                                                        799.9707
#> 12                                                        799.9707
#> 13                                                        799.9707
#> 14                                                        799.9707
#> 15                                                        800.0560
#> 16                                                        800.0560
#> 17                                                        800.0560
#> 18                                                        800.0560
#> 19                                                        800.0560
#> 20                                                        800.0560
#> 21                                                        800.0560
#> 22                                                        799.9856
#> 23                                                        799.9856
#> 24                                                        799.9856
#> 25                                                        799.9856
#> 26                                                        799.9856
#> 27                                                        799.9856
#> 28                                                        799.9856
```
