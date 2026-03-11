# Check response curve data for common issues

Checks to make sure an
[`exdf`](https://eloch216.github.io/PhotoGEA/reference/exdf.md) object
representing multiple response curves meets basic expectations.

## Usage

``` r
check_response_curve_data(
    exdf_obj,
    identifier_columns,
    expected_npts = 0,
    driving_column = NULL,
    driving_column_tolerance = 1.0,
    col_to_ignore_for_inf = 'gmc',
    constant_col = list(),
    error_on_failure = TRUE,
    print_information = TRUE
  )
```

## Arguments

- exdf_obj:

  An `exdf` object representing multiple response curves.

- identifier_columns:

  A vector or list of strings representing the names of columns in
  `exdf_obj` that, taken together, uniquely identify each curve. This
  often includes names like `plot`, `event`, `replicate`, etc.

- expected_npts:

  A numeric vector of length 1 or 2 specifying conditions for the number
  of points in each curve. If `expected_npts` is set to a negative
  number, then this check will be skipped. See below for more details.

- driving_column:

  The name of a column that is systematically varied to produce each
  curve; for example, in a light response curve, this would typically by
  `Qin`. If `driving_column` is `NULL`, then this check will be skipped.

- driving_column_tolerance:

  An absolute tolerance for the deviation of each value of
  `driving_column` away from its mean across all the curves; the
  `driving_column_tolerance` can be set to `Inf` to disable this check.

- col_to_ignore_for_inf:

  Any columns to ignore while checking for infinite values. Mesophyll
  conductance (`gmc`) is often set to infinity intentionally so should
  be ignored when performing this check. To completely disable this
  check, set `col_to_ignore_for_inf` to `NULL`.

- constant_col:

  A list of named numeric elements, where the name indicates a column of
  `exdf_obj` that should be constant, and the value indicates whether
  the column's values must be identical or whether they must lie within
  a specified numeric range. If `constant_col` is an empty list, then
  this check will be skipped. See below for more details.

- error_on_failure:

  A logical value indicating whether to send an error message when an
  issue is detected. See details below.

- print_information:

  A logical value indicating whether to print additional information to
  the R terminal when an issue is detected. See details below.

## Details

**Basic Behavior:**

This function makes a few basic checks to ensure that the response curve
data includes the expected information and does not include any
mistakes. If no problems are detected, this function will be silent with
no return value. If a problem is detected, then the user will be
notified in one or more ways:

- If `error_on_failure` is `TRUE`, then this function will throw an
  error with a short message. If `print_information` is also `TRUE`,
  then additional information will be printed to the R terminal.

- If `error_on_failure` is `FALSE` and `print_information` is also
  `FALSE`, then this function will throw a warning with a short message.

- If `error_on_failure` is `FALSE` and `print_information` is true,
  information about the problem will be printed to the R terminal.

This function will (optionally) perform several checks:

- Checking for infinite values: If `col_to_ignore_for_inf` is not
  `NULL`, no numeric columns in `exdf_obj` should have infinite values,
  with the exception of columns designated in `col_to_ignore_for_inf`.

- Checking required columns: All elements of `identifier_columns` should
  be present as columns in `exdf_obj`. If `driving_column` is not
  `NULL`, it should also be present as a column in `exdf_obj`. If
  `constant_col` is not empty, then these columns must also be present
  in `exdf_obj`.

- Checking the number of points in each curve: The general idea is to
  ensure that each curve has the expected number of points. Several
  options can be specified via the value of `expected_npts`, as
  discussed below.

- Checking the driving column: If `driving_column` is not `NULL`, then
  each curve should have the same sequence of values in this column. To
  allow for small variations, a nonzero `driving_column_tolerance` can
  be specified.

- Checking the constant columns: If `constant_col` is not empty, then
  each specified column should either be constant, or only vary by a
  specified amount. See details below.

By default, most of these are not performed (except the simplest ones
like checking for infinite values or checking that key columns are
present). This enables an "opt-in" use style, where users can specify
just the checks they wish to make.

**More Details:**

There are several options for checking the number of points in each
curve:

- If `expected_npts` is a single negative number, no check will be
  performed.

- If `expected_npts` is 0, then each curve is expected to have the same
  number of points.

- If `expected_npts` is a single positive number, then each curve is
  expected to have that many points. For example, if `expected_npts` is
  7, then each curve must have 7 points.

- If `expected_npts` is a pair of positive numbers, then each curve is
  expected to have a number of points lying within the range defined by
  `expected_npts`. For example, if `expected_npts` is `c(6, 8)`, then
  each curve must have no fewer than 6 points and no more than 8 points.

- If `expected_npts` is a pair of numbers, one of which is zero and one
  of which is positive, then the positive number specifies a range; each
  curve must differ from the average number of points by less than the
  range. For example, if `expected_npts` is `c(0, 3)`, then every curve
  must have a number of points within 3 of the average number of points.

There are two options for checking columns that should be constant:

- A value of `NA` indicates that all values of that column must be
  exactly identical; this check applies for numeric and character
  columns.

- A numeric value indicates that the range of values of that column must
  be smaller than the specified range; this range applies for numeric
  columns only.

For example, setting `constant_col = list(species = NA, Qin = 10)` means
that each curve must have only a single value of the `species` column,
and that the value of the `Qin` column cannot vary by more than 10
across each curve.

**Use Cases:**

Using `check_response_curve_data` is not strictly necessary, but it can
be helpful both to you and to anyone else reading your analysis code.
Here are a few typical use cases:

- **Average response curves:** It is common to calculate and plot
  average response curves, either manually or by using
  [`xyplot_avg_rc`](https://eloch216.github.io/PhotoGEA/reference/xyplot_avg_rc.md).
  But, it only makes sense to do this if each curve followed the same
  sequence of the driving variable. In this case,
  `check_response_curve_data` can be used to confirm that each curve
  used the same values of `CO2_r_sp` (for an A-Ci curve) or `Qin` (for
  an A-Q curve).

- **Removing recovery points:** It is common to wish to remove one or
  more recovery points from a set of curves. The safest way to do this
  is to confirm that all the curves use the same sequence of setpoints;
  then you can be sure that, for example, points 9 and 10 are the
  recovery points in every curve.

- **Making a statement of expectations:** If you measured a set of A-Ci
  curves where each curve has 16 points and used the same sequence of
  `CO2_r` setpoints, you could record this somewhere in your notes. But
  it would be even more meaningful to use `check_response_curve_data` in
  your script with `expected_npts` set to 16. If this check passes, then
  it means not only that your claim is correct, but also that the
  identifier columns are being interpreted properly.

- **Checking identifiers:** If the data set includes some identifying
  metadata, such as a species or location, it may be helpful to confirm
  that each curve has a single value of these "identifier" columns.
  Otherwise, the data set may be difficult to interpret.

- **Checking measurement conditions:** If the response curves are
  expected to be measured under constant temperature, humidity, light,
  or other environmental variables, it may be helpful to confirm that
  these variables do not vary too much across each individual curve.
  Otherwise, parameter values estimated from the curves may not be
  meaningful.

Sometimes the response curves in a large data set were not all measured
with the same sequence of setpoints. If only a few different sequences
were used, it is possible to split them into groups and separately run
`check_response_curve_data` on each group. This scenario is discussed in
the Frequently Asked Questions vignette.

Even if none of the above situations are relevant to you, it may still
be helpful to run `run check_response_curve_data` but with
`expected_npts` set to 0 and `error_on_failure` set to `FALSE`. With
these settings, if there are curves with different numbers of points,
the function will print the number of points in each curve to the R
terminal, but won't stop the rest of the script from running. This can
be useful for detecting problems with the `curve_identifier` column. For
example, if the longest curves in the set are known to have 17 points,
but `check_response_curve_data` identifies a curve with 34 points, it is
clear that the same identifier was accidentally used for two different
curves.

## Value

The `check_response_curve_data` function does not return anything.

## Examples

``` r
# Read an example Licor file included in the PhotoGEA package and check it.
# This file includes several 7-point light-response curves that can be uniquely
# identified by the values of its 'species' and 'plot' columns. Since these are
# light-response curves, each one follows a pre-set sequence of `Qin` values.
licor_file <- read_gasex_file(
  PhotoGEA_example_file_path('ball_berry_1.xlsx')
)

# Make sure there are no infinite values and that all curves have the same
# number of points
check_response_curve_data(licor_file, c('species', 'plot'))

# Make sure there are no inifinite values and that all curves have 7 points
check_response_curve_data(licor_file, c('species', 'plot'), expected_npts = 7)

# Make sure there are no infinite values, that all curves have 7 points, and
# that the values of the `Qin` column follow the same sequence in all curves
# (to within 1.0 micromol / m^2 / s)
check_response_curve_data(
  licor_file,
  c('species', 'plot'),
  expected_npts = 7,
  driving_column = 'Qin',
  driving_column_tolerance = 1.0
)

# Make sure that there are no infinite values and that all curves have between
# 8 and 10 points; this will intentionally fail
check_response_curve_data(
  licor_file,
  c('species', 'plot'),
  expected_npts = c(8, 10),
  error_on_failure = FALSE
)
#>   species plot npts
#> 1 soybean   1a    7
#> 2 soybean   1b    7
#> 3 tobacco    2    7
#> 4 soybean    5    7
#> 
#> The following curves have fewer than 8 points or greater than 10 points: soybean.1a, soybean.1b, tobacco.2, soybean.5

# Split the data set by `species` and make sure all curves have similar numbers
# of points (within 3 of the mean value); this will intentionally fail
check_response_curve_data(
  licor_file,
  'species',
  expected_npts = c(0, 3),
  error_on_failure = FALSE
)
#>   species npts
#> 1 soybean   21
#> 2 tobacco    7
#> 
#> The following curves have atypical numbers of points (more than 3 away from the mean value of 14): soybean, tobacco

# Split the data set by `species` and make sure all curves have a constant value
# of `plot` and a limited range of `TLeafCnd`; this will intentionally fail
check_response_curve_data(
  licor_file,
  'species',
  constant_col = list(plot = NA, TleafCnd = 0.001),
  error_on_failure = FALSE
)
#>   species npts
#> 1 soybean   21
#> 2 tobacco    7
#> 
#> Not all curves have the same number of points
#> 
#> Curve ID: species = soybean:
#>     The `plot` column takes 3 values: `1a`, `5`, `1b`
#>     The `TleafCnd` column range (4.0666) exceeds the limit (0.001)
#> Curve ID: species = tobacco:
#>     The `TleafCnd` column range (1.8677) exceeds the limit (0.001)
#> 
#> 
#> One or more curves has non-constant values of a column that should be constant.
```
