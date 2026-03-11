# Remove specific points from an exdf object

Removes all points from an `exdf` object that satisfy a set of
conditions.

## Usage

``` r
remove_points(exdf_obj, ..., method = 'remove')
```

## Arguments

- exdf_obj:

  An `exdf` object.

- ...:

  Each optional argument should be a list of named elements that specify
  points to be removed from `exdf_obj`. For example,
  `list(species = 'soybean', plot = c('1a', '1b'))` specifies the set of
  points where (1) `species` is `'soybean'` and (2) `plot` is `'1a'` or
  `'1b'`.

- method:

  Specify whether to remove points (`'remove'`) or designate them as
  being excluded from subsequent fits (`'exclude'`); see below for more
  details.

## Value

This function returns an `exdf` object formed from `exdf_obj`, where the
result depends on the value of `method`.

When `method` is `'remove'`, the returned object is a modified copy of
`exdf_obj` where all rows that meet the conditions specified by the
optional arguments have been removed.

When `method` is `'exclude'`, the returned object is a modified copy of
`exdf_obj` with a new column called `include_when_fitting`. The value of
this column is `FALSE` for all rows that meet the conditions specified
by the optional arguments, and `TRUE` otherwise. Points where this
column is `FALSE` will not be used for fitting by
[`fit_c3_aci`](https://eloch216.github.io/PhotoGEA/reference/fit_c3_aci.md)
or other fitting functions.

## See also

[`exdf`](https://eloch216.github.io/PhotoGEA/reference/exdf.md)

## Examples

``` r
# Create an exdf object by reading a Licor Excel file
licor_file <- read_gasex_file(
  PhotoGEA_example_file_path('ball_berry_1.xlsx')
)

# Print the number of points in the data set
nrow(licor_file)
#> [1] 28

# Remove the following:
# - All points where `obs` is 28 (1 point)
# - All points where `species` is `soybean` and `plot` is `1a` or `1b` (14 points)
licor_file_2 <- remove_points(
  licor_file,
  list(obs = 28),
  list(species = 'soybean', plot = c('1a', '1b')),
  method = 'remove'
)

# There should now be 15 fewer points remaining in the data set
nrow(licor_file_2)
#> [1] 13

# We can also specify the same points for exclusion rather than removal:
licor_file_3 <- remove_points(
  licor_file,
  list(obs = 28),
  list(species = 'soybean', plot = c('1a', '1b')),
  method = 'exclude'
)

print(licor_file_3[, c('species', 'plot', 'include_when_fitting')])
#>    species plot include_when_fitting
#> 1  soybean   1a                FALSE
#> 2  soybean   1a                FALSE
#> 3  soybean   1a                FALSE
#> 4  soybean   1a                FALSE
#> 5  soybean   1a                FALSE
#> 6  soybean   1a                FALSE
#> 7  soybean   1a                FALSE
#> 8  soybean    5                 TRUE
#> 9  soybean    5                 TRUE
#> 10 soybean    5                 TRUE
#> 11 soybean    5                 TRUE
#> 12 soybean    5                 TRUE
#> 13 soybean    5                 TRUE
#> 14 soybean    5                 TRUE
#> 15 soybean   1b                FALSE
#> 16 soybean   1b                FALSE
#> 17 soybean   1b                FALSE
#> 18 soybean   1b                FALSE
#> 19 soybean   1b                FALSE
#> 20 soybean   1b                FALSE
#> 21 soybean   1b                FALSE
#> 22 tobacco    2                 TRUE
#> 23 tobacco    2                 TRUE
#> 24 tobacco    2                 TRUE
#> 25 tobacco    2                 TRUE
#> 26 tobacco    2                 TRUE
#> 27 tobacco    2                 TRUE
#> 28 tobacco    2                FALSE

# The number of points where `include_when_fitting` is TRUE should be the same
# as the number of remaining rows when using the `remove` method
sum(licor_file_3[, 'include_when_fitting'])
#> [1] 13
```
