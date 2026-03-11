# Exclude outliers from a data set

Excludes outliers from a data set using the "1.5 interquartile range"
rule.

## Usage

``` r
exclude_outliers(x, col_for_analysis, INDICES, method = 'exclude')

  # S3 method for class 'data.frame'
exclude_outliers(x, col_for_analysis, INDICES, method = 'exclude')

  # S3 method for class 'exdf'
exclude_outliers(x, col_for_analysis, INDICES, method = 'exclude')
```

## Arguments

- x:

  A data table

- col_for_analysis:

  The name of a column of `x` that should be used to determine outliers.

- INDICES:

  A factor or list of factors that each `nrow(x)` elements.

- method:

  Specify whether to remove rows from `x` (`'remove'`) or to replace
  outlier values of `col_for_analysis` with `NA` (`'exclude'`).

## Value

This function returns an object formed from `x`, where the results
depend on on the value of `method`.

When `method` is `'remove'`, the returned object is a modified copy of
`x` where all rows in which the value of `col_for_analysis` is an
outlier have been removed.

When `method` is `'exclude'`, the returned object is a modified copy of
`x` where all outlier values of `col_for_analysis` have been replaced
with `NA`.

## Details

`exclude_outliers` is generic, with methods defined for data frames and
`exdf` objects. This function uses a simple rule to detect outliers,
where any point that deviates from the mean by more than `1.5 * IQR`,
where `IQR` is the interquartile range, is said to be an outlier. This
method is also sometimes referred to as "Tukey's Fences," as seen in the
[Wikipedia page about outliers](https://en.wikipedia.org/wiki/Outlier).

For data sets with extreme outliers, it may be necessary to exclude
outliers more than once to actually remove them all.

## See also

[`exdf`](https://eloch216.github.io/PhotoGEA/reference/exdf.md)

## Examples

``` r
# Read a Licor file included with the PhotoGEA package; this file includes
# several light response curves that can be identified by the 'species' and
# 'plot' columns.
licor_file <- read_gasex_file(
  PhotoGEA_example_file_path('ball_berry_1.xlsx')
)

# Remove points from each response curve in the data where the leaf temperature
# is determined to be an outlier
licor_file_clean <- exclude_outliers(
  licor_file,
  'TleafCnd',
  list(licor_file[, 'species'], licor_file[, 'plot']),
  method = 'remove'
)

# Check to see how many points remain after removing outliers
str(list('original' = nrow(licor_file), 'clean' = nrow(licor_file_clean)))
#> List of 2
#>  $ original: int 28
#>  $ clean   : int 24
```
