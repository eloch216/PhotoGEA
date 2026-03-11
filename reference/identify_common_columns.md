# Identify columns that are common to multiple objects

Checks whether the input arguments have the same columns

## Usage

``` r
identify_common_columns(...)

  # S3 method for class 'data.frame'
identify_common_columns(...)

  # S3 method for class 'exdf'
identify_common_columns(...)
```

## Arguments

- ...:

  One or more R objects that have column names.

## Value

A character vector of the column names that are common to all the input
objects.

## Details

`identify_common_columns` is generic, with methods defined for data
frames and `exdf` objects. In the case of `exdf` objects, a column will
only be considered common if it has the same name, units, and category
in all of the input objects.

## See also

[`exdf`](https://eloch216.github.io/PhotoGEA/reference/exdf.md)

## Examples

``` r
# Here we create two exdf objects with the same column names and units, but
# where the categories of one column are not the same in both objects
exdf_1 <- exdf(
  data.frame(A = c(3, 2, 7, 9), B = c(4, 5, 1, 8)),
  data.frame(A = 'm', B = 's'),
  data.frame(A = 'Cat1', B = 'Cat2')
)

exdf_2 <- exdf(
  data.frame(A = c(3, 2, 7, 9), B = c(4, 5, 1, 8)),
  data.frame(A = 'm', B = 's'),
  data.frame(A = 'Cat1', B = 'Cat3')
)

# Calling `identify_common_columns` on the exdf objects will only identify one
# common column (A) because the category for column B is not common to all the
# exdf objects.
identify_common_columns(exdf_1, exdf_2)
#> [1] "A"

# Calling `identify_common_columns` on the main_data data frames will identify
# two common columns because unit and category information will not be
# considered here.
identify_common_columns(exdf_1$main_data, exdf_2$main_data)
#> [1] "A" "B"
```
