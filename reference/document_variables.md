# Document exdf columns by specifying units and categories

Adds new columns to a table-like object, and sets/modifies the units or
categories of columns in an `exdf` object.

## Usage

``` r
document_variables(x, ...)

  # S3 method for class 'data.frame'
document_variables(x, ...)

  # S3 method for class 'exdf'
document_variables(x, ...)
```

## Arguments

- x:

  A table-like R object such as a data frame or an `exdf`.

- ...:

  Each optional argument should be a character vector with three
  elements that describe a column, where the first element is the
  category, the second is the name, and the third is the units. For
  example, `c('GasEx', 'A', 'micromol m^(-2) s^(-1)')` specifies that
  the category and units for the `A` column are `GasEx` and
  `micromol m^(-2) s^(-1)`, respectively. If the column name is not in
  `x`, it will be added with all values initialized to `NA`. Categories
  and units will be ignored when `x` is a data frame.

## Value

An object based on `x` with new and/or modified columns.

## See also

[`exdf`](https://eloch216.github.io/PhotoGEA/reference/exdf.md)

## Examples

``` r
# Create a simple exdf object with two columns (`A` and `B`) and default values
# for its units and categories.
simple_exdf <- exdf(data.frame(A = c(3, 2, 7, 9), B = c(4, 5, 1, 8)))

print(simple_exdf)
#> 
#> Converting an `exdf` object to a `data.frame` before printing
#> 
#>   A [NA] (NA) B [NA] (NA)
#> 1           3           4
#> 2           2           5
#> 3           7           1
#> 4           9           8

# Specify units and categories for the `A` and `B` columns, and add a new `C`
# column.
document_variables(
  simple_exdf,
  c('cat1', 'A', 'm'), # The category of `A` is `cat1` and its units are `m`
  c('cat2', 'B', 's'), # The category of `B` is `cat2` and its units are `s`
  c('cat3', 'C', 'g')  # The category of `C` is `cat3` and its units are `g`
)
#> 
#> Converting an `exdf` object to a `data.frame` before printing
#> 
#>   A [cat1] (m) B [cat2] (s) C [cat3] (g)
#> 1            3            4           NA
#> 2            2            5           NA
#> 3            7            1           NA
#> 4            9            8           NA

# Do the same but for a data frame; in this case columns A and B will not be
# altered, but a new column C will be added (and initialized to NA)
document_variables(
  simple_exdf$main_data,
  c('cat1', 'A', 'm'), # The category of `A` is `cat1` and its units are `m`
  c('cat2', 'B', 's'), # The category of `B` is `cat2` and its units are `s`
  c('cat3', 'C', 'g')  # The category of `C` is `cat3` and its units are `g`
)
#>   A B  C
#> 1 3 4 NA
#> 2 2 5 NA
#> 3 7 1 NA
#> 4 9 8 NA
```
