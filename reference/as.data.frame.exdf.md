# Convert an exdf object to a data frame

Converts an `exdf` object to a data frame by appending the units and
categories to the top of each column in the `exdf` object's `main_data`
data frame. Typically this function is used for displaying the contents
of an `exdf` object; in fact, it is used internally by `View`,
`write.csv`, and other functions. The `main_data` of an `exdf` object
`x` can be accessed directly (without including the units and categories
in the first row) via `x[['main_data']]` as with any other list element.

## Usage

``` r
# S3 method for class 'exdf'
as.data.frame(x, ...)
```

## Arguments

- x:

  An `exdf` object.

- ...:

  Unused.

## Value

A data frame formed from `x`.

## See also

[`exdf`](https://eloch216.github.io/PhotoGEA/reference/exdf.md)

## Examples

``` r
simple_exdf <- exdf(data.frame(A = 1), data.frame(A = 'u'), data.frame(A = 'c'))
as.data.frame(simple_exdf) # Includes units and categories in the first rows
#>   A
#> 1 c
#> 2 u
#> 3 1
simple_exdf[['main_data']] # Just returns the main data
#>   A
#> 1 1
```
