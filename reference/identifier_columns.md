# Find columns that have a single value across all rows

Identifies columns that have a single value across all rows and returns
them.

## Usage

``` r
identifier_columns(x)

  # S3 method for class 'data.frame'
identifier_columns(x)

  # S3 method for class 'exdf'
identifier_columns(x)
```

## Arguments

- x:

  A table-like R object such as a data frame or an `exdf`.

## Details

`identifier_columns` is generic, with methods defined for data frames
and `exdf` objects.

`identifier_columns` gets the names and values of any columns in a
table-like object that have a single unique value. If the object
represents a set of data from one replicate, then these special columns
are taken to be "identifiers" that describe the replicate. This function
is often used inside fitting functions that are passed to
[`by.exdf`](https://eloch216.github.io/PhotoGEA/reference/by.exdf.md) as
its `FUN` input argument. For example, see the code for
[`fit_ball_berry`](https://eloch216.github.io/PhotoGEA/reference/fit_ball_berry.md)
by typing
[`PhotoGEA::fit_ball_berry`](https://eloch216.github.io/PhotoGEA/reference/fit_ball_berry.md)
in the R terminal.

## Value

The return value will be a subset of `x`, restricted to only include
columns whose values are constant. Only one row will be returned.

## See also

[`exdf`](https://eloch216.github.io/PhotoGEA/reference/exdf.md)

## Examples

``` r
# Create a simple exdf object
simple_exdf <- exdf(
  data.frame(A = c(3, 2, 7, 9), species = c('a', 'a', 'a', 'a'), plot = c(1, 1, 1, 1)),
  data.frame(A = 'm', species = '', plot = ''),
  data.frame(A = 'Cat1', species = '', plot = '')
)

# Find its identifier columns
identifier_columns(simple_exdf)
#> 
#> Converting an `exdf` object to a `data.frame` before printing
#> 
#>   species [] () plot [] ()
#> 1             a          1

# Apply the data frame method to the exdf object's main data frame
identifier_columns(simple_exdf$main_data)
#>   species plot
#> 1       a    1
```
