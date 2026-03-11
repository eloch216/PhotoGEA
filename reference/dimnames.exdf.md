# Retrieve or set the dimension names of an exdf object

Returns or sets the dimension names of an `exdf` object's `main_data`.
When setting names, the column names of the `exdf` object's `units` and
`categories` are also set. Also enables `colnames` and `rownames` for
`exdf` objects.

## Usage

``` r
# S3 method for class 'exdf'
dimnames(x)

  # S3 method for class 'exdf'
dimnames(x) <- value
```

## Arguments

- x:

  An `exdf` object.

- value:

  A possible value for `dimnames(x)`

## Value

Returns `dimnames(x[['main_data']])`.

## See also

[`exdf`](https://eloch216.github.io/PhotoGEA/reference/exdf.md)

## Examples

``` r
simple_exdf <- exdf(data.frame(A = 1), data.frame(A = 'u'), data.frame(A = 'c'))

dimnames(simple_exdf)
#> [[1]]
#> [1] "1"
#> 
#> [[2]]
#> [1] "A"
#> 
dimnames(simple_exdf[['main_data']]) # An equivalent command
#> [[1]]
#> [1] "1"
#> 
#> [[2]]
#> [1] "A"
#> 

colnames(simple_exdf) <- "B"
rownames(simple_exdf) <- 2

colnames(simple_exdf)
#> [1] "B"
rownames(simple_exdf)
#> [1] "2"
```
