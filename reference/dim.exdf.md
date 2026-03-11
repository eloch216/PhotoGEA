# Retrieve the dimension of an exdf object

Returns the dimensions of an `exdf` object's `main_data`. Also enables
`nrow` and `ncol` for `exdf` objects.

## Usage

``` r
# S3 method for class 'exdf'
dim(x)
```

## Arguments

- x:

  An `exdf` object.

## Value

Returns `dim(x[['main_data']])`.

## See also

[`exdf`](https://eloch216.github.io/PhotoGEA/reference/exdf.md)

## Examples

``` r
simple_exdf <- exdf(data.frame(A = 1), data.frame(A = 'u'), data.frame(A = 'c'))

dim(simple_exdf)
#> [1] 1 1
dim(simple_exdf[['main_data']]) # An equivalent command
#> [1] 1 1

nrow(simple_exdf)
#> [1] 1
ncol(simple_exdf)
#> [1] 1
```
