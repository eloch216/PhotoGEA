# Length of an exdf object

Returns the length of an `exdf` object's `main_data`.

## Usage

``` r
# S3 method for class 'exdf'
length(x)
```

## Arguments

- x:

  An `exdf` object.

## Value

Returns `length(x[['main_data']])`.

## See also

[`exdf`](https://eloch216.github.io/PhotoGEA/reference/exdf.md)

## Examples

``` r
simple_exdf <- exdf(data.frame(A = 1), data.frame(A = 'u'), data.frame(A = 'c'))
length(simple_exdf)
#> [1] 1
length(simple_exdf[['main_data']]) # An equivalent command
#> [1] 1
```
