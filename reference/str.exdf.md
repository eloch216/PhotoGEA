# Display the structure of an exdf object

Displays the structure of an `exdf` object's `main_data`. Each column is
described by its name, unit, and category formatted like
`name [category] (units)`.

## Usage

``` r
# S3 method for class 'exdf'
str(object, ...)
```

## Arguments

- object:

  An `exdf` object.

- ...:

  Additional arguments to be passed to `str`.

## Value

None.

## See also

[`exdf`](https://eloch216.github.io/PhotoGEA/reference/exdf.md)

## Examples

``` r
simple_exdf <- exdf(data.frame(A = 1), data.frame(A = 'u'), data.frame(A = 'c'))
str(simple_exdf)
#> 
#> Converting an `exdf` object to a `data.frame` before printing
#> 
#> 'data.frame':    1 obs. of  1 variable:
#>  $ A [c] (u): num 1
```
