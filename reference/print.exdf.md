# Print the contents of an exdf object

Prints the contents of an `exdf` object's `main_data`. Each column is
described by its name, unit, and category formatted like
`name [category] (units)`.

## Usage

``` r
# S3 method for class 'exdf'
print(x, ...)
```

## Arguments

- x:

  An `exdf` object.

- ...:

  Additional arguments to be passed to `print`.

## Value

None.

## See also

[`exdf`](https://eloch216.github.io/PhotoGEA/reference/exdf.md)

## Examples

``` r
simple_exdf <- exdf(data.frame(A = 1), data.frame(A = 'u'), data.frame(A = 'c'))
print(simple_exdf)
#> 
#> Converting an `exdf` object to a `data.frame` before printing
#> 
#>   A [c] (u)
#> 1         1
```
