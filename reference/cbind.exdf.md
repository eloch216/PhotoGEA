# Combine exdf objects by columns or rows

Combines one or more `exdf` objects by the columns or rows of their
`main_data`. For `rbind`, errors will occur if column names are not the
same in all of the `exdf` objects, and if all units and categories are
not identical.

## Usage

``` r
# S3 method for class 'exdf'
cbind(..., deparse.level = 1)

  # S3 method for class 'exdf'
rbind(
    ...,
    deparse.level = 1,
    make.row.names = TRUE,
    stringsAsFactors = FALSE
  )
```

## Arguments

- ...:

  Two or more `exdf` objects.

- deparse.level:

  See associated documentation for the generic versions of
  [`cbind`](https://rdrr.io/r/base/cbind.html) and
  [`rbind`](https://rdrr.io/r/base/cbind.html).

- make.row.names:

  See associated documentation for the generic version of `rbind`.

- stringsAsFactors:

  See associated documentation for the generic version of `rbind`.

## Value

Returns a new `exdf` object.

## See also

[`exdf`](https://eloch216.github.io/PhotoGEA/reference/exdf.md)

## Examples

``` r
# Make some simple exdf objects. 1 and 2 have the same number of rows but
# different columns, while 1 and 3 have the same columns but different rows.
simple_exdf_1 <- exdf(data.frame(A = 1), data.frame(A = 'au'), data.frame(A = 'ac'))
simple_exdf_2 <- exdf(data.frame(B = 2), data.frame(B = 'bu'), data.frame(B = 'bc'))
simple_exdf_3 <- exdf(data.frame(A = 2), data.frame(A = 'au'), data.frame(A = 'ac'))

cbind(simple_exdf_1) # will just return simple_exdf_1
#> 
#> Converting an `exdf` object to a `data.frame` before printing
#> 
#>   A [ac] (au)
#> 1           1
cbind(simple_exdf_1, simple_exdf_2)
#> 
#> Converting an `exdf` object to a `data.frame` before printing
#> 
#>   A [ac] (au) B [bc] (bu)
#> 1           1           2

rbind(simple_exdf_1) # will just return simple_exdf_1
#> 
#> Converting an `exdf` object to a `data.frame` before printing
#> 
#>   A [ac] (au)
#> 1           1
rbind(simple_exdf_1, simple_exdf_3)
#> 
#> Converting an `exdf` object to a `data.frame` before printing
#> 
#>   A [ac] (au)
#> 1           1
#> 2           2
```
