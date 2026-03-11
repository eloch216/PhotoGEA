# Is an object an exdf?

Checks whether an object is an `exdf` object.

## Usage

``` r
is.exdf(x, consistency_check = FALSE)
```

## Arguments

- x:

  An R object.

- consistency_check:

  A logical value indicating whether to perform additional consistency
  checks.

## Details

The default version of `is.exdf` simply checks to see if `'exdf'` is in
`class(x)`.

If `consistency_check` is `TRUE`, then additional checks will be
performed to make sure the object has three elements named `main_data`,
`units`, and `categories`; that these elements are data frames with the
same column names; and that `units` and `categories` each have one row.
These requirements are all part of the definition of an `exdf` object,
but these checks require additional time so they are not always desired.

## Value

A logical (TRUE / FALSE) value indicating whether the object is an
`exdf` object.

## See also

[`exdf`](https://eloch216.github.io/PhotoGEA/reference/exdf.md)

## Examples

``` r
# Test a simple exdf object
simple_exdf <- exdf(data.frame(A = 1), data.frame(A = 'u'), data.frame(A = 'c'))
is.exdf(simple_exdf)
#> [1] TRUE
is.exdf(simple_exdf, TRUE)
#> [1] TRUE

# Test an object that is clearly not an exdf
not_an_exdf <- 2
is.exdf(not_an_exdf)
#> [1] FALSE
is.exdf(not_an_exdf, TRUE)
#> [1] FALSE

# Test an object that claims to be an exdf but does not meet all of the
# requirements
fake_exdf <- not_an_exdf
class(fake_exdf) <- c('exdf', class(fake_exdf))
is.exdf(fake_exdf)
#> [1] TRUE
is.exdf(fake_exdf, TRUE)
#> Warning: x must have elements named `main_data`, `units`, and `categories`
#> [1] FALSE
```
