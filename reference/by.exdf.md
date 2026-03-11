# Apply a function to an exdf object split by one or more factors

Divides an `exdf` object into groups defined by one or more factors and
applies a function to each group.

## Usage

``` r
# S3 method for class 'exdf'
by(data, INDICES, FUN, ...)
```

## Arguments

- data:

  An `exdf` object.

- INDICES:

  A factor or a list of factors.

- FUN:

  A function whose first input argument is an `exdf` object.

- ...:

  Additional arguments to be passed to `FUN`.

## Value

Splits `data` into chunks `x` by the values of the `INDICES` and calls
`FUN(x, ...)` for each chunk; returns a list where each element is the
output from each call to `FUN`.

## See also

[`exdf`](https://eloch216.github.io/PhotoGEA/reference/exdf.md)

## Examples

``` r
# Read a Licor file, split it into chunks according to the `species` column,
# and count the number of measurements for each species
licor_file <- read_gasex_file(
  PhotoGEA_example_file_path('ball_berry_1.xlsx')
)

by(licor_file, licor_file[, 'species'], nrow)
#> $soybean
#> [1] 21
#> 
#> $tobacco
#> [1] 7
#> 
```
