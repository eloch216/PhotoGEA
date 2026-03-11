# Divide an exdf object into groups

Divides an `exdf` object into groups defined by one or more factors.

## Usage

``` r
# S3 method for class 'exdf'
split(x, f, drop = FALSE, lex.order = FALSE, ...)
```

## Arguments

- x:

  An `exdf` object.

- f:

  A factor or a list of factors.

- drop:

  A logical value indicating whether levels of `f` that do not occur
  should be dropped.

- lex.order:

  A logical value passed to
  [`interaction`](https://rdrr.io/r/base/interaction.html).

- ...:

  Additional arguments to be passed to the default method of
  [`split`](https://rdrr.io/r/base/split.html).

## Value

Returns a list of `exdf` objects created by splitting `x` along the
values of `f`.

## See also

[`exdf`](https://eloch216.github.io/PhotoGEA/reference/exdf.md)

## Examples

``` r
# Read a Licor file, select just a few columns, and then split it by the value
# of the `plot` column
licor_file <- read_gasex_file(
  PhotoGEA_example_file_path('ball_berry_1.xlsx')
)

licor_file <- licor_file[, c('plot', 'species', 'Qin', 'A', 'gsw'), TRUE]

split(
  licor_file,
  list(licor_file[,'species'], licor_file[,'plot']),
  drop = TRUE
)
#> $soybean.1a
#> 
#> Converting an `exdf` object to a `data.frame` before printing
#> 
#>   plot [UserDefCon] (NA) species [UserDefCon] (NA)
#> 1                     1a                   soybean
#> 2                     1a                   soybean
#> 3                     1a                   soybean
#> 4                     1a                   soybean
#> 5                     1a                   soybean
#> 6                     1a                   soybean
#> 7                     1a                   soybean
#>   Qin [LeafQ] (micromol m^(-2) s^(-1)) A [GasEx] (micromol m^(-2) s^(-1))
#> 1                             2000.030                          35.400473
#> 2                             1099.800                          22.996635
#> 3                              799.983                          25.430106
#> 4                              649.951                          23.290614
#> 5                              500.135                          16.617929
#> 6                              350.106                          12.240074
#> 7                              200.141                           5.679174
#>   gsw [GasEx] (mol m^(-2) s^(-1))
#> 1                       0.3881885
#> 2                       0.3044869
#> 3                       0.2804279
#> 4                       0.2181023
#> 5                       0.2123801
#> 6                       0.1848656
#> 7                       0.1685880
#> 
#> $soybean.1b
#> 
#> Converting an `exdf` object to a `data.frame` before printing
#> 
#>    plot [UserDefCon] (NA) species [UserDefCon] (NA)
#> 15                     1b                   soybean
#> 16                     1b                   soybean
#> 17                     1b                   soybean
#> 18                     1b                   soybean
#> 19                     1b                   soybean
#> 20                     1b                   soybean
#> 21                     1b                   soybean
#>    Qin [LeafQ] (micromol m^(-2) s^(-1)) A [GasEx] (micromol m^(-2) s^(-1))
#> 15                             1999.850                          47.950460
#> 16                             1099.890                          30.250627
#> 17                              799.873                          25.472460
#> 18                              650.070                          26.383035
#> 19                              499.982                          17.009070
#> 20                              350.060                          14.666531
#> 21                              200.070                           9.276031
#>    gsw [GasEx] (mol m^(-2) s^(-1))
#> 15                       0.6974767
#> 16                       0.5627892
#> 17                       0.4941281
#> 18                       0.4676789
#> 19                       0.4235752
#> 20                       0.1587479
#> 21                       0.2030845
#> 
#> $tobacco.2
#> 
#> Converting an `exdf` object to a `data.frame` before printing
#> 
#>    plot [UserDefCon] (NA) species [UserDefCon] (NA)
#> 22                      2                   tobacco
#> 23                      2                   tobacco
#> 24                      2                   tobacco
#> 25                      2                   tobacco
#> 26                      2                   tobacco
#> 27                      2                   tobacco
#> 28                      2                   tobacco
#>    Qin [LeafQ] (micromol m^(-2) s^(-1)) A [GasEx] (micromol m^(-2) s^(-1))
#> 22                             2000.020                           29.61454
#> 23                             1099.980                           25.49228
#> 24                              800.046                           22.21817
#> 25                              650.001                           16.34589
#> 26                              500.036                           12.10635
#> 27                              350.073                           12.04469
#> 28                              200.236                            7.84987
#>    gsw [GasEx] (mol m^(-2) s^(-1))
#> 22                      0.29002835
#> 23                      0.26095999
#> 24                      0.24715976
#> 25                      0.22395230
#> 26                      0.22348130
#> 27                      0.20608729
#> 28                      0.08301745
#> 
#> $soybean.5
#> 
#> Converting an `exdf` object to a `data.frame` before printing
#> 
#>    plot [UserDefCon] (NA) species [UserDefCon] (NA)
#> 8                       5                   soybean
#> 9                       5                   soybean
#> 10                      5                   soybean
#> 11                      5                   soybean
#> 12                      5                   soybean
#> 13                      5                   soybean
#> 14                      5                   soybean
#>    Qin [LeafQ] (micromol m^(-2) s^(-1)) A [GasEx] (micromol m^(-2) s^(-1))
#> 8                              1999.980                           45.78075
#> 9                              1099.960                           36.21825
#> 10                              799.910                           32.76940
#> 11                              650.108                           22.92448
#> 12                              499.921                           22.06819
#> 13                              349.991                           17.09159
#> 14                              200.029                           10.76825
#>    gsw [GasEx] (mol m^(-2) s^(-1))
#> 8                        0.8771680
#> 9                        0.6695119
#> 10                       0.5599051
#> 11                       0.5758128
#> 12                       0.3805223
#> 13                       0.3462828
#> 14                       0.3181253
#> 
```
