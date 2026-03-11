# Identify C3 Limiting Processes

Identify limiting processes in a C3 curve, typically the result of a
fit. It is rate for users to call this function directly because it is
used internally by
[`fit_c3_aci`](https://eloch216.github.io/PhotoGEA/reference/fit_c3_aci.md)
and
[`fit_c3_variable_j`](https://eloch216.github.io/PhotoGEA/reference/fit_c3_variable_j.md).

## Usage

``` r
identify_c3_limiting_processes(
    data_table,
    a_column_name = 'A_fit',
    ac_column_name = 'Ac',
    aj_column_name = 'Aj',
    ap_column_name = 'Ap',
    tol = 1e-3
)
```

## Arguments

- data_table:

  A table-like R object such as a data frame or an `exdf`.

- a_column_name:

  The name of the column in `data_table` that contains the modeled net
  CO2 assimilation rate in `micromol m^(-2) s^(-1)`.

- ac_column_name:

  The name of the column in `data_table` that contains the modeled
  Rubisco-limited net CO2 assimilation rate in `micromol m^(-2) s^(-1)`.

- aj_column_name:

  The name of the column in `data_table` that contains the modeled
  RuBP-regeneration-limited net CO2 assimilation rate in
  `micromol m^(-2) s^(-1)`.

- ap_column_name:

  The name of the column in `data_table` that contains the modeled
  TPU-limited net CO2 assimilation rate in `micromol m^(-2) s^(-1)`.

- tol:

  A relative tolerance factor used to identify when two rates are equal.

## Details

For a C3 leaf, `An` is given by either `Ac`, `Aj`, or `Ap`. See the
documentation for
[`calculate_c3_assimilation`](https://eloch216.github.io/PhotoGEA/reference/calculate_c3_assimilation.md)
for more information.

This function first identifies points where `An = Ac`, `An = Aj`, and
`An = Ap`. The results are stored in columns called `Ac_limiting`,
`Aj_limiting`, and `Ap_limiting`, where a value of `TRUE` indicates that
the corresponding process is limiting.

Then, the overall limiting state is specified in the `limiting_process`
column. For example, points where `An` equals `Ac` but not `Aj` or `Ap`
are designated by `limiting_process = 'Ac'`, and likewise for the other
potential limiting processes. If more than one process is limiting for a
point, `limiting_process` is set to `'co-limited'`.

## Value

An `exdf` object based on `licor_exdf` that includes new columns as
described above: `Ac_limiting`, `Aj_limiting`, `Ap_limiting`, and
`limiting_process`. The categories of these new columns are set to
`identify_c3_limiting_processes` to indicate that they were created
using this function.

## Examples

``` r
# Identify limiting processes in an example curve
example_curve <- exdf(
  data.frame(
    A_fit = c(1.0, 2.0, 3.0, 4.0, 4.0),
    Ac    = c(1.0, 2.0, 5.0, 8.0, 9.0),
    Aj    = c(2.0, 2.5, 3.0, 4.0, 8.0),
    Ap    = c(NA,  NA,  4.0, 4.0, 4.0)
  ),
  units = data.frame(
    A_fit = 'micromol m^(-2) s^(-1)',
    Ac    = 'micromol m^(-2) s^(-1)',
    Aj    = 'micromol m^(-2) s^(-1)',
    Ap    = 'micromol m^(-2) s^(-1)',
    stringsAsFactors = FALSE
  )
)

identify_c3_limiting_processes(example_curve)
#> 
#> Converting an `exdf` object to a `data.frame` before printing
#> 
#>   A_fit [NA] (micromol m^(-2) s^(-1)) Ac [NA] (micromol m^(-2) s^(-1))
#> 1                                   1                                1
#> 2                                   2                                2
#> 3                                   3                                5
#> 4                                   4                                8
#> 5                                   4                                9
#>   Aj [NA] (micromol m^(-2) s^(-1)) Ap [NA] (micromol m^(-2) s^(-1))
#> 1                              2.0                               NA
#> 2                              2.5                               NA
#> 3                              3.0                                4
#> 4                              4.0                                4
#> 5                              8.0                                4
#>   Ac_limiting [identify_c3_limits] () Aj_limiting [identify_c3_limits] ()
#> 1                                TRUE                               FALSE
#> 2                                TRUE                               FALSE
#> 3                               FALSE                                TRUE
#> 4                               FALSE                                TRUE
#> 5                               FALSE                               FALSE
#>   Ap_limiting [identify_c3_limits] () limiting_process [identify_c3_limits] ()
#> 1                               FALSE                                       Ac
#> 2                               FALSE                                       Ac
#> 3                               FALSE                                       Aj
#> 4                                TRUE                   co-limited (Aj and Ap)
#> 5                                TRUE                                       Ap

# This function also works for data frames
identify_c3_limiting_processes(example_curve$main_data)
#>   A_fit Ac  Aj Ap Ac_limiting Aj_limiting Ap_limiting       limiting_process
#> 1     1  1 2.0 NA        TRUE       FALSE       FALSE                     Ac
#> 2     2  2 2.5 NA        TRUE       FALSE       FALSE                     Ac
#> 3     3  5 3.0  4       FALSE        TRUE       FALSE                     Aj
#> 4     4  8 4.0  4       FALSE        TRUE        TRUE co-limited (Aj and Ap)
#> 5     4  9 8.0  4       FALSE       FALSE        TRUE                     Ap
```
