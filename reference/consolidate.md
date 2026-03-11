# Consolidate a list of lists

Consolidates a list of lists into a regular list by combining like-named
elements.

## Usage

``` r
consolidate(x)

  # S3 method for class 'data.frame'
consolidate(x)

  # S3 method for class 'exdf'
consolidate(x)
```

## Arguments

- x:

  A list of lists `list_1`, `list_2`, ..., `list_N`, where each sub-list
  `list_i` has elements named `name_1`, `name_2`, ..., `name_M`.

## Value

A list with elements named `name_1`, `name_2`, ..., `name_M`, where each
element was created by combining all elements of `x` with the same name
using `rbind`; for example, the element with name `name_1` will be
created by calling
`rbind(list_1$name_1, list_2$name_1, ..., list_N$name_1)`. Before
calling `rbind`, each element will be limited to the columns that are
common to all elements with the same name.

## Details

`consolidate` is generic, with methods defined for nested lists of data
frames and `exdf` objects.

## See also

[`exdf`](https://eloch216.github.io/PhotoGEA/reference/exdf.md)

## Examples

``` r
# Example 1: Create a nested list of data frames and then consolidate them into
# a regular list by combining the like-named elements
nested_df_list <- list(
  list_1 = list(
    name_1 = data.frame(A = c(1, 2), B = c(0, 0)),
    name_2 = data.frame(A = c(3, 4), B = c(0, 0)),
    name_3 = data.frame(A = c(5, 6), B = c(0, 0))
  ),
  list_2 = list(
    name_1 = data.frame(A = c(7, 8), B = c(0, 0)),
    name_2 = data.frame(A = c(9, 10), B = c(0, 0)),
    name_3 = data.frame(A = c(11, 12), B = c(0, 0))
  ),
  list_3 = list(
    name_1 = data.frame(A = c(13, 14), B = c(0, 0)),
    name_2 = data.frame(A = c(15, 16), B = c(0, 0)),
    name_3 = data.frame(A = c(17, 18), B = c(0, 0))
  )
)

str(nested_df_list)
#> List of 3
#>  $ list_1:List of 3
#>   ..$ name_1:'data.frame':   2 obs. of  2 variables:
#>   .. ..$ A: num [1:2] 1 2
#>   .. ..$ B: num [1:2] 0 0
#>   ..$ name_2:'data.frame':   2 obs. of  2 variables:
#>   .. ..$ A: num [1:2] 3 4
#>   .. ..$ B: num [1:2] 0 0
#>   ..$ name_3:'data.frame':   2 obs. of  2 variables:
#>   .. ..$ A: num [1:2] 5 6
#>   .. ..$ B: num [1:2] 0 0
#>  $ list_2:List of 3
#>   ..$ name_1:'data.frame':   2 obs. of  2 variables:
#>   .. ..$ A: num [1:2] 7 8
#>   .. ..$ B: num [1:2] 0 0
#>   ..$ name_2:'data.frame':   2 obs. of  2 variables:
#>   .. ..$ A: num [1:2] 9 10
#>   .. ..$ B: num [1:2] 0 0
#>   ..$ name_3:'data.frame':   2 obs. of  2 variables:
#>   .. ..$ A: num [1:2] 11 12
#>   .. ..$ B: num [1:2] 0 0
#>  $ list_3:List of 3
#>   ..$ name_1:'data.frame':   2 obs. of  2 variables:
#>   .. ..$ A: num [1:2] 13 14
#>   .. ..$ B: num [1:2] 0 0
#>   ..$ name_2:'data.frame':   2 obs. of  2 variables:
#>   .. ..$ A: num [1:2] 15 16
#>   .. ..$ B: num [1:2] 0 0
#>   ..$ name_3:'data.frame':   2 obs. of  2 variables:
#>   .. ..$ A: num [1:2] 17 18
#>   .. ..$ B: num [1:2] 0 0

consolidated_df_list <- consolidate(nested_df_list)

str(consolidated_df_list)
#> List of 3
#>  $ name_1:'data.frame':  6 obs. of  2 variables:
#>   ..$ A: num [1:6] 1 2 7 8 13 14
#>   ..$ B: num [1:6] 0 0 0 0 0 0
#>  $ name_2:'data.frame':  6 obs. of  2 variables:
#>   ..$ A: num [1:6] 3 4 9 10 15 16
#>   ..$ B: num [1:6] 0 0 0 0 0 0
#>  $ name_3:'data.frame':  6 obs. of  2 variables:
#>   ..$ A: num [1:6] 5 6 11 12 17 18
#>   ..$ B: num [1:6] 0 0 0 0 0 0

# Example 2: Create a nested list of `exdf` objects and then consolidate them
# into a regular list by combining the like-named elements. Here, some of the
# elements have columns not present in the others (for example,
# `nested_exdf_list$list_3$name_1`). However, these "extra" columns are removed
# before calling `rbind` and they do not appear in `consolidated_exdf_list`.
nested_exdf_list <- list(
  list_1 = list(
    name_1 = exdf(data.frame(A = c(1, 2), B = c(0, 0))),
    name_2 = exdf(data.frame(A = c(3, 4), B = c(0, 0))),
    name_3 = exdf(data.frame(A = c(5, 6), B = c(0, 0)))
  ),
  list_2 = list(
    name_1 = exdf(data.frame(A = c(7, 8), B = c(0, 0))),
    name_2 = exdf(data.frame(A = c(9, 10), B = c(0, 0))),
    name_3 = exdf(data.frame(A = c(11, 12), B = c(0, 0)))
  ),
  list_3 = list(
    name_1 = exdf(data.frame(A = c(13, 14), B = c(0, 0), C = c(-1, -2))),
    name_2 = exdf(data.frame(A = c(15, 16), B = c(0, 0), C = c(-1, -2))),
    name_3 = exdf(data.frame(A = c(17, 18), B = c(0, 0), C = c(-1, -2)))
  )
)

str(nested_exdf_list)
#> List of 3
#>  $ list_1:List of 3
#>   ..$ name_1:
#> Converting an `exdf` object to a `data.frame` before printing
#> 
#> 'data.frame':    2 obs. of  2 variables:
#>   .. ..$ A [NA] (NA): num [1:2] 1 2
#>   .. ..$ B [NA] (NA): num [1:2] 0 0
#>   ..$ name_2:
#> Converting an `exdf` object to a `data.frame` before printing
#> 
#> 'data.frame':    2 obs. of  2 variables:
#>   .. ..$ A [NA] (NA): num [1:2] 3 4
#>   .. ..$ B [NA] (NA): num [1:2] 0 0
#>   ..$ name_3:
#> Converting an `exdf` object to a `data.frame` before printing
#> 
#> 'data.frame':    2 obs. of  2 variables:
#>   .. ..$ A [NA] (NA): num [1:2] 5 6
#>   .. ..$ B [NA] (NA): num [1:2] 0 0
#>  $ list_2:List of 3
#>   ..$ name_1:
#> Converting an `exdf` object to a `data.frame` before printing
#> 
#> 'data.frame':    2 obs. of  2 variables:
#>   .. ..$ A [NA] (NA): num [1:2] 7 8
#>   .. ..$ B [NA] (NA): num [1:2] 0 0
#>   ..$ name_2:
#> Converting an `exdf` object to a `data.frame` before printing
#> 
#> 'data.frame':    2 obs. of  2 variables:
#>   .. ..$ A [NA] (NA): num [1:2] 9 10
#>   .. ..$ B [NA] (NA): num [1:2] 0 0
#>   ..$ name_3:
#> Converting an `exdf` object to a `data.frame` before printing
#> 
#> 'data.frame':    2 obs. of  2 variables:
#>   .. ..$ A [NA] (NA): num [1:2] 11 12
#>   .. ..$ B [NA] (NA): num [1:2] 0 0
#>  $ list_3:List of 3
#>   ..$ name_1:
#> Converting an `exdf` object to a `data.frame` before printing
#> 
#> 'data.frame':    2 obs. of  3 variables:
#>   .. ..$ A [NA] (NA): num [1:2] 13 14
#>   .. ..$ B [NA] (NA): num [1:2] 0 0
#>   .. ..$ C [NA] (NA): num [1:2] -1 -2
#>   ..$ name_2:
#> Converting an `exdf` object to a `data.frame` before printing
#> 
#> 'data.frame':    2 obs. of  3 variables:
#>   .. ..$ A [NA] (NA): num [1:2] 15 16
#>   .. ..$ B [NA] (NA): num [1:2] 0 0
#>   .. ..$ C [NA] (NA): num [1:2] -1 -2
#>   ..$ name_3:
#> Converting an `exdf` object to a `data.frame` before printing
#> 
#> 'data.frame':    2 obs. of  3 variables:
#>   .. ..$ A [NA] (NA): num [1:2] 17 18
#>   .. ..$ B [NA] (NA): num [1:2] 0 0
#>   .. ..$ C [NA] (NA): num [1:2] -1 -2

consolidated_exdf_list <- consolidate(nested_exdf_list)

str(consolidated_exdf_list)
#> List of 3
#>  $ name_1:
#> Converting an `exdf` object to a `data.frame` before printing
#> 
#> 'data.frame':    6 obs. of  2 variables:
#>   ..$ A [NA] (NA): num [1:6] 1 2 7 8 13 14
#>   ..$ B [NA] (NA): num [1:6] 0 0 0 0 0 0
#>  $ name_2:
#> Converting an `exdf` object to a `data.frame` before printing
#> 
#> 'data.frame':    6 obs. of  2 variables:
#>   ..$ A [NA] (NA): num [1:6] 3 4 9 10 15 16
#>   ..$ B [NA] (NA): num [1:6] 0 0 0 0 0 0
#>  $ name_3:
#> Converting an `exdf` object to a `data.frame` before printing
#> 
#> 'data.frame':    6 obs. of  2 variables:
#>   ..$ A [NA] (NA): num [1:6] 5 6 11 12 17 18
#>   ..$ B [NA] (NA): num [1:6] 0 0 0 0 0 0
```
