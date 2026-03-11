# Extended data frame

An "extended data frame" (`exdf`) is an object similar to a data frame,
but which also contains information about the units and categories of
each column.

## Usage

``` r
exdf(
    main_data = data.frame(),
    units = NULL,
    categories = NULL,
    ...
  )
```

## Arguments

- main_data:

  A data frame.

- units:

  A data frame with the same columns as `main_data` (or a subset of the
  columns in `main_data`) but with just one row, where each entry
  describes the units for the corresponding column of `main_data`. If
  `units` is `NULL`, it will be initialized with `NA` for each column.
  The units of any columns in `main_data` that are not present in
  `units` will also be initialized to `NA`.

- categories:

  A data frame with the same columns as `main_data` (or a subset of the
  columns in `main_data`) but with just one row, where each entry
  describes the category for the corresponding column of `main_data`. If
  `categories` is `NULL`, it will be initialized with `NA` for each
  column. The categories of any columns in `main_data` that are not
  present in `catgories` will also be initialized to `NA`.

- ...:

  Any additional properties to include as entries in the resulting
  `exdf` object; these must be passed as named arguments.

## Details

The `exdf` class was originally created as a way to represent the
contents of a Licor Excel file in an R structure. In Licor Excel files,
each column has a name, units, and a category; for example, the column
for values of net assimilation rate is called `A`, has units of
`micromol / m^2 / s`, and is categorized as a `GasEx` variable.

From a technical point of view, an `exdf` object is simply a list with
three required elements: `main_data`, `units`, and `categories`. Each of
these should be a data frame with the same column names, as described
above. It is also possible for an `exdf` object to have additional
entries such as a `filename` that stores the name of the file that was
used to create the `exdf`.

Several S3 methods have been defined for `exdf` objects, following the
general guidance from [Advanced R on S3
classes](http://adv-r.had.co.nz/S3.md):

- [`is.exdf`](https://eloch216.github.io/PhotoGEA/reference/is.exdf.md)

- [`as.data.frame.exdf`](https://eloch216.github.io/PhotoGEA/reference/as.data.frame.exdf.md)

- [`print.exdf`](https://eloch216.github.io/PhotoGEA/reference/print.exdf.md)

- [`str.exdf`](https://eloch216.github.io/PhotoGEA/reference/str.exdf.md)

- [`length.exdf`](https://eloch216.github.io/PhotoGEA/reference/length.exdf.md)

- [`dim.exdf`](https://eloch216.github.io/PhotoGEA/reference/dim.exdf.md)

- [`dimnames.exdf`](https://eloch216.github.io/PhotoGEA/reference/dimnames.exdf.md)

- [`[.exdf`](https://eloch216.github.io/PhotoGEA/reference/extract.exdf.md)

- `[<-.exdf`

- [`rbind.exdf`](https://eloch216.github.io/PhotoGEA/reference/cbind.exdf.md)

- [`cbind.exdf`](https://eloch216.github.io/PhotoGEA/reference/cbind.exdf.md)

- [`split.exdf`](https://eloch216.github.io/PhotoGEA/reference/split.exdf.md)

- [`by.exdf`](https://eloch216.github.io/PhotoGEA/reference/by.exdf.md)

Note that the column names of `main_data`, `units`, and `categories`
must be unique; the
[`make.unique`](https://rdrr.io/r/base/make.unique.html) function can be
useful for ensuring this.

## Value

An `exdf` object as described above.

## Examples

``` r
# Example 1: Creating a simple exdf object with two columns (`A` and `B`) and
# default values for its units and categories. There are four values of each
# variable.
exdf(data.frame(A = c(3, 2, 7, 9), B = c(4, 5, 1, 8)))
#> 
#> Converting an `exdf` object to a `data.frame` before printing
#> 
#>   A [NA] (NA) B [NA] (NA)
#> 1           3           4
#> 2           2           5
#> 3           7           1
#> 4           9           8

# Example 2: Creating a simple exdf object with two columns (`A` and `B`) that
# have units of `m` and `s`, respectively, and categories of `Cat1` and `Cat2`,
# respectively. There are four values of each variable.
exdf(
  data.frame(A = c(3, 2, 7, 9), B = c(4, 5, 1, 8)),
  data.frame(A = 'm', B = 's'),
  data.frame(A = 'Cat1', B = 'Cat2')
)
#> 
#> Converting an `exdf` object to a `data.frame` before printing
#> 
#>   A [Cat1] (m) B [Cat2] (s)
#> 1            3            4
#> 2            2            5
#> 3            7            1
#> 4            9            8
```
