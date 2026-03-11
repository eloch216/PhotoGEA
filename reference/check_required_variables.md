# Make sure required variables exist

Checks whether the input table has the required variables.

## Usage

``` r
check_required_variables(x, required_variables, check_NA = TRUE)

  # S3 method for class 'data.frame'
check_required_variables(x, required_variables, check_NA = TRUE)

  # S3 method for class 'exdf'
check_required_variables(x, required_variables, check_NA = TRUE)
```

## Arguments

- x:

  A table-like R object such as a data frame or an `exdf`.

- required_variables:

  A set of variables that must each be included in `x` as columns.

- check_NA:

  A logical value indicating whether to check for columns that are all
  `NA`; see below.

## Details

`check_required_variables` is generic, with methods defined for data
frames and `exdf` objects.

When `x` is an `exdf`, the `required_variables` input argument must be a
list of named strings, where the name of each element specifies the name
of a column that must be included in `x`, while the value of each column
specifies the corresponding units for that column. If the value is `NA`,
no unit checking will be performed.

When `x` is a `data.frame`, the `required_variables` input argument can
be specified as a list (as if `x` were an `exdf` object) or as a
character vector specifying the names of columns that should be included
in `x`.

The required variables will be checked as follows:

- If any required variable columns are missing from the table, an
  informative error message will be thrown.

- If `check_NA` is `TRUE` and any required variable columns are entirely
  `NA`, an informative error message will be thrown.

- If any required variable colums have incorrect units, an informative
  error message will be thrown. (Only applies to `exdf` objects.)

Otherwise, `check_required_variables` will have no output and produce no
messages.

This function is used internally by many other functions from the
`PhotoGEA` package to check for important columns and make sure they
have the correct units. For example, see the code for
[`apply_gm`](https://eloch216.github.io/PhotoGEA/reference/apply_gm.md)
by typing
[`PhotoGEA::apply_gm`](https://eloch216.github.io/PhotoGEA/reference/apply_gm.md)
in the R terminal.

## Value

The `check_required_variables` function does not return anything.

## See also

[`exdf`](https://eloch216.github.io/PhotoGEA/reference/exdf.md)

## Examples

``` r
# Create a simple exdf object
simple_exdf <- exdf(
  data.frame(A = c(3, 2, 7, 9), B = c(4, 5, 1, 8)),
  data.frame(A = 'm', B = 's', stringsAsFactors = FALSE),
  data.frame(A = 'Cat1', B = 'Cat2', stringsAsFactors = FALSE)
)

# Confirm that columns named `A` and `B` are in the object, and that they have
# units of `m` and `s`, respectively.
check_required_variables(simple_exdf, list(A = 'm', B = 's'))

# Confirm that columns named `A` and `B` are in the object, but only check units
# for the `A` column.
check_required_variables(simple_exdf, list(A = 'm', B = NA))

# Use the data frame method on `simple_exdf$main_data` to confirm that columns
# named `A` and `B` are present
check_required_variables(simple_exdf$main_data, c('A', 'B'))
```
