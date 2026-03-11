# Set values, units, and categories for a column in a table

Sets the value, units, and/or category of a new or existing column of a
table-like object.

## Usage

``` r
set_variable(
    data_table,
    name,
    units = NULL,
    category = NULL,
    value = NA,
    id_column = NULL,
    value_table = NULL
  )
```

## Arguments

- data_table:

  A table-like R object such as a data frame or an `exdf`.

- name:

  The name of the column to be added to `data_table`.

- units:

  The units of the column to be added to `data_table`.

- category:

  The category of the column to be added to `data_table`.

- value:

  The value of the column to be added to `data_table`.

- id_column:

  The name of an identifier column in `data_table`.

- value_table:

  A list of named elements, where the name of each element is a possible
  value of the `id_column` and the value of each element is the
  corresponding value that the `name` column should take.

## Details

There are two main "modes" for setting the value of the new column: it
can be set to a fixed value (using the `value` input argument), or it
can be set according to the values of another column (using the
`id_column` and `value_table` input arguments). The latter method is
useful when different values must be specified for different treatments
within the data set.

In greater detail, this function attempts to set the value of a new or
existing column in an `exdf` object according to the following rules:

- The value of the `name` column of `data_table` will be set to `value`;
  this assignment follows the usual rules; in other words, `value` could
  be a single value or a vector of length `nrow(data_table)`.

- If `units` and `categories` are both `NULL`, the units and category
  will not be specified. In this case, if the `name` column already
  exists, its units and category will remain the same; if the `name`
  column is new, it will be initialized with `NA` for its units and
  category.

- If either `units` \_or\_ `category` is not `NULL`, the units and
  category for the `name` column \_will\_ be specified. In this case, if
  one of `units` or `category` \_is\_ `NULL`, its value will be set to
  `NA`.

- If `id_column` is not `NULL`, then the `value_table` will be used to
  set different values of the `name` column for each specified value of
  `id_column`. For example, if `id_column` is `species` and
  `value_table = list(soybean = 1, tobacco = 2)`, then the `name` column
  will be set to `1` when `species` is `'soybean'` and `2` when
  `species` is `'tobacco'`. For any other values of species (such as
  `'maize'`), the value of `name` will still be `value`. \*\*Note\*\*:
  values of the `id_column` will be converted using `as.character`
  before making comparisons.

For other table-like objects, such as data frames, only the values will
be set, and the units and categories will be ignored.

## Value

An object based on `data_table` with new and/or modified columns.

## See also

[`exdf`](https://eloch216.github.io/PhotoGEA/reference/exdf.md)

## Examples

``` r
# Create a simple exdf object with two columns (`A` and `B`) and default values
# for its units and categories.
simple_exdf <- exdf(data.frame(A = c(3, 2, 7, 9), B = c(4, 5, 1, 8)))

print(simple_exdf)
#> 
#> Converting an `exdf` object to a `data.frame` before printing
#> 
#>   A [NA] (NA) B [NA] (NA)
#> 1           3           4
#> 2           2           5
#> 3           7           1
#> 4           9           8

# Add a new column called 'C' with units 'u1' and category 'cat1' whose value is
# 1000.
simple_exdf <- set_variable(simple_exdf, 'C', 'u1', 'cat1', 1000)

# Set the value of the 'B' column to 2000 when 'A' is 3, to 3000 when 'A' is 9,
# and to 4000 for all other values of 'A'. Do not modify its units or category.
simple_exdf <- set_variable(
  simple_exdf,
  'B',
  value = 4000,
  id_column = 'A',
  value_table = list('3' = 2000, '9' = 3000)
)

print(simple_exdf)
#> 
#> Converting an `exdf` object to a `data.frame` before printing
#> 
#>   A [NA] (NA) B [NA] (NA) C [cat1] (u1)
#> 1           3        2000          1000
#> 2           2        4000          1000
#> 3           7        4000          1000
#> 4           9        3000          1000

# Take the same operations, but using a data frame instead
simple_df <- data.frame(A = c(3, 2, 7, 9), B = c(4, 5, 1, 8))

simple_df <- set_variable(simple_exdf$main_data, 'C', 'u1', 'cat1', 1000)

simple_df <- set_variable(
  simple_df,
  'B',
  value = 4000,
  id_column = 'A',
  value_table = list('3' = 2000, '9' = 3000)
)

print(simple_df)
#>   A    B    C
#> 1 3 2000 1000
#> 2 2 4000 1000
#> 3 7 4000 1000
#> 4 9 3000 1000

# As a more realistic example, load a Licor file and set different values of
# mesophyll conductance for each species in the data set.
licor_file <- read_gasex_file(
  PhotoGEA_example_file_path('ball_berry_1.xlsx')
)

licor_file <- set_variable(
  licor_file,
  'gmc',
  'mol m^(-2) s^(-1) bar^(-1)',
  '',
  id_column = 'species',
  value_table = list(soybean = 0.9, tobacco = 1.1)
)

print(licor_file[, c('species', 'gmc'), TRUE])
#> 
#> Converting an `exdf` object to a `data.frame` before printing
#> 
#>    species [UserDefCon] (NA) gmc [] (mol m^(-2) s^(-1) bar^(-1))
#> 1                    soybean                                 0.9
#> 2                    soybean                                 0.9
#> 3                    soybean                                 0.9
#> 4                    soybean                                 0.9
#> 5                    soybean                                 0.9
#> 6                    soybean                                 0.9
#> 7                    soybean                                 0.9
#> 8                    soybean                                 0.9
#> 9                    soybean                                 0.9
#> 10                   soybean                                 0.9
#> 11                   soybean                                 0.9
#> 12                   soybean                                 0.9
#> 13                   soybean                                 0.9
#> 14                   soybean                                 0.9
#> 15                   soybean                                 0.9
#> 16                   soybean                                 0.9
#> 17                   soybean                                 0.9
#> 18                   soybean                                 0.9
#> 19                   soybean                                 0.9
#> 20                   soybean                                 0.9
#> 21                   soybean                                 0.9
#> 22                   tobacco                                 1.1
#> 23                   tobacco                                 1.1
#> 24                   tobacco                                 1.1
#> 25                   tobacco                                 1.1
#> 26                   tobacco                                 1.1
#> 27                   tobacco                                 1.1
#> 28                   tobacco                                 1.1
```
