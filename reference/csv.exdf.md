# Read and write CSV files representing an exdf object

Functions for reading and writing CSV files that represent an
[`exdf`](https://eloch216.github.io/PhotoGEA/reference/exdf.md) object.

## Usage

``` r
read.csv.exdf(file, ...)

  write.csv.exdf(x, file, ...)
```

## Arguments

- file:

  The name of the file which the data are to be read from; to be passed
  to [`read.csv`](https://rdrr.io/r/utils/read.table.html).

- ...:

  Additional arguments to be passed to
  [`read.csv`](https://rdrr.io/r/utils/read.table.html) or
  [`write.csv`](https://rdrr.io/r/utils/write.table.html). Note that
  some arguments cannot be specified; an error message will be sent if a
  used attempts to set one of these forbidden arguments.

- x:

  An `exdf` object to be written to a CSV file.

## Details

An `exdf` object can be written to a CSV file by directly calling
[`write.csv`](https://rdrr.io/r/utils/write.table.html), but this
approach causes some column names to be unintentionally modified. For
example, any spaces will be replaced by periods. This can potentially
cause problems when reloading the data from the CSV file.

Instead, it is preferred to use `write.csv.exdf`, which will not modify
any column names. When writing the CSV file, it is saved with the column
names in the first row, the categories in the second row, the units in
the third row, and the data in the remaining rows.

The resulting file can be read using `read.csv.exdf`. Here, the names,
categories, and units are read from the first three rows of the
specified file, and the data values from the remaining rows. An
[`exdf`](https://eloch216.github.io/PhotoGEA/reference/exdf.md) object
is then created from this information.

## Value

The `write.csv.exdf` function does not return anything. The
`read.csv.exdf` function returns an
[`exdf`](https://eloch216.github.io/PhotoGEA/reference/exdf.md) object
representing the contents of `file`.

## Examples

``` r
# Read a CSV file included with the PhotoGEA package; this file was created
# using `write.csv.exdf`.
licor_file <- read.csv.exdf(
  PhotoGEA_example_file_path('ball_berry_1.csv')
)

# Now rewrite this to a temporary CSV file
tf <- tempfile(fileext = ".csv")
tf
#> [1] "/tmp/RtmpRcVtEp/file1a4051d46371.csv"

write.csv.exdf(licor_file, tf)
```
