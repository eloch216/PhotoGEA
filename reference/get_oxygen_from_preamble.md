# Extract oxygen information from a Licor file

Ensures that the oxygen percentage (contained in the preamble of a Licor
file) is properly documented (including units).

## Usage

``` r
get_oxygen_from_preamble(licor_exdf)
```

## Arguments

- licor_exdf:

  An `exdf` object representing data from a photosynthetic gas exchange
  measurement system. It should contain a column called `Oxygen`; this
  should automatically be the case if `licor_exdf` was created by
  [`read_gasex_file`](https://eloch216.github.io/PhotoGEA/reference/read_gasex_file.md).

## Details

Licor LI-6800 log files include the oxygen concentration as an entry in
the preamble, which is automatically converted to a column by
[`read_licor_6800_Excel`](https://eloch216.github.io/PhotoGEA/reference/read_licor_6800_Excel.md)
and
[`read_licor_6800_plaintext`](https://eloch216.github.io/PhotoGEA/reference/read_licor_6800_plaintext.md).
However, the Licor log files do not specify units for this value; they
should be `percent`.

This function adds the correct units to the `Oxygen` column, or creates
one (with values initialized to `NA`) if the `Oxygen` column is not
present.

Also note that this function is called automatically by
[`read_licor_6800_Excel`](https://eloch216.github.io/PhotoGEA/reference/read_licor_6800_Excel.md)
and
[`read_licor_6800_plaintext`](https://eloch216.github.io/PhotoGEA/reference/read_licor_6800_plaintext.md)
whenever `get_oxygen` is `TRUE`; this happens by default, so it is rare
for a user to call this function directly.

## Value

An `exdf` object based on `licor_exdf` that includes the `Oxygen` column
units.

## Examples

``` r
# Example: Read data from a Licor log file and get the oxygen information from
# the preamble

# Read the file without automatically specifying oxygen units
licor_data <- read_gasex_file(
  PhotoGEA_example_file_path('licor_for_gm_site11.xlsx'),
  get_oxygen = FALSE
)

# Here we can see that the oxygen column is included, but without units
print(licor_data$categories$Oxygen)
#> [1] "SysConst"
print(licor_data$units$Oxygen)
#> [1] "NA"
print(licor_data[, 'Oxygen'])
#>  [1] 21 21 21 21 21 21 21 21 21 21 21 21

# Include the oxygen units
licor_data <- get_oxygen_from_preamble(licor_data)

# The units have changed, but the category and values have not
print(licor_data$categories$Oxygen)
#> [1] "SysConst"
print(licor_data$units$Oxygen)
#> [1] "percent"
print(licor_data[, 'Oxygen'])
#>  [1] 21 21 21 21 21 21 21 21 21 21 21 21
```
