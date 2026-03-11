# Extract TDL valve information from file name

Determines the TDL valve number from a photosynthetic gas exchange
system log file name.

## Usage

``` r
get_sample_valve_from_filename(
    exdf_obj,
    reference_table = NULL
  )
```

## Arguments

- exdf_obj:

  An `exdf` object representing data from a photosynthetic gas exchange
  measurement system. The `exdf_obj$file_name` field must be defined and
  contain the file name; this will automatically be the case if
  `exdf_obj` was created by
  [`read_gasex_file`](https://eloch216.github.io/PhotoGEA/reference/read_gasex_file.md).

- reference_table:

  An optional list of named elements, where the name of each element is
  a Licor sample line valve number (as a character) and the value of
  each element is the corresponding Licor reference line valve number.

## Details

When making combined gas exchange and isotope discrimination
measurements using a portable photosynthetic gas exchange system (such
as a Licor LI-6800) coupled with a tunable diode laser (TDL) absorption
spectroscopy system, the TDL's gas handling system cycles through
several gas lines (or sites) by opening and closing valves. When
analyzing such data, a key step is to identify which TDL valve numbers
correspond to the sample and reference gas lines of the Licor.

At UIUC, there is a convention for designating the sample line valve
numbers in the Licor file names, where `"siteNN"` or `"site NN"` means
that the Licor's sample line is valve `NN` in the TDL data file. The
`get_sample_valve_from_filename` function extracts the valve number from
the file name and stores it in a new column in `exdf_obj` called
`valve_number_s`.

Optionally, it is also possible to specify the reference line valve
number corresponding to each sample line valve number using the
`reference_table` input argument. Reference line valve numbers will be
stored in the `valve_number_r` column.

## Value

An `exdf` object based on `exdf_obj` that includes the Licor sample line
valve number as a new column called `valve_number_s` and (optionally)
the Licor reference line valve number as a new column called
`valve_number_r`.

## Examples

``` r
## In this example we load a gas exchange data file and determine the TDL valve
## numbers from its file name

# Read the gas exchange data
licor_data <- read_gasex_file(
  PhotoGEA_example_file_path('licor_for_gm_site11.xlsx'),
)

# Get TDL valve information from Licor file name; for this TDL system, the
# reference valve is 12 when the sample valve is 11
licor_data <- get_sample_valve_from_filename(licor_data, list('11' = 12))

# View the results
licor_data[, c('obs', 'valve_number_s', 'valve_number_r')]
#>    obs valve_number_s valve_number_r
#> 1    1             11             12
#> 2    2             11             12
#> 3    3             11             12
#> 4    4             11             12
#> 5    5             11             12
#> 6    6             11             12
#> 7    7             11             12
#> 8    8             11             12
#> 9    9             11             12
#> 10  10             11             12
#> 11  11             11             12
#> 12  12             11             12
```
