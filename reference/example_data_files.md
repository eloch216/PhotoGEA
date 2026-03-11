# Example data files

The PhotoGEA package includes several data files that can be used to
demonstrate different functions and analysis techniques.

## Details

The following files are included with the package:

- `ball_berry_1.xlsx` and `ball_berry_2.xlsx`: Two log files created by
  Licor Li-6800 portable gas exchange measurement systems. These log
  files each contain several Ball-Berry curves. Several user constants
  were defined in these logs that can be used to identify individual
  curves or subsets of curves: `species`, `plot`, and `instrument`.
  These files are used in the "Analyzing Ball-Berry Data" vignette and
  in other examples.

- `ball_berry_1.csv`: A CSV version of `ball_berry_1.xlsx`, which was
  created by reading the Excel file with
  [`read_gasex_file`](https://eloch216.github.io/PhotoGEA/reference/read_gasex_file.md)
  and then saving it using
  [`write.csv.exdf`](https://eloch216.github.io/PhotoGEA/reference/csv.exdf.md).
  This can be done as follows:
  `tmp <- read_gasex_file(PhotoGEA_example_file_path('ball_berry_1.xlsx')); write.csv.exdf(tmp, 'ball_berry_1.csv')`

- `c3_aci_1.xlsx` and `c3_aci_2.xlsx`: Two log files created by Licor
  Li-6800 portable gas exchange measurement systems. These log files
  each contain several C3 CO2 response (or A-Ci) curves. Several user
  constants were defined in these logs that can be used to identify
  individual curves or subsets of curves: `species`, `plot`, and
  `instrument`. These files are used in the "Analyzing C3 A-Ci Curves"
  vignette and in other examples. The `Remarks` sheet of `c3_aci_2.xlsx`
  was deleted from the original version, and an updated value of Oxygen
  was added after observation 10, as a test for
  [`read_licor_6800_Excel`](https://eloch216.github.io/PhotoGEA/reference/read_licor_6800_Excel.md).

- `c4_aci_1.xlsx` and `c4_aci_2.xlsx`: Two log files created by Licor
  Li-6800 portable gas exchange measurement systems. These log files
  each contain several C4 CO2 response (or A-Ci) curves. Several user
  constants were defined in these logs that can be used to identify
  individual curves or subsets of curves: `species`, `plot`, and
  `instrument`. These files are used in the "Analyzing C4 A-Ci Curves"
  vignette and in other examples.

- `tdl_sampling_1.dat` and `tdl_sampling_2.dat`: Two log files created
  by a Campbell Scientific CR3000 data logger, representing data from a
  tunable diode laser (TDL) system. These files are used in the
  "Analyzing TDL Data" vignette and in other examples.

- `plaintext_licor_file`: A log file created by a Licor Li-6800 portable
  gas exchange measurement system. This file contains several CO2
  response (or A-Ci) curves. Several user constants were defined in this
  log that can be used to identify individual curves or subsets of
  curves: `species`, `plot`, and `instrument`. Several updated values of
  Oxygen were added as a test of
  [`read_licor_6800_plaintext`](https://eloch216.github.io/PhotoGEA/reference/read_licor_6800_plaintext.md).

- `plaintext_licor_file_v2`: A log file based on `plaintext_licor_file`
  that has two separate `[Data]` and `[Header]` sections, as if the log
  file had been closed and reopened halfway through the measurement
  sequence. It also has an extra blank line at the end. Several updated
  values of Oxygen were added as a test of
  [`read_licor_6800_plaintext`](https://eloch216.github.io/PhotoGEA/reference/read_licor_6800_plaintext.md).

- `licor_for_gm_site11.xlsx`, `licor_for_gm_site13.xslsx`, and
  `tdl_for_gm`: Two Licor Li-6800 log files and a CR3000 TDL log file,
  respectively. These files are used as an example of loading and
  processing combined gas exchange and isotope discrimination
  measurements. Each Licor log file includes 6 points measured with the
  `CO2_r` setpoint set to 715 ppm and 6 points with the setpoint set to
  450 ppm.

Since none of these data files have been published, noise has been added
to the original data. Thus, they are similar to real measurements, but
no useful conclusions can be drawn from them.

After installing \`PhotoGEA\`, copies of these files will be stored in
the R package directory (in the `PhotoGEA/extdata` subdirectory). This
location will be unique to your computer, but full paths to these files
can be obtained using the
[`PhotoGEA_example_file_path`](https://eloch216.github.io/PhotoGEA/reference/PhotoGEA_example_file_path.md)
function.

## Examples

``` r
# Print full paths to the example files
PhotoGEA_example_file_path('ball_berry_1.xlsx')
#> [1] "/home/runner/work/_temp/Library/PhotoGEA/extdata/ball_berry_1.xlsx"
PhotoGEA_example_file_path('ball_berry_2.xlsx')
#> [1] "/home/runner/work/_temp/Library/PhotoGEA/extdata/ball_berry_2.xlsx"
PhotoGEA_example_file_path('c3_aci_1.xlsx')
#> [1] "/home/runner/work/_temp/Library/PhotoGEA/extdata/c3_aci_1.xlsx"
PhotoGEA_example_file_path('c3_aci_2.xlsx')
#> [1] "/home/runner/work/_temp/Library/PhotoGEA/extdata/c3_aci_2.xlsx"
PhotoGEA_example_file_path('c4_aci_1.xlsx')
#> [1] "/home/runner/work/_temp/Library/PhotoGEA/extdata/c4_aci_1.xlsx"
PhotoGEA_example_file_path('c4_aci_2.xlsx')
#> [1] "/home/runner/work/_temp/Library/PhotoGEA/extdata/c4_aci_2.xlsx"
PhotoGEA_example_file_path('licor_for_gm_site11.xlsx')
#> [1] "/home/runner/work/_temp/Library/PhotoGEA/extdata/licor_for_gm_site11.xlsx"
PhotoGEA_example_file_path('licor_for_gm_site13.xlsx')
#> [1] "/home/runner/work/_temp/Library/PhotoGEA/extdata/licor_for_gm_site13.xlsx"
PhotoGEA_example_file_path('plaintext_licor_file')
#> [1] "/home/runner/work/_temp/Library/PhotoGEA/extdata/plaintext_licor_file"
PhotoGEA_example_file_path('plaintext_licor_file_v2')
#> [1] "/home/runner/work/_temp/Library/PhotoGEA/extdata/plaintext_licor_file_v2"
PhotoGEA_example_file_path('tdl_for_gm.dat')
#> [1] "/home/runner/work/_temp/Library/PhotoGEA/extdata/tdl_for_gm.dat"
PhotoGEA_example_file_path('tdl_sampling_1.dat')
#> [1] "/home/runner/work/_temp/Library/PhotoGEA/extdata/tdl_sampling_1.dat"
PhotoGEA_example_file_path('tdl_sampling_2.dat')
#> [1] "/home/runner/work/_temp/Library/PhotoGEA/extdata/tdl_sampling_2.dat"
```
