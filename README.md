## licor-processing-and-analysis

This repository contains the source code for an R package called **PhotoGEA**
(short for **photo**synthetic **g**as **e**xchange **a**nalysis), which is a
suite of tools for loading, processing, and analyzing photosynthetic gas
exchange data.

This package is a work in progress and is not yet fully documented.

### R package installation steps
Obtain the (unzipped) source code from GitHub and install within R using the
following set of commands:

```
setwd('path_to_unzipped_directory')
install.packages('PhotoGEA', repos=NULL, type='SOURCE')
```

### Using example scripts
Several example scripts are provided in the `example_scripts` directory:
- `extract_licor_data_for_gm.R`
- `plot_gm_analysis.R`
- `plot_response_curve_analysis.R`

To run one of these scripts, set the working directory to `example_scripts` and
use the `source` command to execute the code in the script.
