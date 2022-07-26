## licor-processing-and-analysis

This repository contains the source code for an R package called **PhotoGEA**
(short for **photo**synthetic **g**as **e**xchange **a**nalysis), which is a
suite of tools for loading, processing, and analyzing photosynthetic gas
exchange data.

This package is a work in progress and is not yet fully documented.

### R package installation steps

Obtain the (unzipped) source code from GitHub and install from the command line
using the following set of commands, which assume the source code is contained
in a directory called `licor-processing-and-analysis`, which is a subdirectory
of `path_to_unzipped_directory`:

```
cd path_to_unzipped_directory
R CMD build licor-processing-and-analysis
R CMD INSTALL output_file
```

where `output_file` is the name of the file produced by `R CMD build`; its name
will be something like `PhotoGEA_0.2.0.tar.gz`, although the details will depend
on your operating system. This is the preferred installation method because it
will build and install the package vignettes.

Alternatively, PhotoGEA can be installed from within R using these commands:

```
setwd('path_to_unzipped_directory')
install.packages('licor-processing-and-analysis', repos=NULL, type='SOURCE')
```

However, the vignettes will not automatically be built using this method.

### Viewing the vignettes

To see the available vignettes, type `browseVignettes("PhotoGEA")` from within
R. (This command will only work if you follow the preferred installation method
described above.) The vignettes describe important use cases for the PhotoGEA
package and can be used as the basis for your own scripts.

### Using example scripts

Additionally, several example scripts are provided in the `example_scripts`
directory. To run one of these scripts, set the working directory to
`example_scripts` and use the `source` command to execute the code in the
script.

### License

The PhotoGEA R package and associated example scripts are licensed under the MIT
license.
