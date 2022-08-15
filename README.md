## PhotoGEA

**PhotoGEA** (short for **photo**synthetic **g**as **e**xchange **a**nalysis) is
an R package that provides a suite of tools for loading, processing, and
analyzing photosynthetic gas exchange data.

This package is a work in progress and is not yet fully documented. As
documentation proceeds, functions are subject to modification or removal without
notice.

### Installing the R Package

The easiest way to install `PhotoGEA` is to type the following from within the R
terminal:

```
remotes::install_github('eloch216/PhotoGEA')
```

Note that this method requires the `remotes` package, which can be installed
from within R by typing `install.packages('remotes')`.

### Learning to Use PhotoGEA

The best way to learn about using `PhotoGEA` is to visit the
[PhotoGEA website](https://eloch216.github.io/PhotoGEA/index.html), which
includes articles that describe several important use cases. Currently, the
following articles are available:

- [Analyzing Ball-Berry Data](https://eloch216.github.io/PhotoGEA/articles/analyzing_ball_berry_data.html)
- [Analyzing TDL Data](https://eloch216.github.io/PhotoGEA/articles/analyzing_tdl_data.html)

### Example Scripts

Several example scripts are provided in the `example_scripts` directory of the
[source code repository](https://github.com/eloch216/PhotoGEA). To run one of
these scripts, set the working directory to a folder that contains a local copy
of the script and use the `source` command to execute the code in the script. No
guarantees are made that these scripts will run on your machine or be
compatible with your data, but they may be a useful source of ideas.

### License

The `PhotoGEA` R package, its documentation, and its associated example scripts
are licensed under the MIT license.
