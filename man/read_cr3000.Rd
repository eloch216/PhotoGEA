\name{read_cr3000}

\alias{read_cr3000}

\title{Reading a CR3000 data file}

\description{
  Tool for reading output files created by Campbell Scientific CR3000 data
  loggers and storing their contents in \code{\link{exdf}} objects.
}

\usage{
  read_cr3000(
    file_name,
    rows_to_skip = 1,
    variable_name_row = 2,
    variable_unit_row = 3,
    data_start_row = 5,
    \dots
  )
}

\arguments{
  \item{file_name}{
    A relative or absolute path to a \code{.dat} file containing TDL data.
  }

  \item{rows_to_skip}{
    The number of rows to skip at the beginning of the file; the first row in a
    TDL file typically has fewer columns than the others, which causes problems
    when storing it as a table.
  }

  \item{variable_name_row}{
    The row number in the TDL file containing the names of the variables
    (\code{RECORD}, \code{Conc12C_Avg}, etc).
  }

  \item{variable_unit_row}{
    The row number in the TDL file containing the units of the variables
    (\code{ppm}, \code{V}, etc).
  }

  \item{data_start_row}{
    The first row number of the table containing the measured data.
  }

  \item{\dots}{
    Additional arguments to be passed to \code{\link[utils]{read.csv}}.
  }
}

\value{
  An \code{exdf} object that fully includes all the data from the CR3000 output
  file. In addition to the elements described in the documentation for
  \code{\link{read_gasex_file}}, the following "extra" elements are also
  included:
  \itemize{
    \item \code{rows_to_skip}: A copy of the input argument with the same name
    \item \code{variable_name_row}: A copy of the input argument with the same
          name.
    \item \code{variable_unit_row}: A copy of the input argument with the same
          name.
    \item \code{data_start_row}: A copy of the input argument with the same
          name.
  }
}

\examples{
# Example: reading a TDL file that is included with the PhotoGEA package.
tdl_file <- PhotoGEA:::read_cr3000(
  PhotoGEA_example_file_path('tdl_sampling_1.dat')
)

tdl_file$file_name # A record of where the data came from
str(tdl_file)      # View the contents of the exdf object's main_data
}

\seealso{\code{\link{read_gasex_file}}}

\concept{exdf}
