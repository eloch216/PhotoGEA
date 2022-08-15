\name{read_tdl_file}

\alias{read_tdl_file}

\title{Reading a TDL data file}

\description{
  Tool for reading tunable diode laser (TDL) output files created by Campbell
  Scientific CR3000 data loggers and storing their contents in
  \code{\link{exdf}} objects.
}

\usage{
  read_tdl_file(
    file_name,
    rows_to_skip,
    variable_name_row,
    variable_unit_row,
    data_start_row,
    timestamp_colname
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

  \item{timestamp_colname}{
    The name of the column that contains the timestamp of each measurement;
    typically, this is \code{'TIMESTAMP'}.
  }
}

\details{
  When reading a TDL data file, it is useful to identify the timestamp column
  so its values can be properly interpreted as \code{\link{POSIXlt}} objects.
}

\value{
  An \code{exdf} object that fully includes all the data from the TDL data file.
  In addition to the required elements of an \code{exdf} object, the following
  "extra" elements are also included:
  \itemize{
    \item \code{file_name}: A copy of the input argument with the same name.
    \item \code{rows_to_skip}: A copy of the input argument with the same name
    \item \code{variable_name_row}: A copy of the input argument with the same
          name.
    \item \code{variable_unit_row}: A copy of the input argument with the same
          name.
    \item \code{data_start_row}: A copy of the input argument with the same
          name.
    \item \code{timestamp_colname}: A copy of the input argument with the same
          name.
  }
}

\examples{
# Example: reading a TDL file that is included with the PhotoGEA package. Note
# that although named input arguments are not required, here we supply them for
# clarity.
tdl_file <- read_tdl_file(
  system.file("extdata", "tdl_sampling_1.dat", package = "PhotoGEA"),
  rows_to_skip = 1,
  variable_name_row = 2,
  variable_unit_row = 3,
  data_start_row = 5,
  timestamp_colname = 'TIMESTAMP'
)

tdl_file$file_name # A record of where the data came from
str(tdl_file)      # View the contents of the exdf object's main_data
}

\concept{exdf}
