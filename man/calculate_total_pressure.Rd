\name{calculate_total_pressure}

\alias{calculate_total_pressure}

\title{Calculate the total pressure in bar}

\description{
  Calculates the total pressure in \code{bar}. Licor gas exchange measurement
  systems report both the abient air pressure (\code{Pa}) and the chamber
  overpressure (\code{DeltaPcham}) in \code{kPa}; the total pressure in the
  chamber is therefore given by the sum of these two columns. This function can
  accomodate alternative column names for the variables taken from Licor log
  files in case they change at some point in the future. This function also
  checks the units of each required column and will produce an error if any
  units are incorrect.
}

\usage{
  calculate_total_pressure(
    exdf_obj,
    pa_column_name = 'Pa',
    deltapcham_column_name = 'DeltaPcham'
  )
}

\arguments{
  \item{exdf_obj}{An \code{exdf} object that contains pressure measurements.}

  \item{pa_column_name}{
    The name of the column in \code{exdf_obj} that contains the ambient air
    pressure in \code{kPa}.
  }

  \item{deltapcham_column_name}{
    The name of the column in \code{exdf_obj} that contains the chamber
    overpressure in \code{kPa}.
  }
}

\details{
  If \code{deltapcham_column_name} is \code{NA}, this function will simply
  convert the values in the \code{pa_column_name} to units of \code{bar}.
  Otherwise, the values from the \code{pa_column_name} and
  \code{deltapcham_column_name} columns will be added together and converted to
  \code{bar}.
}

\value{
  An \code{exdf} object based on \code{exdf_obj} that includes the total
  pressure values in a new column called \code{total_pressure}. The category of
  this new column is \code{calculate_total_pressure} to indicate that it was
  created using this function.
}

\examples{
# Read an example Licor file included in the PhotoGEA package and calculate the
# total pressure.
licor_file <- read_gasex_file(
  system.file('extdata', 'ball_berry_1.xlsx', package = 'PhotoGEA', mustWork = TRUE)
)

licor_file <- calculate_total_pressure(licor_file)

licor_file$units$total_pressure      # View the units of the new `total_pressure` column
licor_file$categories$total_pressure # View the category of the new `total_pressure` column
licor_file[, 'total_pressure']       # View the values of the new `total_pressure` column
}

\concept{exdf}
