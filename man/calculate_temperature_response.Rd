\name{calculate_temperature_response}

\alias{calculate_temperature_response}

\title{Calculate temperature-dependent parameter values}

\description{
  Calculate leaf-temperature-dependent values of various parameters using
  various temperature response functions.
}

\usage{
  calculate_temperature_response(
    exdf_obj,
    temperature_response_parameters,
    tleaf_column_name = 'TleafCnd'
  )
}

\arguments{
  \item{exdf_obj}{
    An \code{exdf} object representing data from a Licor gas exchange
    measurement system.
  }

  \item{temperature_response_parameters}{
    A list, where each element describes the temperature response of a parameter
    value. The name of each element must be the name of the parameter. Each
    element must be a list itself, whose named elements must include the type of
    temperature response function to use (\code{type}), thee units of the
    parameter (\code{units}), and the values of necessary temperature response
    parameters. See below for more details.
  }

  \item{tleaf_column_name}{
    The name of the column in \code{exdf_obj} that contains the leaf temperature
    in units of \code{degrees C}.
  }
}

\details{
  Some key photosynthetic parameters are known to vary with temperature
  according to well-established temperature response functions such as the
  Arrhenius equation. The \code{calculate_temperature_response} function can be
  used to calculate such temperature-dependent parameter values at leaf
  temperature.

  Depending on the \code{type} value supplied in each element of
  \code{temperature_response_parameters}, one of several possible functions will
  be used to calculate the temperature response:
  \itemize{
    \item When \code{type} is \code{'Arrhenius'}, the
          \code{\link{calculate_temperature_response_arrhenius}} function will
          be used.
    \item When \code{type} is \code{'Gaussian'}, the
          \code{\link{calculate_temperature_response_gaussian}} function will
          be used.
  }

  Values of \code{type} are not case-sensitive.
}

\value{
  An \code{exdf} object based on \code{exdf_obj} that includes one new column
  for each element of \code{temperature_response_parameters}, where the
  temperature-dependent values of these new columns are determined using the
  temperature values specified by the \code{tleaf_column_name} column. The
  category of each of these new columns is \code{calculate_temperature_response}
  to indicate that they were created using this function.
}

\examples{
# Read an example Licor file included in the PhotoGEA package
licor_file <- read_gasex_file(
  PhotoGEA_example_file_path('ball_berry_1.xlsx')
)

# In this example we will calculate temperature-dependent values of two
# parameters:
#
# - The `Kc` parameter (in units of `micromol mol^(-1)`) will be calculated
#   using an Arrhenius function with scaling constant `c` = 38.05 and activation
#   energy `Ea` = 79.43 kJ / mol.
#
# - The `Jmax` parameter (in units of `micromol m^(-2) s^(-1)) will be
#   using a Gaussian function with optimal temperature `t_opt` = 43 degrees C
#   and width `sigma` = 16 degrees C.
#
# So the `temperature_response_parameters` list will contain two elements,
# defined as follows:

trp <- list(
  Kc = list(
    type = 'Arrhenius',
    c = 38.05,
    Ea = 79.43,
    units = 'micromol mol^(-1)'
  ),
  Jmax = list(
    type = 'Gaussian',
    optimum_rate = 4,
    t_opt = 43,
    sigma = 16,
    units = 'micromol m^(-2) s^(-1)'
  )
)

# Now we can calculate the values of Kc and Jmax at the measured leaf
# temperatures recorded in the log file
licor_file <- calculate_temperature_response(licor_file, trp)

licor_file$units$Kc      # View the units of the new `Kc` column
licor_file$categories$Kc # View the category of the new `Kc` column
licor_file[,'Kc']        # View the values of the new `Kc` column

licor_file$units$Jmax      # View the units of the new `Jmax` column
licor_file$categories$Jmax # View the category of the new `Jmax` column
licor_file[,'Jmax']        # View the values of the new `Jmax` column
}

\concept{exdf}
