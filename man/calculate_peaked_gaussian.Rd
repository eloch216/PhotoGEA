\name{calculate_peaked_gaussian}

\alias{calculate_peaked_gaussian}

\title{Calculate temperature-dependent values using Gaussian equations}

\description{
  Calculate leaf-temperature-dependent values of various parameters using
  Gaussian equations. This function can accomodate alternative names for the
  leaf temperature column. It also checks the units of the leaf temperature
  column and will produce an error of the units are not correct.
}

\usage{
  calculate_peaked_gaussian(
    exdf_obj,
    peaked_gaussian_parameters,
    tleaf_column_name = 'TleafCnd'
  )
}

\arguments{
  \item{exdf_obj}{
    An \code{exdf} object representing data from a Licor gas exchange
    measurement system.
  }

  \item{peaked_gaussian_parameters}{
    A list of named lists. Each list element should describe the optimal
    temperature in \code{degrees C} (\code{t_opt}), the "width" in
    \code{degrees C} (\code{sigma}), and the units (\code{units}) for a variable
    that follows a peaked Gaussian temperature dependence. The name of each list
    element should be the corresponding name of the variable.
  }

  \item{tleaf_column_name}{
    The name of the column in \code{exdf_obj} that contains the leaf temperature
    in units of \code{degrees C}.
  }
}

\details{
  A Gaussian equation is sometimes used to model the temperature dependence of
  a biochemical rate parameter. Typically this is expressed by

  (1) \code{rate = optimal_rate * exp(-(T - T_opt)^2 / sigma^2)}

  where \code{optimal_rate} is the highest rate which occurs at the optimal
  temperature \code{T_opt}, \code{T} is the current temperature, and
  \code{sigma} represents the "width" of the peak. More technically, it can be
  described as the difference in temperature away from the optimal value at
  which the rate falls to 37 percent (\code{1/e}) of its maximum.

  The \code{calculate_peaked_gaussian} function only calculates the temperature
  response part of equation (1):

  (2) \code{exp(-(T - T_opt)^2 / sigma^2)}

  Thus, it returns dimensionless scaling factors that are normalized to the
  optimal rate.

  In \code{calculate_peaked_gaussian}, the optimal temperature (\code{t_opt}),
  width (\code{sigma}), and units (\code{units}) for a variable must be specified
  as elements of a list, which itself is a named element of
  \code{peaked_gaussian_parameters}. For example, if a variable called
  \code{Jmax} has \code{t_opt = 43}, {sigma = 26}, and units of
  \code{micromol mol^(-1)}, the \code{peaked_gaussian_parameters} argument could
  be specified as follows:
  \code{list(Jmax_norm = list(t_opt = 43, sigma = 26, units = 'normalized to Jmax at its optimal temperature'))}.
  It is rare to directly specify these parameters directly; instead, it is more
  typical to use one of the pre-set values such as
  \code{\link{c4_peaked_gaussian_von_caemmerer}}.
}

\value{
  An \code{exdf} object based on \code{exdf_obj} that includes one new column
  for each element of \code{peaked_gaussian_parameters}, where the
  temperature-dependent values of these new columns are determined using the
  temperature values specified by the \code{tleaf_column_name} column. The
  category of each of these new columns is \code{calculate_peaked_gaussian} to
  indicate that they were created using this function.
}

\examples{
# Read an example Licor file included in the PhotoGEA package
licor_file <- read_gasex_file(
  PhotoGEA_example_file_path('ball_berry_1.xlsx')
)

licor_file <- calculate_peaked_gaussian(licor_file, c4_peaked_gaussian_von_caemmerer)

licor_file$units$Jmax_norm      # View the units of the new `Jmax_norm` column
licor_file$categories$Jmax_norm # View the category of the new `Jmax_norm` column
licor_file[,'Jmax_norm']        # View the values of the new `Jmax_norm` column
}

\concept{exdf}
