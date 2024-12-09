\name{calculate_temperature_response_gaussian}

\alias{calculate_temperature_response_gaussian}

\title{Calculate temperature-dependent values using Gaussian equations}

\description{
  Calculate leaf-temperature-dependent values of various parameters using
  Gaussian equations. It is rare for users to call this function directly;
  instead, it is used internally by
  \code{\link{calculate_temperature_response}}.
}

\usage{
  calculate_temperature_response_gaussian(
    exdf_obj,
    gaussian_parameters,
    tleaf_column_name = 'TleafCnd'
  )
}

\arguments{
  \item{exdf_obj}{
    An \code{exdf} object representing data from a Licor gas exchange
    measurement system.
  }

  \item{gaussian_parameters}{
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

  \code{rate = optimal_rate * exp(-(T - T_opt)^2 / sigma^2)}

  where \code{optimal_rate} is the highest rate which occurs at the optimal
  temperature \code{T_opt}, \code{T} is the current temperature, and
  \code{sigma} represents the "width" of the peak. More technically, it can be
  described as the difference in temperature away from the optimal value at
  which the rate falls to 37 percent (\code{1/e}) of its maximum.

  In \code{calculate_temperature_response_gaussian}, the optimal rate
  (\code{optimal_rate}), optimal temperature (\code{t_opt}),
  width (\code{sigma}), and units (\code{units}) for a variable must be
  specified as elements of a list, which itself is a named element of
  \code{gaussian_parameters}. For example, if a variable called
  \code{Jmax} has \code{optimal_rate = 1}, \code{t_opt = 43}, {sigma = 26}, and
  units of \code{micromol mol^(-1)}, the \code{gaussian_parameters} argument
  could be specified as follows:
  \code{list(Jmax = list(optimal_rate = 1, t_opt = 43, sigma = 26, units =
  'micromol mol^(-1)'))}.

  It is rare to specify these parameters directly; instead, it is more typical
  to use one of the pre-set values such as those included in
  \code{\link{c4_temperature_param_vc}}.
}

\value{
  An \code{exdf} object based on \code{exdf_obj} that includes one new column
  for each element of \code{gaussian_parameters}, where the
  temperature-dependent values of these new columns are determined using the
  temperature values specified by the \code{tleaf_column_name} column. The
  category of each of these new columns is
  \code{calculate_temperature_response_gaussian} to indicate that they were
  created using this function.
}

\examples{
# Read an example Licor file included in the PhotoGEA package
licor_file <- read_gasex_file(
  PhotoGEA_example_file_path('ball_berry_1.xlsx')
)

licor_file <- calculate_temperature_response_gaussian(
  licor_file,
  list(J_norm = c4_temperature_param_vc$J_norm)
)

licor_file$units$J_norm      # View the units of the new `J_norm` column
licor_file$categories$J_norm # View the category of the new `J_norm` column
licor_file[,'J_norm']        # View the values of the new `J_norm` column
}

\concept{exdf}
