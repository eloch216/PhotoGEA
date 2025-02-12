\name{calculate_temperature_response_polynomial}

\alias{calculate_temperature_response_polynomial}

\title{Calculate temperature-dependent values using polynomial equations}

\description{
  Calculate leaf-temperature-dependent values of various parameters using
  polynomial equations. It is rare for users to call this function directly;
  instead, it is used internally by
  \code{\link{calculate_temperature_response}}.
}

\usage{
  calculate_temperature_response_polynomial(
    exdf_obj,
    polynomial_parameters,
    tleaf_column_name = 'TleafCnd'
  )
}

\arguments{
  \item{exdf_obj}{
    An \code{exdf} object representing data from a Licor gas exchange
    measurement system.
  }

  \item{polynomial_parameters}{
    A list of named lists. Each list element should describe the polynomial
    coefficients (\code{coef}) and units (\code{units}) for a variable that
    follows a polynomial temperature dependence. The name of each list element
    should be the corresponding name of the variable.
  }

  \item{tleaf_column_name}{
    The name of the column in \code{exdf_obj} that contains the leaf temperature
    in units of \code{degrees C}.
  }
}

\details{
  Polynomial equations are often used to calculate the temperature dependence
  of the rates of chemical reactions. For example, a second-order polynomial
  could be given as follows:

  (1) \code{rate = R_0 + R_1 * T + R_2 * T^2}

  where \code{R_0}, \code{R_1}, and \code{R_2} are the zeroth, first, and second
  order coefficients and \code{T} is the temperature. Higher order polynomials
  can also be defined, where an order-\code{N} polynomial is given by

  (2) \code{rate = R_0 + R_1 * T + R_2 * T^2 + ... + R_N * T^N}

  In general, an order-\code{N} polynomial has \code{N} coefficients, although
  some of them may be zero.

  In \code{calculate_temperature_response_polynomial}, the coefficients
  (\code{coef}) and units (\code{units}) for a variable must be specified as
  elements of a list, which itself is a named element of
  \code{polynomial_parameters}. The coefficients must be specified as a numeric
  vector, where the \code{i}th element represents the \code{i}th coefficient.
  For example, if a dimensionless variable called \code{theta} is calculated
  according to \code{theta = 0.352 + 0.022 * T - 3.4e-4 * T^2}, the
  \code{polynomial_parameters} argument could be supplied as follows:
  \code{list(theta = list(coef = c(0.352, 0.022, -3.4e-4), units = 'dimensionless'))}.

  It is rare to directly specify the polynomial parameters; instead, it is more
  typical to use one of the pre-set values such as those included in
  \code{\link{jmax_temperature_param_bernacchi}}.
}

\value{
  An \code{exdf} object based on \code{exdf_obj} that includes one new column
  for each element of \code{polynomial_parameters}, where the
  temperature-dependent values of these new columns are determined using the
  temperature values specified by the \code{tleaf_column_name} column. The
  category of each of these new columns is
  \code{calculate_temperature_response_polynomial} to indicate that they were
  created using this function.
}

\examples{
# Read an example Licor file included in the PhotoGEA package
licor_file <- read_gasex_file(
  PhotoGEA_example_file_path('ball_berry_1.xlsx')
)

licor_file <- calculate_temperature_response_polynomial(
  licor_file,
  list(theta = list(coef = c(0.352, 0.022, -3.4e-4), units = 'dimensionless'))
)

licor_file$units$theta      # View the units of the new `theta` column
licor_file$categories$theta # View the category of the new `theta` column
licor_file[,'theta']        # View the values of the new `theta` column
}

\concept{exdf}
