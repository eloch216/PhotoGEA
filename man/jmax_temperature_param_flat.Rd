\name{jmax_temperature_param_flat}

\docType{data}

\alias{jmax_temperature_param_flat}

\title{Jmax-related temperature response parameters from Bernacchi et al.}

\description{
  Parameters that describe a flat temperature response (in other words,
  no dependence on temperature) for Jmax-related photosynthetic
  parameters, intended to be passed to the
  \code{\link{calculate_temperature_response}} function.
}

\usage{jmax_temperature_param_flat}

\format{
  List with 2 named elements that each represent a variable whose
  temperature-dependent value can be calculated using a polynomial equation:
  \itemize{
    \item \code{alpha_j_norm}: The apparent quantum efficiency of electron
          transport (\code{alpha_j}) normalized to its value at 25 degrees C.
    \item \code{theta_j_norm}: The empirical curvature parameter normalized to
          its value at 25 degrees C.
  }

  In turn, each of these elements is a list with 3 named elements:
  \itemize{
    \item \code{type}: the type of temperature response.
    \item \code{coef}: the polynomial coefficients.
    \item \code{units}: the units of the corresponding variable.
  }
}

\source{
  Here, the polynomial coefficients (\code{coef}) are all set to 1, speciying a
  zeroth-order polynomial equal to 1, which means that the values will not
  depend on temperature.
}

\keyword{datasets}

\concept{temperature_response_parameters}
