\name{jmax_temperature_param_bernacchi}

\docType{data}

\alias{jmax_temperature_param_bernacchi}

\title{Jmax-related temperature response parameters from Bernacchi et al.}

\description{
  Parameters describing the temperature response of Jmax-related photosynthetic
  parameters, intended to be passed to the
  \code{\link{calculate_temperature_response}} function.
}

\usage{jmax_temperature_param_bernacchi}

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
  Polynomial coefficients were obtained from Bernacchi et al. (2003).
  Here, we use the values determined from plants grown at 25 degrees C
  (Table 2). The coefficients given in the paper are used to calculate the
  values of \code{alpha_j} and \code{theta_j} at leaf temperature. Here we
  normalize by the values of \code{alpha_j} and \code{theta_j} at 25 degrees C,
  which are 0.6895 and 0.97875, respectively.

  References:

  \itemize{
    \item Bernacchi, C. J., Pimentel, C. & Long, S. P. "In vivo temperature
          response functions of parameters required to model RuBP-limited
          photosynthesis" Plant, Cell & Environment 26, 1419–1430 (2003)
          [\doi{10.1046/j.0016-8025.2003.01050.x}].
  }
}

\keyword{datasets}

\concept{temperature_response_parameters}
