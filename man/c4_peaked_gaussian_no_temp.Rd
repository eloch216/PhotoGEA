\name{c4_peaked_gaussian_no_temp}

\docType{data}

\alias{c4_peaked_gaussian_no_temp}

\title{C4 peaked Gaussian parameters with no temperature response}

\description{
  Peaked Gaussian parameters that describe a flat temperature response of
  important C4 photosynthetic parameters, intended to be passed to the
  \code{\link{calculate_peaked_gaussian}} function.
}

\usage{c4_peaked_gaussian_no_temp}

\format{
  List with 1 named element that represents a variable whose
  temperature-dependent value can be calculated using the peaked Gaussian
  equation:
  \itemize{
    \item \code{Jmax_norm}: The maximum rate of RuBP regneration \code{Jmax}
          normalized to its value at its optimal temperature.
  }

  In turn, each of these elements is a list with 3 named elements:
  \itemize{
    \item \code{t_opt}: the optimal temperature (where the variable takes its
          highest value) in \code{degrees C}.
    \item \code{sigma}: the difference in temperature from the optimal value at
          which the variable's value falls to 37 percent (\code{1/e}) of its
          highest value, in \code{degrees C}.
    \item \code{units}: the units of the variable.
  }
}

\source{
  Here, the \code{sigma} values are set to infinity (\code{Inf}), which ensures
  that the normalization factors are 1 for all temperatures.

  References:

  \itemize{
    \item von Caemmerer, S. "Updating the steady-state model of C4
          photosynthesis" Journal of Experimental Botany 72, 6003–6017 (2021)
          [\doi{10.1093/jxb/erab266}].
  }
}

\keyword{datasets}

\concept{temperature_response_parameters}
