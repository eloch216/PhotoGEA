\name{c3_temperature_param_flat}

\docType{data}

\alias{c3_temperature_param_flat}

\title{C3 temperature response parameters for a flat response}

\description{
  Parameters that specify a flat temperature response (in other words,
  no dependence on temperature) for important C3 photosynthetic parameters,
  intended to be passed to the \code{\link{calculate_temperature_response}}
  function.
}

\usage{c3_temperature_param_flat}

\format{
  List with 8 named elements that each represent a variable whose
  temperature-dependent value can be calculated using the Arrhenius equation:
  \itemize{
    \item \code{Kc}: The Michaelis-Menten constant for rubisco carboxylation.
    \item \code{Ko}: The Michaelis-Menten constant for rubisco oxygenation.
    \item \code{Gamma_star}: The value of CO2 partial pressure in the
          chloroplast at which no assimilation occurs in the absence of
          respiration.
    \item \code{Vcmax_norm}: The maximum rate of rubisco carboxylation
          (\code{Vcmax}) normalized to its value at 25 degrees C.
    \item \code{J_norm}: The electron transport rate (\code{J}) normalized to
          its value at 25 degrees C.
    \item \code{RL_norm}: The respiration rate (\code{RL}) normalized to its
          value at 25 degrees C.
    \item \code{gmc_norm}: The mesophyll conductance (\code{gmc}) normalized to
          its value at 25 degrees C.
    \item \code{Tp_norm}: The maximum rate of triose phosphate utilization
          (\code{Tp}) normalized to its value at 25 degrees C.
  }

  In turn, each of these elements is a list with 4 named elements:
  \itemize{
    \item \code{type}: the type of temperature response (\code{'Arrhenius'})
    \item \code{c}: the (dimensionless) Arrhenius scaling factor.
    \item \code{Ea}: the activation energy in \code{kJ / mol}.
    \item \code{units}: the units of the corresponding variable.
  }
}

\source{
  Some of these parameters (\code{RL}, \code{Vcmax}, \code{J}, \code{Tp}, and
  \code{gmc}) are normalized to their values at 25 degrees C. These variables
  include \code{_norm} in their names to indicate this.

  The remaining parameters (\code{Gamma_star}, \code{Kc}, and \code{Ko}) are not
  normalized because they are assumed to not vary significantly between species.

  Here, the activation energy values (\code{Ea}) are all set to 0, which means
  that the values will not depend on temperature. The Arrhenius scaling factors
  \code{c} are chosen to reproduce the parameter values at 25 degrees C as
  specified in Sharkey et al. (2007). (See
  \code{\link{c3_temperature_param_sharkey}}.)

  References:

  \itemize{
    \item Sharkey, T. D., Bernacchi, C. J., Farquhar, G. D. & Singsaas, E. L.
          "Fitting photosynthetic carbon dioxide response curves for C3 leaves"
          Plant, Cell & Environment 30, 1035–1040 (2007)
          [\doi{10.1111/j.1365-3040.2007.01710.x}].
  }
}

\keyword{datasets}

\concept{temperature_response_parameters}
