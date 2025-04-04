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
  List with 11 named elements that each represent a variable whose
  temperature-dependent value can be calculated using an Arrhenius equation or a
  polynomial equation:
  \itemize{
    \item \code{Gamma_star_at_25}: The value of chloroplastic CO2 concentration
          at which CO2 gains from Rubisco carboxylation are exactly balanced by
          CO2 losses from Rubisco oxygenation (\code{Gamma_star}) at 25
          degrees C.

    \item \code{Gamma_star_norm}: \code{Gamma_star} normalized to its value at
          25 degrees C.

    \item \code{gmc_norm}: The mesophyll conductance to CO2 diffusion
          (\code{gmc}) normalized to its value at 25 degrees C.

    \item \code{J_norm}: The electron transport rate (\code{J}) normalized to
          its value at 25 degrees C.

    \item \code{Kc_at_25}: The Michaelis-Menten constant for rubisco
          carboxylation (\code{Kc}) at 25 degrees C.

    \item \code{Kc_norm}: \code{Kc} normalized to its value at 25 degrees C.

    \item \code{Ko_at_25}: The Michaelis-Menten constant for rubisco oxygenation
          (\code{Ko}) at 25 degrees C.

    \item \code{Ko_norm}: \code{Ko} normalized to its value at 25 degrees C.

    \item \code{RL_norm}: The rate of non-photorespiratory CO2 release in the
          light (\code{RL}) normalized to its value at 25 degrees C.

    \item \code{Tp_norm}: The maximum rate of triose phosphate utilization
          (\code{Tp}) normalized to its value at 25 degrees C.

    \item \code{Vcmax_norm}: The maximum rate of rubisco carboxylation
          (\code{Vcmax}) normalized to its value at 25 degrees C.
  }

  In turn, each of these elements is a list with at least 2 named elements:
  \itemize{
    \item \code{type}: the type of temperature response
    \item \code{units}: the units of the corresponding variable.
  }
}

\source{
  Many of these parameters are normalized to their values at 25 degrees C. These
  variables include \code{_norm} in their names to indicate this.

  Here, the activation energy values (\code{Ea}) are all set to 0, which means
  that the values will not depend on temperature. Some parameters are specified
  at 25 degrees C; these values were obtained from Sharkey et al. (2007). (See
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
