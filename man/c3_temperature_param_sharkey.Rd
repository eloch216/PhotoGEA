\name{c3_temperature_param_sharkey}

\docType{data}

\alias{c3_temperature_param_sharkey}

\title{C3 temperature response parameters from Sharkey et al.}

\description{
  Parameters describing the temperature response of important C3 photosynthetic
  parameters, intended to be passed to the
  \code{\link{calculate_temperature_response}} function.
}

\usage{c3_temperature_param_sharkey}

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

  Response parameters were obtained from Sharkey et al. (2007). In this
  publication, gas concentrations are expressed as partial pressures (in
  \code{Pa} or \code{kPa}) rather than mole fractions (\code{micromol / mol} or
  \code{mmol / mol}). However, for consistency with
  \code{\link{c3_temperature_param_bernacchi}}, here we prefer to use mole
  fractions.

  To convert a concentration expressed as a partial pressure (\code{P}; in
  \code{Pa}) to a concentration expressed as a mole fraction (\code{C}; in
  \code{micromol / mol}), we need a value for atmospheric pressure; we will use
  the typical value of \code{101325 Pa}. Then \code{C = P / 101325 * 1e6} or
  \code{C = P * cf}, where \code{cf = 1e6 / 101325} is a conversion factor. The
  same correction can be used to convert \code{kPa} to \code{mmol / mol}. The
  value of \code{cf} can be accessed using \code{PhotoGEA:::c_pa_to_ppm}.

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
