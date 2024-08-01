\name{c3_arrhenius_sharkey}

\docType{data}

\alias{c3_arrhenius_sharkey}

\title{C3 Arrhenius parameters from Sharkey et al.}

\description{
  Arrhenius parameters describing the temperature response of important C3
  photosynthetic parameters, intended to be passed to the
  \code{\link{calculate_arrhenius}} function.
}

\usage{c3_arrhenius_sharkey}

\format{
  List with 6 named elements that each represent a variable whose
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
  }

  In turn, each of these elements is a list with 3 named elements:
  \itemize{
    \item \code{c}: the (dimensionless) Arrhenius scaling factor.
    \item \code{Ea}: the activation energy in \code{kJ / mol}.
    \item \code{units}: the units of the corresponding variable.
  }
}

\source{
  Some of these parameters (\code{RL}, \code{Vcmax}, and \code{J}) are
  normalized to their values at 25 degrees C. These variables include
  \code{_norm} in their names to indicate this.

  The remaining parameters (\code{Gamma_star}, \code{Kc}, and \code{Ko}) are not
  normalized because they are assumed to not vary significantly between species.

  Arrhenius parameters were obtained from Sharkey et al. (2007). In this
  publication, gas concentrations are expressed as partial pressures (in
  \code{Pa} or \code{kPa}) rather than mole fractions (\code{micromol / mol} or
  \code{mmol / mol}). However, for consistency with
  \code{\link{c3_arrhenius_bernacchi}}, here we prefer to use mole fractions.

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
