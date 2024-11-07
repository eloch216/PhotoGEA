\name{c3_temperature_param_bernacchi}

\docType{data}

\alias{c3_temperature_param_bernacchi}

\title{C3 temperature response parameters from Bernacchi et al.}

\description{
  Parameters describing the temperature response of important C3 photosynthetic
  parameters, intended to be passed to the
  \code{\link{calculate_temperature_response}} function.
}

\usage{c3_temperature_param_bernacchi}

\format{
  List with 8 named elements that each represent a variable whose
  temperature-dependent value can be calculated using the Arrhenius equation:
  \itemize{
    \item \code{RL_norm}: The respiration rate (\code{RL}) normalized to its
          value at 25 degrees C.
    \item \code{Vcmax_norm}: The maximum rate of rubisco carboxylation
          (\code{Vcmax}) normalized to its value at 25 degrees C.
    \item \code{Vomax_norm}: The maximum rate of rubisco oxygenation
          (\code{Vomax}) normalized to its value at 25 degrees C.
    \item \code{Gamma_star}: The value of CO2 partial pressure in the
          chloroplast at which no assimilation occurs in the absence of
          respiration.
    \item \code{Kc}: The Michaelis-Menten constant for rubisco carboxylation.
    \item \code{Ko}: The Michaelis-Menten constant for rubisco oxygenation.
    \item \code{J_norm}: The electron transport rate (\code{J}) normalized to
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
  Some of these parameters (\code{RL}, \code{Vcmax}, \code{J}, and \code{Tp})
  are normalized to their values at 25 degrees C. \code{Vomax} is normalized to
  the value of \code{Vcmax} at 25 degrees C. These variables include
  \code{_norm} in their names to indicate this.

  The remaining parameters (\code{Gamma_star}, \code{Kc}, and \code{Ko}) are not
  normalized because they are assumed to not vary significantly between species.

  Arrhenius parameters for \code{J} were obtained from Bernacchi et al. (2003).
  Here, we use the values determined from chlorophyll fluorescence measured from
  plants grown at 25 degrees C (Table 1). Although Bernacchi et al. (2003)
  reports values of \code{Jmax}, here we assume that both \code{Jmax} and the
  light-dependent values of \code{J} follow the same temperature response
  function and refer to it as \code{J} for compatibility with
  \code{\link{c3_temperature_param_sharkey}}.

  The Arrhenius parameters for the other variables were obtained from Bernacchi
  et al. (2001).

  References:

  \itemize{
    \item Bernacchi, C. J., Singsaas, E. L., Pimentel, C., Jr, A. R. P. & Long,
    S. P. "Improved temperature response functions for models of Rubisco-limited
    photosynthesis" Plant, Cell & Environment 24, 253–259 (2001)
    [\doi{10.1111/j.1365-3040.2001.00668.x}].

    \item Bernacchi, C. J., Pimentel, C. & Long, S. P. "In vivo temperature
    response functions of parameters required to model RuBP-limited
    photosynthesis" Plant, Cell & Environment 26, 1419–1430 (2003)
    [\doi{10.1046/j.0016-8025.2003.01050.x}].
  }
}

\keyword{datasets}

\concept{temperature_response_parameters}
