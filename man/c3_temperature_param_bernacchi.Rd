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
  List with 12 named elements that each represent a variable whose
  temperature-dependent value can be calculated using an Arrhenius equation,
  Johnson-Eyring-Williams equation, or a polynomial equation:
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

    \item \code{Vomax_norm}: The maximum rate of rubisco oxygenation
          (\code{Vomax}) normalized to \code{Vcmax} at 25 degrees C.
  }

  In turn, each of these elements is a list with at least 2 named elements:
  \itemize{
    \item \code{type}: the type of temperature response
    \item \code{units}: the units of the corresponding variable.
  }
}

\source{
  Many of these parameters are normalized to their values at 25 degrees C.
  \code{Vomax} is normalized to the value of \code{Vcmax} at 25 degrees C. These
  variables include \code{_norm} in their names to indicate this.

  Arrhenius parameters for \code{J} were obtained from Bernacchi et al. (2003).
  Here, we use the values determined from chlorophyll fluorescence measured from
  plants grown at 25 degrees C (Table 1). Although Bernacchi et al. (2003)
  reports values of \code{Jmax}, here we assume that both \code{Jmax} and the
  light-dependent values of \code{J} follow the same temperature response
  function and refer to it as \code{J} for compatibility with
  \code{\link{c3_temperature_param_sharkey}}.

  Johnson-Eyring-Williams parameters for \code{gmc} were obtained from Bernacchi
  et al. (2002).

  The Bernacchi papers from the early 2000s do not specify a temperature
  response for \code{Tp}, so we instead use the Johnson-Eyring-Williams
  response from Sharkey et al. (2007). Another option would be to use a flat
  temperature response; in other words, to assume that \code{Tp} is constant
  with temperature. This could be done with the following code, which takes the
  flat response parameters from \code{\link{c3_temperature_param_flat}}:
  \code{within(c3_temperature_param_bernacchi, {Tp_norm = c3_temperature_param_flat$Tp_norm})}

  The Arrhenius parameters for the other variables were obtained from Bernacchi
  et al. (2001).

  References:

  \itemize{
    \item Bernacchi, C. J., Singsaas, E. L., Pimentel, C., Jr, A. R. P. & Long,
          S. P. "Improved temperature response functions for models of
          Rubisco-limited photosynthesis" Plant, Cell & Environment 24, 253–259
          (2001) [\doi{10.1111/j.1365-3040.2001.00668.x}].

    \item Bernacchi, C. J., Portis, A. R., Nakano, H., von Caemmerer, S. & Long,
          S. P. "Temperature Response of Mesophyll Conductance. Implications for
          the Determination of Rubisco Enzyme Kinetics and for Limitations to
          Photosynthesis in Vivo" Plant Physiology 130, 1992–1998 (2002)
          [\doi{10.1104/pp.008250}].

    \item Bernacchi, C. J., Pimentel, C. & Long, S. P. "In vivo temperature
          response functions of parameters required to model RuBP-limited
          photosynthesis" Plant, Cell & Environment 26, 1419–1430 (2003)
          [\doi{10.1046/j.0016-8025.2003.01050.x}].

    \item Sharkey, T. D., Bernacchi, C. J., Farquhar, G. D. & Singsaas, E. L.
          "Fitting photosynthetic carbon dioxide response curves for C3 leaves"
          Plant, Cell & Environment 30, 1035–1040 (2007)
          [\doi{10.1111/j.1365-3040.2007.01710.x}].
  }
}

\keyword{datasets}

\concept{temperature_response_parameters}
