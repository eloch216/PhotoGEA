\name{c4_temperature_param_flat}

\docType{data}

\alias{c4_temperature_param_flat}

\title{C4 temperature response parameters for a flat response}

\description{
  Parameters that specify a flat temperature response (in other words, no
  dependence on temperature) for important C4 photosynthetic parameters,
  intended to be passed to the \code{\link{calculate_temperature_response}}
  function.
}

\usage{c4_temperature_param_flat}

\format{
  List with 10 named elements that each represent a variable whose
  temperature-dependent value can be calculated using either an Arrhenius or
  Gaussian equation:
  \itemize{
    \item \code{Vcmax_norm}: The maximum rate of rubisco carboxylation
          (\code{Vcmax}) normalized to its value at 25 degrees C.
    \item \code{Vpmax_norm}: The maximum rate of PEP carboxylase activity
          (\code{Vpmax}) normalized to its value at 25 degrees C.
    \item \code{RL_norm}: The respiration rate (\code{RL}) normalized to the
          value of \code{Vcmax} at 25 degrees C.
    \item \code{Kc}: The Michaelis-Menten constant for rubisco carboxylation.
    \item \code{Ko}: The Michaelis-Menten constant for rubisco oxygenation.
    \item \code{Kp}: The Michaelis-Menten constant of PEP carboxylase.
    \item \code{gamma_star}: Half the reciprocal of rubisco specificity.
    \item \code{ao}: The ratio of solubility and diffusivity of O2 to CO2.
    \item \code{gmc_norm}: The mesophyll conductance to CO2 diffusion normalized
          to its value at 25 degrees C.
    \item \code{J_norm}: The electron transport rate \code{J} normalized to its
          value at 25 degrees C.
  }

  Each of these is a list with 4 named elements:
  \itemize{
    \item \code{type}: the type of temperature response (\code{'Arrhenius'})
    \item \code{c}: the (dimensionless) Arrhenius scaling factor.
    \item \code{Ea}: the activation energy in \code{kJ / mol}.
    \item \code{units}: the units of the corresponding variable.
  }
}

\source{
  Some of these parameters (\code{Vcmax}, \code{Vpmax}, \code{RL}, \code{gmc},
  and \code{J}) are normalized to their values at 25 degrees C. These
  variables include \code{_norm} in their names to indicate this.

  The remaining parameters (\code{Kc}, \code{Ko}, \code{Kp}, \code{gamma_star},
  \code{ao}, and \code{gmc}) are not normalized because they are assumed to not
  vary significantly between species.

  Here, the activation energy values (\code{Ea}) are all set to 0, which means
  that the values will not depend on temperature. The Arrhenius scaling factors
  \code{c} are chosen to reproduce the parameter values at 25 degrees C as
  specified in von Caemmerer (2021).
  (See \code{\link{c4_temperature_param_vc}}.)

  References:

  \itemize{
    \item von Caemmerer, S. "Updating the steady-state model of C4
          photosynthesis" Journal of Experimental Botany 72, 6003–6017 (2021)
          [\doi{10.1093/jxb/erab266}].
  }
}

\keyword{datasets}

\concept{temperature_response_parameters}
