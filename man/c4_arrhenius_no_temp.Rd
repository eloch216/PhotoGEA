\name{c4_arrhenius_no_temp}

\docType{data}

\alias{c4_arrhenius_no_temp}

\title{C4 Arrhenius parameters with no temperature dependence}

\description{
  Arrhenius parameters that specify a flat temperature response (in other words,
  no dependence on temperature) for important C4 photosynthetic parameters,
  intended to be passed to the \code{\link{calculate_arrhenius}} function.
}

\usage{c4_arrhenius_no_temp}

\format{
  List with 9 named elements that each represent a variable whose
  temperature-dependent value can be calculated using the Arrhenius equation:
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
    \item \code{gmc}: The mesophyll conductance to CO2 diffusion.
  }

  In turn, each of these elements is a list with 3 named elements:
  \itemize{
    \item \code{c}: the (dimensionless) Arrhenius scaling factor.
    \item \code{Ea}: the activation energy in \code{kJ / mol}.
    \item \code{units}: the units of the corresponding variable.
  }
}

\source{
  Some of these parameters (\code{Vcmax}, \code{Vpmax}, and \code{RL}) are
  normalized to their values at 25 degrees C. These variables include
  \code{_norm} in their names to indicate this.

  The remaining parameters (\code{Kc}, \code{Ko}, \code{Kp}, \code{gamma_star},
  \code{ao}, and \code{gmc}) are not normalized because they are assumed to not
  vary significantly between species.

  Here, the activation energy values (\code{Ea}) are all set to 0, which means
  that the values will not depend on temperature. The Arrhenius scaling factors
  \code{c} are chosen to reproduce the parameter values at 25 degrees C as
  specified in von Caemmerer (2021).
  (See \code{\link{c4_arrhenius_von_caemmerer}}.)

  References:

  \itemize{
    \item von Caemmerer, S. "Updating the steady-state model of C4
          photosynthesis" Journal of Experimental Botany 72, 6003–6017 (2021)
          [\doi{10.1093/jxb/erab266}].
  }
}

\keyword{datasets}

\concept{temperature_response_parameters}
