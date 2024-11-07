\name{c4_temperature_param_vc}

\docType{data}

\alias{c4_temperature_param_vc}

\title{C4 temperature response parameters from von Caemmerer}

\description{
  Temperature response parameters describing the temperature response of
  important C4 photosynthetic parameters, intended to be passed to the
  \code{\link{calculate_temperature_response}} function.
}

\usage{c4_temperature_param_vc}

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
    \item \code{gmc}: The mesophyll conductance to CO2 diffusion.
    \item \code{Jmax_norm}: The maximum rate of RuBP regneration \code{Jmax}
          normalized to its value at its optimal temperature.
  }

  The \code{Jmax_norm} parameter is calculated using a Gaussian function and
  hence its corresponding list element is itself a list with 4 named elements:
  \itemize{
    \item \code{type}: the type of temperature response (\code{'Gaussian'})
    \item \code{optimum_rate}: the largest value this parameter can take.
    \item \code{t_opt}: the temperature where the optimum occurs in
          \code{degrees C}.
    \item \code{sigma}: the width of the Gaussian in \code{degrees C}.
    \item \code{units}: the units of the corresponding variable.
  }

  Each of the remaining elements is a list with 4 named elements:
  \itemize{
    \item \code{type}: the type of temperature response (\code{'Arrhenius'})
    \item \code{c}: the (dimensionless) Arrhenius scaling factor.
    \item \code{Ea}: the activation energy in \code{kJ / mol}.
    \item \code{units}: the units of the corresponding variable.
  }
}

\source{
  Some of these parameters (\code{Vcmax}, \code{Vpmax}, \code{RL}, and
  \code{Jmax}) are normalized to their values at 25 degrees C. These variables
  include \code{_norm} in their names to indicate this.

  The remaining parameters (\code{Kc}, \code{Ko}, \code{Kp}, \code{gamma_star},
  \code{ao}, and \code{gmc}) are not normalized because they are assumed to not
  vary significantly between species.

  Here, the Arrhenius scaling factors (\code{c}; dimensionless) and activation
  energy values (\code{Ea}; kJ / mol) are obtained from von Caemmerer (2021). In
  that publication, the overall scaling for each parameter is specified by its
  value at 25 degrees C; the scaling factors are determined from this
  information as described in the documentation for
  \code{\link{calculate_temperature_response_arrhenius}}.

  The Gaussian parameters (\code{t_opt} and \code{sigma}) for \code{Jmax_norm}
  are also obtained from von Caemmerer (2021). The value of \code{optimum_rate}
  is chosen such that \code{Jmax_norm} is equal to 1 at a temperature of 25
  degrees C.

  References:

  \itemize{
    \item von Caemmerer, S. "Updating the steady-state model of C4
          photosynthesis" Journal of Experimental Botany 72, 6003–6017 (2021)
          [\doi{10.1093/jxb/erab266}].
  }
}

\keyword{datasets}

\concept{temperature_response_parameters}
