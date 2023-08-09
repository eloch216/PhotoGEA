\name{calculate_ternary_correction}

\alias{calculate_ternary_correction}

\title{Calculate ternary correction factor}

\description{
  Calculates the ternary correction factor \code{t} that is used in many carbon
  isotope discrimination calculations.
}

\usage{
  calculate_ternary_correction(
    exdf_obj,
    ci_column_name = 'Ci',
    co2_s_column_name = 'CO2_s',
    csurface_column_name = 'Csurface',
    e_column_name = 'E',
    gtc_column_name = 'gtc'
  )
}

\arguments{
  \item{exdf_obj}{
    An \code{exdf} object containing photosynthetic gas exchange data.
  }

  \item{ci_column_name}{
    The name of the column in \code{exdf_obj} that contains the intercellular
    CO2 concentration in \code{micromol mol^(-1)}.
  }

  \item{co2_s_column_name}{
    The name of the column in \code{exdf_obj} that contains the sample line
    (incoming air) CO2 concentration in \code{micromol mol^(-1)}.
  }

  \item{csurface_column_name}{
    The name of the column in \code{exdf_obj} that contains the CO2
    concentration at the leaf surface in \code{micromol mol^(-1)}.
  }

  \item{e_column_name}{
    The name of the column in \code{exdf_obj} that contains the leaf
    transpiration rate in \code{mol m^(-2) s^(-1)}.
  }

  \item{gtc_column_name}{
    The name of the column in \code{exdf_obj} that contains the total
    conductance to CO2 diffusion across the boundary layer and stomata in series
    in \code{mol m^(-2) s^(-1)}.
  }
}

\details{
  During photosynthetic gas exchange, there are separate fluxes of CO2 and H2O
  flowing in and out of the leaf. These gases interact with each other and with
  air, forming a ternary mixture. These interactions must be taken into account
  when modeling carbon isotope discrimination. Typically this is done via
  \code{t}, a ternary correction factor first introduced by Farquhar and
  Cernusak (2012). Here we calculate \code{t} as described in Equations 9 and 10
  from Ubierna et al. (2018):

  \code{t = alpha_ac * E / (2 * g_ac)}

  and

  \code{a_bar = (a_b * (C_a - C_s) + a_s * (C_s - C_i)) / (C_a - C_i)},

  where \code{E} is the transpiration rate, \code{g_ac} is the total conductance
  to CO2 diffusion across the boundary layer and stomata in series, \code{a_bar}
  is the weighted fractionation across the boundary layer and stomata in series,
  \code{a_b} is the fractionation during diffusion through the boundary layer,
  \code{a_s} is the fractionation during diffusion through the stomata,
  \code{C_a} is the ambient CO2 concentration (in wet air), \code{C_s} is the
  CO2 concentration (in wet air) at the leaf surface, and \code{C_i} is the CO2
  concentration (in wet air) in the intercellular spaces.

  \code{alpha_ac} is the overall fractionation during diffusion through air;
  \code{alpha_ac} and \code{a_bar} are related according to an un-numbered
  equation in Ubierna et al. (2018) that appears just after Equation 9:

  \code{alpha_ac = 1 + a_bar}

  References:

  Farquhar, G. D. and Cernusak, L. A. "Ternary effects on the gas exchange of
  isotopologues of carbon dioxide." Plant, Cell & Environment 35, 1221–1231
  (2012) [\doi{10.1111/j.1365-3040.2012.02484.x}].

  Ubierna, N., Holloway-Phillips, M.-M. and Farquhar, G. D. "Using Stable Carbon
  Isotopes to Study C3 and C4 Photosynthesis: Models and Calculations." in
  Photosynthesis: Methods and Protocols (ed. Covshoff, S.) 155–196 (Springer,
  2018) [\doi{10.1007/978-1-4939-7786-4_10}].
}

\value{
  An \code{exdf} object based on \code{exdf_obj} that includes values of
  \code{t}, \code{a_bar}, and \code{alpha_ac} calculated as described above.
  The category of each new column is \code{calculate_ternary_correction} to
  indicate that it was created using this function.
}

\concept{exdf}
