\name{calculate_gm_ubierna}

\alias{calculate_gm_ubierna}

\title{Calculate mesophyll conductance to CO2 diffusion}

\description{
  Calculates mesophyll conductance to CO2 diffusion (\code{gmc}) from combined
  gas exchange and isotope discrimination measurements as described in Ubierna
  et al. (2018). This function can accomodate alternative colum names for the
  variables taken from \code{exdf_obj}; it also checks the units of each
  required column and will produce an error if any units are incorrect.
}

\usage{
  calculate_gm_ubierna(
    exdf_obj,
    e = -3,
    f = 11,
    a_bar_column_name = 'a_bar',
    a_column_name = 'A',
    ci_column_name = 'Ci',
    co2_s_column_name = 'CO2_s',
    csurface_column_name = 'Csurface',
    delta_c13_r_column_name = 'delta_C13_r',
    delta_obs_tdl_column_name = 'Delta_obs_tdl',
    gamma_star_column_name = 'Gamma_star',
    rd_column_name = 'Rd',
    total_pressure_column_name = 'total_pressure',
    t_column_name = 't'
)
}

\arguments{
  \item{exdf_obj}{An \code{exdf} object.}

  \item{e}{
    The isotopic fractionation during day respiration in \code{ppt}.
  }

  \item{f}{
    The isotopic fractionation during photorespiration in \code{ppt}.
  }

  \item{a_bar_column_name}{
    The name of the column in \code{exdf_obj} that contains the weighted
    isotopic fractionation across the boundary layer and stomata in \code{ppt}.
    Values of \code{a_bar} are typically calculated using
    \code{\link{calculate_ternary_correction}}.
  }

  \item{a_column_name}{
    The name of the column in \code{exdf_obj} that contains the net CO2
    assimilation rate in \code{micromol m^(-2) s^(-1)}.
  }

  \item{ci_column_name}{
    The name of the column in \code{exdf_obj} that contains the intercellular
    CO2 concentration in \code{micromol mol^(-1)}.
  }

  \item{co2_s_column_name}{
    The name of the column in \code{exdf_obj} that contains the CO2
    concentration in the sample line (outgoing air) in \code{micromol mol^(-1)}.
  }

  \item{csurface_column_name}{
    The name of the column in \code{exdf_obj} that contains the CO2
    concentration at the leaf surface in \code{micromol mol^(-1)}. Values of
    \code{Csurface} are typically calculated using
    \code{\link{calculate_gas_properties}}.
  }

  \item{delta_c13_r_column_name}{
    The name of the column in \code{exdf_obj} that contains the CO2 isotope
    ratio in the reference line (incoming air) in \code{ppt}.
  }

  \item{delta_obs_tdl_column_name}{
    The name of the column in \code{exdf_obj} that contains the observed isotope
    discrimination values in \code{ppt}.
  }

  \item{gamma_star_column_name}{
    The name of the column in \code{exdf_obj} that contains the CO2 compensation
    point in the absence of day respiration in \code{micromol mol^(-1)}.  Values
    of \code{Gamma_star} are typically calculated using
    \code{\link{calculate_gamma_star}} or \code{\link{calculate_arrhenius}}.
  }

  \item{rd_column_name}{
    The name of the column in \code{exdf_obj} that contains the rate of day
    respiration in \code{micromol m^(-2) s^(-1)}.
  }

  \item{total_pressure_column_name}{
    The name of the column in \code{exdf_obj} that contains the total pressure
    in \code{bar}.
  }

  \item{t_column_name}{
    The name of the column in \code{exdf_obj} that contains the ternary
    correction factor (\code{dimensionless}). Values of \code{t} are typically
    calculated using \code{\link{calculate_ternary_correction}}
  }
}

\details{
  This function uses the comprehensive model for photosynthetic discrimination
  against 13C in C3 plants, as described in Ubierna et al. (2018). In
  particular, the following equations from that source are implemented in the
  code:

  \itemize{
    \item Isotopic fractionation due to day respiration (\code{e_prime}) is
          calculated using Equations 28 and 30.
    \item Isotopic discrimination due to photorespiration (\code{Delta_f}), due
          to day respiration (\code{Delta_e}), and that would occur if Ci = Cc
          in the absence of any respiratory fractionation (\code{Delta_i}) are
          calculated using Equations 34, 33, and 31, respectively.
    \item Mesophyll conductance to CO2 diffusion (\code{gmc}) is calculated
          using Equation 44. This equation is broken up into two factors called
          \code{Delta_difference} and \code{equation_top} which are separately
          returned in the output from \code{calculate_gm_ubierna}.
  }

  References:

  Ubierna, N., Holloway-Phillips, M.-M. and Farquhar, G. D. "Using Stable Carbon
  Isotopes to Study C3 and C4 Photosynthesis: Models and Calculations." in
  Photosynthesis: Methods and Protocols (ed. Covshoff, S.) 155–196 (Springer,
  2018) [\doi{10.1007/978-1-4939-7786-4_10}].
}

\value{
  An \code{exdf} object based on \code{exdf_obj} that includes the following
  additional columns, calculated as described above: \code{e_prime},
  \code{Delta_i}, \code{Delta_e}, \code{Delta_f}, \code{Delta_difference},
  \code{equation_top}, and \code{gmc}. The category for each of these new
  columns is \code{calculate_gm_ubierna} to indicate that they were created
  using this function.
}

\concept{exdf}
