\name{calculate_gas_properties}

\alias{calculate_gas_properties}

\title{calculate_gas_properties}

\description{
  Calculates gas properties that are typically not included in Licor files. This
  function can accomodate alternative column names for the variables taken from
  the Licor file in case they change at some point in the future. This function
  also checks the units of each required column and will produce an error if any
  units are incorrect.
}

\usage{
  calculate_gas_properties(
    licor_exdf,
    a_column_name = 'A',
    ca_column_name = 'Ca',
    deltapcham_column_name = 'DeltaPcham',
    e_column_name = 'E',
    gbw_column_name = 'gbw',
    gsw_column_name = 'gsw',
    h2o_s_column_name = 'H2O_s',
    pa_column_name = 'Pa',
    tleaf_column_name = 'TleafCnd'
  )
}

\arguments{
  \item{licor_exdf}{
    An \code{exdf} object representing data from a Licor gas exchange
    measurement system.
  }

  \item{a_column_name}{
    The name of the column in \code{licor_exdf} that contains the net
    assimilation in \code{micromol m^(-2) s^(-1)}.
  }

  \item{ca_column_name}{
    The name of the column in \code{licor_exdf} that contains the ambient CO2
    concentration in the chamber in \code{micromol mol^(-1)}.
  }

  \item{deltapcham_column_name}{
    The name of the column in \code{licor_exdf} that contains the chamber
    overpressure in \code{kPa}.
  }

  \item{e_column_name}{
    The name of the column in \code{licor_exdf} that contains the transpiration
    rate in \code{mol m^(-2) s^(-1)}.
  }

  \item{gbw_column_name}{
    The name of the column in \code{licor_exdf} that contains the boundary layer
    conductance to water vapor in \code{mol m^(-2) s^(-1)}.
  }

  \item{gsw_column_name}{
    The name of the column in \code{licor_exdf} that contains the stomatal
    conductance to water vapor in \code{mol m^(-2) s^(-1)}.
  }

  \item{h2o_s_column_name}{
    The name of the column in \code{licor_exdf} that contains the sample cell
    H2O concentration in \code{mmol mol^(-1)}.
  }

  \item{pa_column_name}{
    The name of the column in \code{licor_exdf} that contains the atmospheric
    pressure in \code{kPa}.
  }

  \item{tleaf_column_name}{
    The name of the column in \code{licor_exdf} that contains the leaf
    temperature in \code{degrees C}.
  }
}

\details{
  By default, a Licor file provides the following gas concentrations and
  conductances:
  \itemize{
    \item Water vapor conductance to diffusion through the stomata (\code{gsw}).
    \item Water vapor conductance to diffusion through the boundary layer
          (\code{gbw}).
    \item Water vapor conductance to diffusion from the leaf's intercellular
          spaces to the ambient air; in other words, the total conductance to
          water vapor (\code{gtw}).
    \item Water vapor concentration in the sample cell (\code{H2O_s}).
    \item CO2 conductance to diffusion from the leaf's intercellular spaces to
          the ambient air; in other words, the total conductance to CO2
          (\code{gtc}).
    \item CO2 concentration in the sample cell, corrected for any chamber leaks
         (\code{Ca}).
    \item CO2 concentration in the leaf's intercellular spaces (\code{Ci}).
  }

  However, it is sometimes helpful to know the "missing" conductances and
  concentrations, for example, when calculating mesophyll conductances or
  Ball-Berry parameters. This function adds these missing values, along with a
  few related water vapor properties:
  \itemize{
    \item Water vapor concentration at the sample surface (\code{H2O_surf}).
    \item Water vapor concentration in the leaf's intercellular spaces
          (\code{H2O_i}).
    \item Saturation water vapor pressure at the leaf temperature
          (\code{SVPleaf}).
    \item Relative humidity at the leaf surface (\code{RHleaf}).
    \item CO2 conductance to diffusion through the stomata (\code{gsc}).
    \item CO2 conductance to diffusion through the boundary layer (\code{gbc}).
    \item CO2 concentration at the leaf surface (\code{Cs}).
  }

  \emph{Equations used for these calculations}

  The equations used to calculate these quantities can be found in the Licor
  Li-6800 manual (Appendix C), which relies heavily on Appendix 2 of the
  following paper: von Caemmerer, S. & Farquhar, G. D. "Some relationships
  between the biochemistry of photosynthesis and the gas exchange of leaves"
  Planta \strong{153}, 376–387 (1981) [\doi{10.1007/BF00384257}]

  Equation C-79 in the Licor manual describes the total flow of water vapor from
  the leaf interior to the ambient air using \code{gtw}, \code{H2O_i},
  \code{H2O_s}, and the transpiration rate \code{E}:

  (1) \code{gtw = E * (1000 - (H2O_i + H2O_s) / 2) / (H2O_i - H2O_s)}

  In steady-state conditions, the flux of H2O molecules across any portion of
  the gas flow is identical to \code{E}, so we can also apply this equation to
  the flow of water vapor from the leaf surface to the ambient air:

  (2) \code{gbw = E * (1000 - (H2O_surf + H2O_s) / 2) / (H2O_surf - H2O_s)}

  Equation (2) can be solved for \code{H2O_surf}:

  (3) \code{H2O_surf = (E * (1000 - H2O_s / 2) + gbw * H2O_s) / (gbw + E / 2)}

  Equation C-70 in the Licor manual describes how to calculate saturation water
  vapor pressure from air temperature. At the leaf surface, the air temperature
  should be the same as the leaf temperature (\code{Tleaf}; in degrees C), so we
  can determine \code{SVPleaf} using Equation C-70 as follows:

  (4) \code{SVPleaf = 0.6135 * e^((17.502 * Tleaf) / (240.97 + Tleaf))}

  For gas exchange measurements, we assume that water vapor is saturated in the
  leaf's intecellular spaces, so we can determine \code{H2O_i} from
  \code{SVPleaf} and the relationship between partial pressure and molar gas
  concentration:

  (5) \code{H2O_i = SVPleaf / Pcham = SVPleaf / (Pa + deltaPcham)}

  where \code{Pcham} is th total pressure in the sample chamber, \code{Pa} is
  the atmospheric pressure, and \code{deltaPcham} is the chamber overpressure.
  These are related by \code{Pcham = Pa + deltaPcham}.

  The relative humidity at the leaf surface \code{RHleaf} can be determined from
  \code{H2O_surf} and \code{SVPleaf} using the definitions of relative humidity
  and partial pressure:

  (6) \code{RHleaf = Pwl / SVPleaf = H2O_surf * (Pa + deltaPcham) / SVPleaf}

  where \code{Pwl}, the partial pressure of H2O at the leaf surface, is given by
  \code{H2O_surf * Pcham}.

  The CO2 conductances through the stomata and boundary layer can be determined
  from the corresponding H2O conductances using the ratios of molecular
  diffusivities for the two molecules, as explained in the vicinty of Equation
  C-106 in the Licor manual:

  (7) \code{gsc = gsw / 1.6}

  (8) \code{gbc = gbw / 1.37}

  Equation C-105 in the Licor manual describes the flow of CO2 from the ambient
  air to the intercellular spaces:

  (9) \code{C_i = ((gtc - E / 2) * Ca - A) / (gtc + E / 2)}

  where we have replaced \code{C_s} (the CO2 concentration in the sample
  chamber) with \code{Ca} for clarity. In steady state conditions, the flows of
  H2O and CO2 are identical to \code{E} and \code{A}, respectively, so we can
  also apply this equation to the flow of CO2 from the ambient air to the leaf
  surface:

  (10) \code{Csurface = ((gbc - E / 2) * Ca - A) / (gbc + E / 2)}

  This function uses Equations (3)-(8) and (10) to calculate the desired values.
}

\value{
  An \code{exdf} object based on \code{licor_exdf} that includes the following
  additional columns, calculated as described above: \code{H2O_surf},
  \code{SVPleaf}, \code{H2O_i}, \code{RHleaf}, \code{gsc}, \code{gbc}, and
  \code{Csurface}. The category for each of these new columns is
  \code{calculate_gas_properties} to indicate that they were created using this
  function.
}

\examples{
# Read an example Licor file included in the PhotoGEA package and calculate
# additional gas properties.
licor_file <- read_licor_file(
  system.file("extdata", "ball_berry_1.xlsx", package = "PhotoGEA"),
  c(3, 5, 7, 9, 11, 13), 14, 15, 16, 17, 'time'
)

licor_file <- calculate_gas_properties(licor_file)

licor_file$units$RHleaf      # View the units of the new `RHleaf` column
licor_file$categories$RHleaf # View the category of the new `RHleaf` column
licor_file[,'RHleaf']        # View the values of the new `RHleaf` column
}

\concept{exdf}
