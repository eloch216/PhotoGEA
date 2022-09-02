\name{calculate_c3_assimilation}

\alias{calculate_c3_assimilation}

\title{calculate_c3_assimilation}

\description{
  Calculates C3 assimilation rates using equations from S. von Caemmerer,
  "Biochemical Models of Leaf Photosynthesis" (CSIRO Publishing, 2000)
  [\doi{10.1071/9780643103405}]. This function can accomodate alternative column
  names for the variables taken from Licor files in case they change at some
  point in the future. This function also checks the units of each required
  column and will produce an error if any units are incorrect.
}

\usage{
  calculate_c3_assimilation(
  exdf_obj,
  cc_column_name,
  tleaf_column_name,
  pa_column_name,
  deltapcham_column_name,
  gmc_column_name,
  PTR_FUN,
  POc,
  TPU,
  J,
  Rd,
  Vcmax
)
}

\arguments{
  \item{exdf_obj}{An \code{exdf} object.}

  \item{cc_column_name}{
    The name of the column in \code{exdf_obj} that contains the chloroplastic
    CO2 concentration (typically \code{'Cc'}). The values of column are often
    calculated using the \code{\link{calculate_cc}} function.
  }

  \item{tleaf_column_name}{
    The name of the column in \code{exdf_obj} that contains the leaf
    temperature (typically \code{'TleafCnd'}).
  }

  \item{pa_column_name}{
    The name of the column in \code{exdf_obj} that contains the atmospheric
    pressure (typically \code{'Pa'}).
  }

  \item{deltapcham_column_name}{
    The name of the column in \code{exdf_obj} that contains the chamber
    overpressure (typically \code{'DeltaPcham'}).
  }

  \item{gmc_column_name}{
    The name of the column in \code{exdf_obj} that contains the mesophyll
    conductance to CO2 (typically \code{'gmc'}). The values of this column are
    often set using the \code{\link{set_variable}} function.
  }

  \item{PTR_FUN}{
    A photosynthetic temperature response function such as
    \code{photosynthesis_TRF(temperature_response_parameters_Bernacchi)}.
  }

  \item{POc}{
    The partial pressure of O2 in the chloroplast, expressed in \code{microbar}.
    This is often assumed to be the ambient value. For air measurements, this
    would be 21 percent O2, which is about 210000 microbar at standard
    atmospheric pressure. For low oxygen measurements, this would be 2 percent
    O2, which is about 20000 microbar.
  }

  \item{TPU}{
    The maximum rate of triphosphate utilization, expressed in
    \code{micromol m^(-2) s^(-1)}.
  }

  \item{J}{
    The electron transport rate at 25 degrees C, expressed in
    \code{micromol m^(-2) s^(-1)}. Note that this is _not_ \code{Jmax}, and in
    general will depend on the incident photosynthetically active flux density.
  }

  \item{Rd}{
    The respiration rate at 25 degrees C, expressed in
    \code{micromol m^(-2) s^(-1)}.
  }

  \item{Vcmax}{
    The maximum rate of rubisco carboxylation at 25 degrees C, expressed in
    \code{micromol m^(-2) s^(-1)}.
  }
}

\details{
  This function follows Chapter 2 from S. von Caemmerer's textbook, referenced
  above. In this model, the net CO2 assimilation rate \code{An} is given by the
  smaller of the RuBP-regeneration-limited rate \code{Ac}, the
  electron-transport-limited rate \code{Aj}, and the phosphate-limited rate
  \code{Ap}, as expressed by Equation 2.27:

  \code{An = min(Ac, Aj, Ap)}

  The RuBP-regeneration-limited rate is calculated using Equation 2.20:

  \code{Ac = (PCc - Gamma_star) * Vcmax_tl / (PCc + Kc * (1.0 + POc / Ko)) - Rd_tl}

  where \code{PCc} is the partial pressure of CO2 in the chloroplast,
  \code{Gamma_star} is the value of \code{Pcc} at which no assimilation occurs
  in the absence of respiration, \code{Vcmax_tl} is the
  leaf-temperature-adjusted maximum rate of rubisco carboxylation, \code{Kc} is
  the Michaelis-Menten constant for CO2, \code{Ko} is the Michaelis-Menten
  constant for O2, and \code{Rd_tl} is the leaf-temperature-adjusted respiration
  rate. Here, the temperature corrections are applied to \code{Vcmax} and
  \code{Rd} using \code{PTR_FUN}, and temperature-dependent values of \code{Kc},
  \code{Ko}, and \code{Gamma_star} are also determined from \code{PTR_FUN}.
  \code{PCc} is calculated from the chloroplastic CO2 concentration \code{Cc}
  using the total pressure (ambient pressure + chamber overpressure).

  The electron-transport-limited rate is calculated using Equation 2.23:

  \code{Aj = (PCc - Gamma_star) * J_tl / (4 * PCc + 8 * Gamma_star) - Rd_tl}

  where \code{J_tl} is the leaf-temperature-adjusted value of \code{J}, as
  determined using \code{PTR_FUN}.

  The phosphate-limited rate is calculated using Equation 2.26:

  \code{Ap = (PCc - Gamma_star) * (3 * TPU) / (PCc - (1 + 3 * alpha / 2) * Gamma_star) - Rd_tl}

  where \code{alpha} is the fraction of glycolate carbon _not_ returned to the
  chloroplast; here, as in many places, \code{alpha} is fixed to be 0.
}

\value{
  An \code{exdf} object with the following columns, calculated as described
  above: \code{Kc}, \code{Ko}, \code{Gamma_star}, \code{Vcmax_tl}, \code{Rd_tl},
  \code{J_tl}, \code{Ac}, \code{Aj}, \code{Ap}, and \code{An}. The category for
  each of these new columns is \code{calculate_c3_assimilation} to indicate that
  they were created using this function.
}

\examples{
# Simulate a C3 A-Cc curve with specified leaf temperature and photosynthetic
# parameters and plot the net assimilation rate along with the different
# enzyme-limited rates
npts <- 101

inputs <- exdf(data.frame(
  Cc = seq(0, 600, length.out = npts),
  Tleaf = 30,
  Pa = 100,
  DeltaPcham = 0.1,
  gmc = 1.0
))

inputs <- document_variables(
  inputs,
  c('', 'Cc',         'micromol mol^(-1)'),
  c('', 'Tleaf',      'degrees C'),
  c('', 'Pa',         'kPa'),
  c('', 'DeltaPcham', 'kPa'),
  c('', 'gmc',        'mol m^(-2) s^(-1) bar^(-1)')
)

assim <- calculate_c3_assimilation(
  inputs,
  'Cc',
  'Tleaf',
  'Pa',
  'DeltaPcham',
  'gmc',
  photosynthesis_TRF(temperature_response_parameters_Sharkey),
  POc = 210000,
  TPU = 12,
  J = 150,
  Rd = 1,
  Vcmax = 120
)

lattice::xyplot(
  Ac + Aj + Ap + An ~ Cc,
  data = cbind(inputs, assim)$main_data,
  type = 'l',
  grid = TRUE,
  auto = TRUE,
  xlab = paste0('Chloroplast CO2 concentration (', inputs$units$Cc, ')'),
  ylab = paste0('Assimilation rate (', assim$units$An, ')')
)
}

\concept{exdf}
