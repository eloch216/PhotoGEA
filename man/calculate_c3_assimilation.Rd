\name{calculate_c3_assimilation}

\alias{calculate_c3_assimilation}

\title{calculate_c3_assimilation}

\description{
  Calculates C3 assimilation rates based on the Farquhar-von-Caemmerer-Berry
  model. This function can accomodate alternative colum names for the variables
  taken from Licor files in case they change at some point in the future. This
  function also checks the units of each required column and will produce an
  error if any units are incorrect.
}

\usage{
  calculate_c3_assimilation(
    exdf_obj,
    TPU,
    J,
    Rd,
    Vcmax,
    POc = 210000,
    curvature = 1,
    cc_column_name = 'Cc',
    pa_column_name = 'Pa',
    deltapcham_column_name = 'DeltaPcham',
    kc_column_name = 'Kc',
    ko_column_name = 'Ko',
    gamma_star_column_name = 'Gamma_star',
    vcmax_norm_column_name = 'Vcmax_norm',
    rd_norm_column_name = 'Rd_norm',
    j_norm_column_name = 'J_norm',
    perform_checks = TRUE,
    return_exdf = TRUE
  )
}

\arguments{
  \item{exdf_obj}{An \code{exdf} object.}

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

  \item{POc}{
    The partial pressure of O2 in the chloroplast, expressed in \code{microbar}.
    This is often assumed to be the ambient value. For air measurements, this
    would be 21 percent O2, which is about 210000 microbar at standard
    atmospheric pressure. For low oxygen measurements, this would be 2 percent
    O2, which is about 20000 microbar.
  }

  \item{curvature}{
    A dimensionless quadratic curvature parameter between 0 and 1 that sets the
    degree of co-limitation between \code{Ac}, \code{Aj}, and \code{Ap}. A value
    of 1 indicates no co-limitation.
  }

  \item{cc_column_name}{
    The name of the column in \code{exdf_obj} that contains the
    chloroplastic CO2 concentration in \code{micromol mol^(-1)}.
  }

  \item{pa_column_name}{
    The name of the column in \code{exdf_obj} that contains the atmospheric
    pressure in \code{kPa}.
  }

  \item{deltapcham_column_name}{
    The name of the column in \code{exdf_obj} that contains the chamber
    overpressure in \code{kPa}.
  }

  \item{kc_column_name}{
    The name of the column in \code{exdf_obj} that contains the
    Michaelis-Menten constant for rubisco carboxylation in
    \code{micromol mol^(-1)}.
  }

  \item{ko_column_name}{
    The name of the column in \code{exdf_obj} that contains the
    Michaelis-Menten constant for rubisco oxygenation in
    \code{mmol mol^(-1)}.
  }

  \item{gamma_star_column_name}{
    The name of the column in \code{exdf_obj} that contains the
    \code{Gamma_star} values in \code{micromol mol^(-1)}.
  }

  \item{vcmax_norm_column_name}{
    The name of the column in \code{exdf_obj} that contains the normalized
    \code{Vcmax} values (with units of
    \code{normalized to Vcmax at 25 degrees C}).
  }

  \item{rd_norm_column_name}{
    The name of the column in \code{exdf_obj} that contains the normalized
    \code{Rd} values (with units of \code{normalized to Rd at 25 degrees C}).
  }

  \item{j_norm_column_name}{
    The name of the column in \code{exdf_obj} that contains the normalized
    \code{J} values (with units of \code{normalized to J at 25 degrees C}).
  }

  \item{perform_checks}{
    A logical value indicating whether to check units for the required columns.
    This should almost always be \code{TRUE}. The option to disable these checks
    is only intended to be used when \code{\link{fit_c3_aci}} calls this
    function, since performing these checks many times repeatedly slows down the
    fitting procedure.
  }

  \item{return_exdf}{
    A logical value indicating whether to return an \code{exdf} object. This
    should almost always be \code{TRUE}. The option to return a vector is mainly
    intended to be used when  \code{\link{fit_c3_aci}} calls this function,
    since creating an \code{exdf} object to return will slows down the fitting
    procedure.
  }
}

\details{
  This function generally follows Chapter 2 from S. von Caemmerer (2000), which
  provides equations for calculating the rubisco-limited assimilation rate
  \code{Ac}, the electron-transport-limited assimilation rate \code{Aj}, and the
  phosphate-limited assimilation rate \code{Ap}.

  The rubisco-regeneration-limited rate is calculated using Equation 2.20:

  \code{Ac = (PCc - Gamma_star) * Vcmax_tl / (PCc + Kc * (1.0 + POc / Ko)) - Rd_tl}

  where \code{PCc} is the partial pressure of CO2 in the chloroplast,
  \code{Gamma_star} is the value of \code{Pcc} at which no assimilation occurs
  in the absence of respiration, \code{Vcmax_tl} is the
  leaf-temperature-adjusted maximum rate of rubisco carboxylation, \code{Kc} is
  the Michaelis-Menten constant for CO2, \code{Ko} is the Michaelis-Menten
  constant for O2, and \code{Rd_tl} is the leaf-temperature-adjusted respiration
  rate. Here, the temperature corrections are applied to \code{Vcmax} and
  \code{Rd} using the normalized values \code{Vcmax_norm} and \code{Rd_norm}.
  \code{PCc} is calculated from the chloroplastic CO2 concentration \code{Cc}
  using the total pressure (ambient pressure + chamber overpressure).

  The electron-transport-limited rate is calculated using Equation 2.23:

  \code{Aj = (PCc - Gamma_star) * J_tl / (4 * PCc + 8 * Gamma_star) - Rd_tl}

  where \code{J_tl} is the leaf-temperature-adjusted value of \code{J},
  determined from the normalized value \code{J_norm}.

  The phosphate-limited rate is calculated using Equation 2.26:

  \code{Ap = (PCc - Gamma_star) * (3 * TPU) / (PCc - (1 + 3 * alpha / 2) * Gamma_star) - Rd_tl}

  where \code{alpha} is the fraction of glycolate carbon _not_ returned to the
  chloroplast; here, as in many places, \code{alpha} is fixed to be 0.

  In the textbook, the net CO2 assimilation rate \code{An} is given by the
  smaller of the rubisco-limited rate \code{Ac}, the electron-transport-limited
  rate \code{Aj}, and the phosphate-limited rate \code{Ap}, as expressed by
  Equation 2.27:

  \code{An = min(Ac, Aj, Ap)}

  This equation determines the assimilation rate for a single chloroplast;
  however, we often wish to calculate assimilation values at the leaf level. In
  this case, some degree of co-limitation between \code{Ac}, \code{Aj}, and
  \code{Ap} should be allowed because chloroplast CO2 concentrations are
  typically not uniform across the leaf. Here we accomplish this co-limitation
  via a quadratic curvature parameter, as discussed in Collatz et al. (1990)
  and Collatz et al. (1991). In particular, we solve the following equation for
  \code{Acj}, the assimilation rate co-limited by rubisco and RuBP saturation:

  \code{curvature^2 * Acj - Acj * (Ac + Aj) + Ac * Aj}

  Then we solve a second equation for \code{An}, the assimilation rate
  co-limited by rubisco, RuBP regeneration, and phosphate utilization:

  \code{curvature^2 * An - An * (Acj + Ap) + Acj * Ap}

  Note that when \code{curvature} is 1, \code{An} will be equal to \code{min(Ac,
  Aj, Ap)} as in the version that applies for a single chloroplast.

  References:
  \itemize{
    \item von Caemmerer, S. "Biochemical Models of Leaf Photosynthesis" (CSIRO
          Publishing, 2000) [\doi{10.1071/9780643103405}].
    \item Collatz, G. J., Ball, J. T., Grivet, C. & Berry, J. A. "Physiological
          and environmental regulation of stomatal conductance, photosynthesis
          and transpiration: a model that includes a laminar boundary layer."
          Agricultural and Forest Meteorology 54, 107–136 (1991)
          [doi:10.1016/0168-1923(91)90002-8].
    \item Collatz, G. J., Berry, J. A., Farquhar, G. D. & Pierce, J. "The
          relationship between the Rubisco reaction mechanism and models of
          photosynthesis." Plant, Cell & Environment 13, 219–225 (1990)
          [doi:10.1111/j.1365-3040.1990.tb01306.x].
  }
}

\value{
  The return value depends on the value of \code{return_exdf}:
  \itemize{
    \item If \code{return_exdf} is \code{TRUE}, the return value is an
          \code{exdf} object with the following columns, calculated as described
          above: \code{TPU}, \code{Vcmax_tl}, \code{Rd_tl}, \code{J_tl},
          \code{Ac}, \code{Aj}, \code{Ap}, and \code{An}. The category for each
          of these new columns is \code{calculate_c3_assimilation} to indicate
          that they were created using this function.
    \item If \code{return_exdf} is \code{FALSE}, the return value is a numeric
          vector containing the calculated values of \code{An}.
  }
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
  DeltaPcham = 0.1
))

inputs <- document_variables(
  inputs,
  c('', 'Cc',         'micromol mol^(-1)'),
  c('', 'Tleaf',      'degrees C'),
  c('', 'Pa',         'kPa'),
  c('', 'DeltaPcham', 'kPa')
)

inputs <- calculate_arrhenius(inputs, c3_arrhenius_sharkey, 'Tleaf')

assim <- calculate_c3_assimilation(inputs, 12, 150, 1, 120)

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
