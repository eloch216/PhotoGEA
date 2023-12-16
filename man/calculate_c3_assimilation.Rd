\name{calculate_c3_assimilation}

\alias{calculate_c3_assimilation}

\title{Calculate C3 assimilation rates}

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
    atp_use = 4.0,
    nadph_use = 8.0,
    alpha = 0.0,
    curvature_cj = 1.0,
    curvature_cjp = 1.0,
    cc_column_name = 'Cc',
    gamma_star_column_name = 'Gamma_star',
    j_norm_column_name = 'J_norm',
    kc_column_name = 'Kc',
    ko_column_name = 'Ko',
    rd_norm_column_name = 'Rd_norm',
    total_pressure_column_name = 'total_pressure',
    vcmax_norm_column_name = 'Vcmax_norm',
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

  \item{atp_use}{
    The number of ATP molecules used per C3 cycle.
  }

  \item{nadph_use}{
    The number of NADPH molecules used per C3 cycle.
  }

  \item{alpha}{
    A dimensionless parameter greater than or equal to 0 and less than or equal
    to 1 that represents the fraction of glycolate carbon not returned to the
    chloroplast.
  }

  \item{curvature_cj}{
    A dimensionless quadratic curvature parameter greater than or equal to 0 and
    less than or equal to 1 that sets the degree of co-limitation between
    \code{Wc} and \code{Wj}. A value of 1 indicates no co-limitation.
  }

  \item{curvature_cjp}{
    A dimensionless quadratic curvature parameter greater than or equal to 0 and
    less than or equal to 1 that sets the degree of co-limitation between
    \code{Wcj} and \code{Wp}. A value of 1 indicates no co-limitation.
  }

  \item{cc_column_name}{
    The name of the column in \code{exdf_obj} that contains the
    chloroplastic CO2 concentration in \code{micromol mol^(-1)}.
  }

  \item{gamma_star_column_name}{
    The name of the column in \code{exdf_obj} that contains the
    \code{Gamma_star} values in \code{micromol mol^(-1)}.
  }

  \item{j_norm_column_name}{
    The name of the column in \code{exdf_obj} that contains the normalized
    \code{J} values (with units of \code{normalized to J at 25 degrees C}).
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

  \item{rd_norm_column_name}{
    The name of the column in \code{exdf_obj} that contains the normalized
    \code{Rd} values (with units of \code{normalized to Rd at 25 degrees C}).
  }

  \item{total_pressure_column_name}{
    The name of the column in \code{exdf_obj} that contains the total pressure
    in \code{bar}.
  }

  \item{vcmax_norm_column_name}{
    The name of the column in \code{exdf_obj} that contains the normalized
    \code{Vcmax} values (with units of
    \code{normalized to Vcmax at 25 degrees C}).
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
    since creating an \code{exdf} object to return will slow down the fitting
    procedure.
  }
}

\details{
  This function generally follows the Farquhar-von-Caemmerer-Berry model as
  described in Lochocki and McGrath (in preparation). In this formulation, the
  steady-state net CO2 assimilation rate \code{An} is calculated according to

  \code{An = (1 - Gamma_star / PCc) * Vc - Rd},

  where \code{Gamma_star} is the CO2 compensation point in the absence of day
  respiration, \code{PCc} is the partial pressure of CO2 in the chloroplast,
  \code{Vc} is the RuBP carboxylation rate, and \code{Rd} is the rate of day
  respiration.

  The model considers three potential values of \code{Vc} that correspond to
  limitations set by three different processes: Rubisco activity, RuBP
  regeneration, and triose phopsphate utilization (TPU). The Rubisco-limited
  carboxylation rate \code{Wc} is given by

  \code{Wc = PCc * Vcmax / (PCc + Kc * (1.0 + POc / Ko))},

  where \code{Vcmax} is the maximum rate of Rubisco carboxylation, \code{Kc} is
  the Michaelis-Menten constant for CO2, \code{Ko} is the Michaelis-Menten
  constant for O2, and \code{POc} is the partial pressure of O2 in the
  chloroplast.

  The RuBP-regeneration-limited carboxylation rate \code{Wj} is given by

  \code{Wj = PCc * J / (atp_use * PCc + nadph_use * Gamma_star)},

  where \code{J} is the RuBP regeneration rate. Here we have allowed
  \code{atp_use} and \code{nadph_use} to be variables rather than taking fixed
  values (as they do in many sources). This is necessary because not all
  descriptions of the FvCB model use the same values.

  The TPU-limited carboxylation rate is given by

  \code{Wp = PCc * 3 * TPU / (PCc - Gamma_star * (1 + 3 * alpha))},

  where \code{TPU} is the maximum rate of triose phosphate utilization and
  \code{0 <= alpha <= 1} is the fraction of glycolate carbon not returned to the
  chloroplast. This equation only applies when
  \code{PCc > Gamma_star * (1 + 3 * alpha)}; for smaller values of \code{PCc},
  TPU cannot limit the RuBP carboxylation rate and \code{Wp = Inf}.

  The actual carboxylation rate is often chosen to be the smallest of the three
  potential rates:

  \code{Vc = min{Wc, Wj, Wp}}.

  However, it is also possible to allow co-limitation between the processes by
  using quadratic mixing equations. In other words, the carboxylation rate
  co-limitated by Rubisco activity and RuBP regeneration \code{Wcj} is given by
  the smaller root of the following quadratic equation:

  \code{curvature_cj * Wcj^2 - Wcj * (Wc + Wj) + Wc * Wj = 0},

  where \code{curvature_cj} described the "curvature" of the mixing; if
  \code{curvature_cj} is 1, this is equivalent to \code{Wcj = min{Wc, Wj}}; if
  \code{curvature_cj} is 0, then \code{Wcj} is the geometric mean of \code{Wc}
  and \code{Wj}. Any value below 1 will produce a smooth curve rather than the
  abrupt transitions that occur when choosing the simple minimum. Likewise, the
  carboxylation rate co-limited by \code{Wcj} and \code{Wp} (called \code{Wcjp};
  in other words, the rate co-limited by all three processes) is given by the
  smaller root of

  \code{curvature_cjp * Wcjp^2 - Wcjp * (Wcj + Wp) + Wcj * Wp = 0}.

  Then, \code{Wcjp} is used in place of \code{Vc} when calculating the net CO2
  assimilation rate. For more information about the quadratic mixing technique,
  see Collatz et al. (1990) and Collatz et al. (1991).

  In the equations above, several of the variables depend on the leaf
  temperature. In particular, the leaf-temperature-adjusted values of
  \code{Vcmax}, \code{J}, and \code{Rd} are determined from their base values
  at 25 degrees C and a temperature-dependent multiplicative factor.

  Also note that \code{PCc} is calculated from the chloroplastic CO2
  concentration \code{Cc} using the total pressure (ambient pressure + chamber
  overpressure).

  In addition to the carboxylation and assimilation rates already mentioned, it
  is also possible to calculate the net CO2 assimilation rates determined by
  Rubisco activity, RuBP regeneration, and TPU as follows:

  \code{Ac = (1 - Gamma_star / PCc) * Wc - Rd}

  \code{Aj = (1 - Gamma_star / PCc) * Wj - Rd}

  \code{Ap = (1 - Gamma_star / PCc) * Wp - Rd}

  References:
  \itemize{
    \item von Caemmerer, S. "Biochemical Models of Leaf Photosynthesis" (CSIRO
          Publishing, 2000) [\doi{10.1071/9780643103405}].
    \item Collatz, G. J., Ball, J. T., Grivet, C. & Berry, J. A. "Physiological
          and environmental regulation of stomatal conductance, photosynthesis
          and transpiration: a model that includes a laminar boundary layer."
          Agricultural and Forest Meteorology 54, 107–136 (1991)
          [\doi{10.1016/0168-1923(91)90002-8}].
    \item Collatz, G. J., Berry, J. A., Farquhar, G. D. & Pierce, J. "The
          relationship between the Rubisco reaction mechanism and models of
          photosynthesis." Plant, Cell & Environment 13, 219–225 (1990)
          [\doi{10.1111/j.1365-3040.1990.tb01306.x}].
  }
}

\value{
  The return value depends on the value of \code{return_exdf}:
  \itemize{
    \item If \code{return_exdf} is \code{TRUE}, the return value is an
          \code{exdf} object with the following columns, calculated as described
          above: \code{TPU}, \code{Vcmax_tl}, \code{Rd_tl}, \code{J_tl},
          \code{Ac}, \code{Aj}, \code{Ap}, \code{An}, and \code{Vc}. The
          category for each of these new columns is
          \code{calculate_c3_assimilation} to indicate that they were created
          using this function.
    \item If \code{return_exdf} is \code{FALSE}, the return value is a list with
          the following named elements: \code{An}, \code{Ac}, \code{Aj}, and
          \code{Ap}. Each element is a numeric vector.
  }
}

\examples{
# Simulate a C3 A-Cc curve with specified leaf temperature and photosynthetic
# parameters and plot the net assimilation rate along with the different
# enzyme-limited rates
inputs <- exdf(data.frame(
  Cc = seq(1, 601, by = 6),
  Tleaf = 30,
  total_pressure = 1
))

inputs <- document_variables(
  inputs,
  c('', 'Cc',             'micromol mol^(-1)'),
  c('', 'Tleaf',          'degrees C'),
  c('', 'total_pressure', 'bar')
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
