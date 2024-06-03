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
    alpha_g,
    alpha_old,
    alpha_s,
    Gamma_star,
    J_at_25,
    Rd_at_25,
    Tp,
    Vcmax_at_25,
    atp_use = 4.0,
    nadph_use = 8.0,
    curvature_cj = 1.0,
    curvature_cjp = 1.0,
    cc_column_name = 'Cc',
    j_norm_column_name = 'J_norm',
    kc_column_name = 'Kc',
    ko_column_name = 'Ko',
    oxygen_column_name = 'oxygen',
    rd_norm_column_name = 'Rd_norm',
    total_pressure_column_name = 'total_pressure',
    vcmax_norm_column_name = 'Vcmax_norm',
    perform_checks = TRUE,
    return_exdf = TRUE
  )
}

\arguments{
  \item{exdf_obj}{An \code{exdf} object.}

  \item{alpha_g}{
    A dimensionless parameter where \code{0 <= alpha_g <= 1}, representing the
    proportion of glycolate carbon taken out of the photorespiratory pathway as
    glycine. \code{alpha_g} is often assumed to be 0. If \code{alpha_g} is not a
    number, then there must be a column in \code{exdf_obj} called \code{alpha_g}
    with appropriate units. A numeric value supplied here will overwrite the
    values in the \code{alpha_g} column of \code{exdf_obj} if it exists.
  }

  \item{alpha_old}{
    A dimensionless parameter where \code{0 <= alpha_old <= 1}, representing the
    fraction of remaining glycolate carbon not returned to the chloroplast
    after accounting for carbon released as CO2. \code{alpha_old} is often
    assumed to be 0. If \code{alpha_old} is not a number, then there must be a
    column in \code{exdf_obj} called \code{alpha_old} with appropriate units.
    A numeric value supplied here will overwrite the values in the
    \code{alpha_old} column of \code{exdf_obj} if it exists.
  }

  \item{alpha_s}{
    A dimensionless parameter where \code{0 <= alpha_s <= 0.75 * (1 - alpha_g)}
    representing the proportion of glycolate carbon taken out of the
    photorespiratory pathway as serine. \code{alpha_s} is often assumed to be 0.
    If \code{alpha_s} is not a number, then there must be a column in
    \code{exdf_obj} called \code{alpha_s} with appropriate units. A numeric
    value supplied here will overwrite the values in the \code{alpha_s} column
    of \code{exdf_obj} if it exists.
  }

  \item{Gamma_star}{
    The CO2 compensation point in the absence of day respiration, expressed in
    \code{micromol mol^(-1)}. If \code{Gamma_star} is not a number, then there
    must be a column in \code{exdf_obj} called \code{Gamma_star} with
    appropriate units. A numeric value supplied here will overwrite the values
    in the \code{Gamma_star} column of \code{exdf_obj} if it exists.
  }

  \item{J_at_25}{
    The electron transport rate at 25 degrees C, expressed in
    \code{micromol m^(-2) s^(-1)}. Note that this is _not_ \code{Jmax}, and in
    general will depend on the incident photosynthetically active flux density.
    If \code{J_at_25} is not a number, then there must be a column in
    \code{exdf_obj} called \code{J_at_25} with appropriate units. A numeric
    value supplied here will overwrite the values in the \code{J_at_25} column
    of \code{exdf_obj} if it exists.
  }

  \item{Rd_at_25}{
    The respiration rate at 25 degrees C, expressed in
    \code{micromol m^(-2) s^(-1)}. If \code{Rd_at_25} is not a number, then there
    must be a column in \code{exdf_obj} called \code{Rd_at_25} with appropriate
    units. A numeric value supplied here will overwrite the values in the
    \code{Rd_at_25} column of \code{exdf_obj} if it exists.
  }

  \item{Tp}{
    The maximum rate of triphosphate utilization, expressed in
    \code{micromol m^(-2) s^(-1)}. If \code{Tp} is not a number, then there
    must be a column in \code{exdf_obj} called \code{Tp} with appropriate
    units. A numeric value supplied here will overwrite the values in the
    \code{Tp} column of \code{exdf_obj} if it exists.
  }

  \item{Vcmax_at_25}{
    The maximum rate of rubisco carboxylation at 25 degrees C, expressed in
    \code{micromol m^(-2) s^(-1)}. If \code{Vcmax_at_25} is not a number, then
    there must be a column in \code{exdf_obj} called \code{Vcmax_at_25} with
    appropriate units. A numeric value supplied here will overwrite the
    values in the \code{Vcmax_at_25} column of \code{exdf_obj} if it exists.
  }

  \item{atp_use}{
    The number of ATP molecules used per C3 cycle.
  }

  \item{nadph_use}{
    The number of NADPH molecules used per C3 cycle.
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

  \item{oxygen_column_name}{
    The name of the column in \code{exdf_obj} that contains the concentration
    of O2 in the ambient air, expressed as a percentage (commonly 21\% or 2\%);
    the units must be \code{percent}.
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
  \strong{The Busch et al. (2018) model:}

  This function generally follows the Farquhar-von-Caemmerer-Berry model as
  described in Busch et al. (2018) with a few modifications described below. In
  this formulation, the steady-state net CO2 assimilation rate \code{An} is
  calculated according to

  \code{An = (1 - Gamma_star_ag / PCc) * Vc - Rd},

  where \code{Gamma_star} is the CO2 compensation point in the absence of day
  respiration, \code{Gamma_star_ag} is the effective value of \code{Gamma_star}
  accounting for glycolate carbon remaining in the cytosol, \code{PCc} is the
  partial pressure of CO2 in the chloroplast, \code{Vc} is the RuBP
  carboxylation rate, and \code{Rd} is the rate of day respiration.
  \code{Gamma_star_ag} is given by

  \code{Gamma_star_ag = (1 - alpha_g) * Gamma_star},

  where \code{alpha_g} is the fraction of glycolate carbon leaving the
  photorespiratory pathway as glycine.

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

  \code{Wj = PCc * J / (4 * PCc + Gamma_star_ag * (8 + 16 * alpha_g + 8 * alpha_s))},

  where \code{J} is the RuBP regeneration rate and \code{alpha_s} is the
  fraction of glycolate carbon leaving the photorespiratory pathway as serine.

  The TPU-limited carboxylation rate is given by

  \code{Wp = PCc * 3 * Tp / (PCc - Gamma_star_ag * (1 + 3 * alpha_g + 4 * alpha_s))},

  where \code{Tp} is the maximum rate of triose phosphate utilization. Note that
  this equation only applies when \code{PCc > Gamma_star_ag * (1 + 3 * alpha_g +
  4 * alpha_s)}; for smaller values of \code{PCc}, TPU cannot limit the RuBP
  carboxylation rate and \code{Wp = Inf}. (Lochocki & McGrath, submitted).

  The actual carboxylation rate is typically chosen to be the smallest of the
  three potential rates:

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
  Rubisco activity, RuBP regeneration, and Tp as follows:

  \code{Ac = (1 - Gamma_star_ag / PCc) * Wc - Rd}

  \code{Aj = (1 - Gamma_star_ag / PCc) * Wj - Rd}

  \code{Ap = (1 - Gamma_star_ag / PCc) * Wp - Rd}

  \strong{The Busch model with nitrogen restrictions}:

  Note that the implementation as described above does not currently facilitate
  the inclusion of nitrogen limitations (Equations 15-21 in Busch et al.
  (2018)).

  \strong{The "old" model:}

  In an older version of the model, \code{alpha_g} and \code{alpha_s} are
  replaced with a single parameter \code{alpha_old}. Most publications refer to
  this simply as \code{alpha}, but here we follow the notation of Busch et al.
  (2018) for clarity. In this version, there is no disctinction between
  \code{Gamma_star_ag} and \code{Gamma_star}. Other differences are as follows.

  he RuBP-regeneration-limited carboxylation rate \code{Wj} is given by

  \code{Wj = PCc * J / (atp_use * PCc + nadph_use * Gamma_star)},

  Here we have allowed \code{atp_use} and \code{nadph_use} to be variables
  rather than taking fixed values (as they do in many sources). This is
  necessary because not all descriptions of the FvCB model use the same values,
  where the different values are due to different assumptions about the energy
  requirements of RuBP regeneration.

  The TPU-limited carboxylation rate is given by

  \code{Wp = PCc * 3 * Tp / (PCc - Gamma_star * (1 + 3 * alpha_old))},

  Note that this equation only applies when \code{PCc > Gamma_star * (1 +
  3 * alpha_old)}; for smaller values of \code{PCc}, TPU cannot limit the RuBP
  carboxylation rate and \code{Wp = Inf}. (Lochocki & McGrath, submitted).

  \strong{Using either version of the model:}

  When using \code{calculate_c3_assimilation}, it is possible to use either
  version of the model. Setting \code{alpha_g} and \code{alpha_s} to zero is
  equivalent to using the older version of the model, while setting
  \code{alpha_old = 0} is equivalent to using the newer version of the model.
  If all \code{alpha} parameters are zero, there is effectively no difference
  between the two versions of the model. Attempting to set a nonzero
  \code{alpha_old} if either \code{alpha_g} or \code{alpha_s} is nonzero is
  forbidden since it would represent a mix between the two models; if such
  values are passed as inputs, then an error will be thrown.

  \strong{References:}

  \itemize{
    \item Busch, Sage, & Farquhar, G. D. "Plants increase CO2 uptake by
          assimilating nitrogen via the photorespiratory pathway." Nature Plants
          4, 46–54 (2018) [\doi{10.1038/s41477-017-0065-x}].
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
    \item Lochocki & McGrath "Widely Used Variants of the
          Farquhar-von-Caemmerer-Berry Model Can Cause Errors in Parameter
          Estimates and Simulations." submitted.
  }
}

\value{
  The return value depends on the value of \code{return_exdf}:
  \itemize{
    \item If \code{return_exdf} is \code{TRUE}, the return value is an
          \code{exdf} object with the following columns, calculated as described
          above: \code{alpha_g}, \code{Gamma_star}, \code{Tp}, \code{Vcmax_tl},
          \code{Rd_tl}, \code{J_tl}, \code{Ac}, \code{Aj}, \code{Ap}, \code{An},
          and \code{Vc}. The category for each of these new columns is
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
  total_pressure = 1,
  oxygen = 21
))

inputs <- document_variables(
  inputs,
  c('', 'Cc',             'micromol mol^(-1)'),
  c('', 'Tleaf',          'degrees C'),
  c('', 'total_pressure', 'bar'),
  c('', 'oxygen',         'percent')
)

inputs <- calculate_arrhenius(inputs, c3_arrhenius_sharkey, 'Tleaf')

assim <- calculate_c3_assimilation(inputs, 0, 0, 0, '', 150, 1, 12, 120)

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
