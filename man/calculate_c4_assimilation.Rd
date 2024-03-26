\name{calculate_c4_assimilation}

\alias{calculate_c4_assimilation}

\title{Calculate C4 assimilation rates}

\description{
  Calculates C4 assimilation rates based on the von Caemmerer (2000) model. This
  function can accomodate alternative colum names for the variables taken from
  Licor files in case they change at some point in the future. This function
  also checks the units of each required column and will produce an error if any
  units are incorrect.
}

\usage{
  calculate_c4_assimilation(
    exdf_obj,
    alpha_psii,
    gbs,
    Rd_at_25,
    Rm_frac,
    Vcmax_at_25,
    Vpmax_at_25,
    Vpr,
    ao_column_name = 'ao',
    gamma_star_column_name = 'gamma_star',
    kc_column_name = 'Kc',
    ko_column_name = 'Ko',
    kp_column_name = 'Kp',
    oxygen_column_name = 'oxygen',
    pcm_column_name = 'PCm',
    rd_norm_column_name = 'Rd_norm',
    total_pressure_column_name = 'total_pressure',
    vcmax_norm_column_name = 'Vcmax_norm',
    vpmax_norm_column_name = 'Vpmax_norm',
    perform_checks = TRUE,
    return_exdf = TRUE
  )
}

\arguments{
  \item{exdf_obj}{An \code{exdf} object.}

  \item{alpha_psii}{
    The fraction of photosystem II activity in the bundle sheath
    (\code{dimensionless}). If \code{alpha_psii} is not a number, then there
    must be a column in \code{exdf_obj} called \code{alpha_psii} with
    appropriate units. A numeric value supplied here will overwrite the values
    in the \code{alpha_psii} column of \code{exdf_obj} if it exists.
  }

  \item{gbs}{
    The bundle sheath conductance to CO2 in \code{mol m^(-2) s^(-1) bar^(-1)}.
    If \code{gbs} is not a number, then there must be a column in
    \code{exdf_obj} called \code{gbs} with appropriate units. A numeric value
    supplied here will overwrite the values in the \code{gbs} column of
    \code{exdf_obj} if it exists.
  }

  \item{Rd_at_25}{
    The total rate of mitochondrial respiration across the mesophyll and bundle
    sheath at 25 degrees C, expressed in \code{micromol m^(-2) s^(-1)}. If
    \code{Rd_at_25} is not a number, then there must be a column in
    \code{exdf_obj} called \code{Rd_at_25} with appropriate units. A numeric
    value supplied here will overwrite the values in the \code{Rd_at_25} column
    of \code{exdf_obj} if it exists.
  }

  \item{Rm_frac}{
    The fraction of the total mitochondrial respiration that occurs in the
    mesophyll. If \code{Rm_frac} is not a number, then there must be a column in
    \code{exdf_obj} called \code{Rm_frac} with appropriate units. A numeric
    value supplied here will overwrite the values in the \code{Rm_frac} column
    of \code{exdf_obj} if it exists.
  }

  \item{Vcmax_at_25}{
    The maximum rate of rubisco carboxylation at 25 degrees C, expressed in
    \code{micromol m^(-2) s^(-1)}. If \code{Vcmax_at_25} is not a number, then
    there must be a column in \code{exdf_obj} called \code{Vcmax_at_25} with
    appropriate units. A numeric value supplied here will overwrite the values
    in the \code{Vcmax_at_25} column of \code{exdf_obj} if it exists.
  }

  \item{Vpmax_at_25}{
    The maximum rate of PEP carboxylase activity at 25 degrees C, expressed in
    \code{micromol m^(-2) s^(-1)}. If \code{Vpmax_at_25} is not a number, then
    there must be a column in \code{exdf_obj} called \code{Vpmax_at_25} with
    appropriate units. A numeric value supplied here will overwrite the values
    in the \code{Vpmax_at_25} column of \code{exdf_obj} if it exists.
  }

  \item{Vpr}{
    The rate of PEP carboxylase regeneration, expressed in
    \code{micromol m^(-2) s^(-1)}. If \code{Vpr} is not a number, then there
    must be a column in \code{exdf_obj} called \code{Vpr} with appropriate
    units. A numeric value supplied here will overwrite the values in the
    \code{Vpr} column of \code{exdf_obj} if it exists.
  }

  \item{ao_column_name}{
    The name of the column in \code{exdf_obj} that contains the dimensionless
    ratio of solubility and diffusivity of O2 to CO2.
  }

  \item{gamma_star_column_name}{
    The name of the column in \code{exdf_obj} that contains the dimensionless
    \code{gamma_star} values.
  }

  \item{kc_column_name}{
    The name of the column in \code{exdf_obj} that contains the Michaelis-Menten
    constant for rubisco carboxylation in \code{microbar}.
  }

  \item{ko_column_name}{
    The name of the column in \code{exdf_obj} that contains the Michaelis-Menten
    constant for rubisco oxygenation in \code{mbar}.
  }

  \item{kp_column_name}{
    The name of the column in \code{exdf_obj} that contains the Michaelis-Menten
    constant for PEP carboxylase carboxylation in \code{microbar}.
  }

  \item{oxygen_column_name}{
    The name of the column in \code{exdf_obj} that contains the concentration
    of O2 in the ambient air, expressed as a percentage (commonly 21\% or 2\%);
    the units must be \code{percent}.
  }

  \item{pcm_column_name}{
    The name of the column in \code{exdf_obj} that contains the partial pressure
    of CO2 in the mesophyll, expressed in \code{microbar}.
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

  \item{vpmax_norm_column_name}{
    The name of the column in \code{exdf_obj} that contains the normalized
    \code{Vpmax} values (with units of
    \code{normalized to Vpmax at 25 degrees C}).
  }

  \item{perform_checks}{
    A logical value indicating whether to check units for the required columns.
    This should almost always be \code{TRUE}. The option to disable these checks
    is only intended to be used when \code{\link{fit_c4_aci}} calls this
    function, since performing these checks many times repeatedly slows down the
    fitting procedure.
  }

  \item{return_exdf}{
    A logical value indicating whether to return an \code{exdf} object. This
    should almost always be \code{TRUE}. The option to return a vector is mainly
    intended to be used when  \code{\link{fit_c4_aci}} calls this function,
    since creating an \code{exdf} object to return will slow down the fitting
    procedure.
  }
}

\details{
  \strong{General Description of the Model}

  This function generally follows Section 4.2.1 from S. von Caemmerer (2000),
  which provides equations for calculating the enzyme-limited net assimilation
  rate \code{An} in a C4 leaf. (These equations are also reproduced in S. von
  Caemmerer (2021), although we use the equation numbers from the 2000 textbook
  here.) The net assimilation rate in this model is calculated according to
  Equation 4.21:

  \code{An = (-b - sqrt(b^2 - 4 * a * c)) / (2 * a)}

  where the parameters \code{a}, \code{b}, and \code{c} are determined by
  Equations 4.22, 4.23, and 4.24, respectively. These equations are fairly long,
  so we do not reproduce them here.

  \strong{An Approximation to the Full Equations}

  The complicated equations above can be approximiated by a simpler one
  (Equation 4.25):

  \code{An = min(Vp + gbs * PCm - Rm, Vcmax - Rd)}

  where \code{Vp} is the rate of PEP carboxylation, \code{gbs} is the bundle
  sheath conductance to CO2, \code{PCm} is the partial pressure of CO2 in the
  mesophyll, \code{Rm} is the rate of mitochondrial respiration occuring in the
  mesophyll, \code{Vcmax} is the maximum rate of Rubisco carboxylation, and
  \code{Rd} is the rate of mitochondrial respiration occurring in the bundle
  sheath and mesophyll. Essentially, the first term in the equation above
  (\code{Vp + gbs * PCm - Rm}) can be thought of as a PEP-carboxylase-limited
  assimilation rate \code{Ap}, while the second term (\code{Vcmax - Rd}) is a
  Rubisco-limited rate \code{Ar}.

  The PEP carboxylation rate \code{Vp} is calculated according to Equation 4.19:

  \code{Vp = min(Pcm * Vpmax / (PCm + Kp), Vpr)}

  where \code{Vpmax} is the maximum rate of PEP carboxylation, \code{Kp} is a
  Michaelis-Menten constant for PEP carboxylation, and \code{Vpr} is the
  carboxylation rate when PEP carboxylase activity is limited by regeneration
  rather than carbon availability. Thus, we can see that the approximation above
  actually calculates the net assimilation rate as the smaller of three separate
  assimilation rates:

  \code{An = min(Apc, Apr, Ar)}

  where \code{Apc = Pcm * Vpmax / (PCm + Kp) + gbs * PCm - Rm} is the rate due
  to carbon-limited PEP carboxylation, \code{Apr = Vpr + gbs * PCm - Rm} is the
  rate due to regeneration-limited PEP carboxylation, and \code{Ar = Vcmax - Rd}
  is the rate due to Rubisco-limited assimilation.

  In the example at the end of this documentation page, we compare \code{Apc},
  \code{Apr}, and \code{Ar} to \code{An} as calculated by Equation 4.21. From
  this example, it is clear that the approximation \code{An = min(Apc, Apr, Ar)}
  is quite accurate for low values of \code{PCm}, but introduces significant
  errors as \code{PCm} increases. Thus, while the approximation can be helpful
  for gaining an intuitive understanding of C4 photosynthesis, it should not be
  used for realistic calculations.

  To be more precise, the approximation is only reliable when \code{Vcmax} is
  much larger than \code{gbs * Kc * (1 + POm / Ko)}, which is rarely the case;
  otherwise, the limiting value of \code{An} at high \code{PCm} will be smaller
  than \code{Ar = Vcmax - Rd}. Conversely, if \code{gbs} and \code{alpha_psii} are
  both set to zero, then the approximation is exact.

  \strong{Limiting Cases of the Approximate Equation}

  The bundle sheath conductance \code{gbs} is generally very small and can be
  ignored in a simple analysis of the above equations. In that case, when
  \code{Pcm} is very high, the approximate equation for \code{An} simplifies
  further to:

  \code{An = min(Vpmax - Rm, Vpr - Rm, Vcmax - Rd)}

  Since respiration costs are also generally much smaller than the maximum
  enzyme activity and regeneration rates, the assimilation rate at high levels
  of CO2 is therefore determined by the smaller of \code{Vpmax}, \code{Vpr}, and
  \code{Vcmax}. As shown in Table 4.1 of the textbook, \code{Vpmax} is typically
  much larger than the other two rates, so high CO2 assimilation in C4 leaves is
  usually limited by either \code{Vpr} or \code{Vcmax}. The exact limiting
  factor can depend on many possible variables, such as the temperature. For
  example, see Wang (2008).

  At lower values of \code{PCm}, net assimilation is determined by CO2-limited
  PEP carboxylation according to:

  \code{An = PCm * Vpmax / Kp - Rm}

  where we have approximated \code{gbs * PCm = 0} and \code{PCm + Kp = Kp}, as
  appropriate for small values of \code{Pcm}. Thus, we can see that for low CO2
  levels, assimilation is linearly related to \code{PCm} with a slope of
  \code{Vpmax / Kp} and intercept of \code{-Rm}.

  \strong{Respiration}

  Table 4.1 from von Caemmerer (2000) suggests that \code{Rd = 0.01 * Vcmax} and
  \code{Rm = 0.5 * Rd}. To allow more flexibility, we allow \code{Rd} to be
  specified independently of \code{Vcmax}, and we also consider the ratio of
  \code{Rm / Rd = Rm_frac} to be a variable (so that \code{Rm} is calculated
  from \code{Rd} according to \code{Rm = Rm_frac * Rd}). If \code{Rm_frac} is
  set to 1, then there is no distinction between \code{Rd} and \code{Rm}.

  \strong{References}

  \itemize{
    \item von Caemmerer, S. "Biochemical Models of Leaf Photosynthesis" (CSIRO
          Publishing, 2000) [\doi{10.1071/9780643103405}].
    \item von Caemmerer, S. "Updating the steady-state model of C4
          photosynthesis." Journal of Experimental Botany 72, 6003–6017 (2021)
          [\doi{10.1093/jxb/erab266}].
    \item Wang, D., Portis, A. R., Jr., Moose, S. P. & Long, S. P. "Cool C4
          Photosynthesis: Pyruvate Pi Dikinase Expression and Activity
          Corresponds to the Exceptional Cold Tolerance of Carbon Assimilation
          in Miscanthus × giganteus." Plant Physiology 148, 557–567 (2008)
          [\doi{10.1104/pp.108.120709}].
  }
}

\value{
  The return value depends on the value of \code{return_exdf}:
  \itemize{
    \item If \code{return_exdf} is \code{TRUE}, the return value is an
          \code{exdf} object with the following columns: \code{alpha_psii},
          \code{gbs}, \code{Rm_frac}, \code{Vcmax_tl}, \code{Vpmax_tl},
          \code{Rd_tl}, \code{Rm_tl}, \code{Vpc}, \code{Vpr}, \code{Vp},
          \code{Apc}, \code{Apr}, \code{Ap}, \code{Ar}, \code{An}, and
          \code{c4_assimilation_msg}. Most of these are calculated as described
          above, while several are copies of the input arguments with the same
          name. The \code{c4_assimilation_msg} is usually blank but may contain
          information about any issues with the inputs. The category for each of
          these new columns is \code{calculate_c4_assimilation} to indicate that
          they were created using this function.
    \item If \code{return_exdf} is \code{FALSE}, the return value is a numeric
          vector containing the calculated values of \code{An}.
  }
}

\examples{
# Simulate a C4 A-Cm curve with specified leaf temperature and photosynthetic
# parameters and plot the net assimilation rate.
npts <- 101

inputs <- exdf(data.frame(
  PCm = seq(0, 500, length.out = npts),
  Tleaf = 25,
  total_pressure = 1,
  oxygen = 21
))

inputs <- document_variables(
  inputs,
  c('', 'PCm',            'microbar'),
  c('', 'Tleaf',          'degrees C'),
  c('', 'total_pressure', 'bar'),
  c('', 'oxygen',         'percent')
)

inputs <- calculate_arrhenius(inputs, c4_arrhenius_von_caemmerer, 'Tleaf')

assim <- calculate_c4_assimilation(inputs, 0, 0.003, 1, 0.5, 40, 200, 80)

# Now we can plot An, Apr, Apc, and Ar. From this plot, we can see that
# replacing the complicated quadratic equation with a simple minimum yields
# very different results. Although this approximation is helpful for
# understanding C4 photosythesis, it should not be used for calculations.
lattice::xyplot(
  Apr + Apc + Ar + An ~ PCm,
  data = cbind(inputs, assim)$main_data,
  type = 'l',
  grid = TRUE,
  auto = TRUE,
  xlab = paste0('Partial pressure of CO2 in the mesophyll (', inputs$units$PCm, ')'),
  ylab = paste0('Net CO2 assimilation rate (', assim$units$An, ')')
)
}

\concept{exdf}
