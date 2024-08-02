\name{calculate_c4_assimilation_hyperbola}

\alias{calculate_c4_assimilation_hyperbola}

\title{Calculate C4 assimilation rates using a hyperbola}

\description{
  Calculates C4 assimilation rates based on an empirical hyperbolic model. This
  function can accomodate alternative colum names for the variables taken from
  Licor files in case they change at some point in the future. This function
  also checks the units of each required column and will produce an error if any
  units are incorrect.
}

\usage{
  calculate_c4_assimilation_hyperbola(
    exdf_obj,
    c4_curvature,
    c4_slope,
    rL,
    Vmax,
    ci_column_name = 'Ci',
    hard_constraints = 0,
    perform_checks = TRUE,
    return_exdf = TRUE
  )
}

\arguments{
  \item{exdf_obj}{An \code{exdf} object.}

  \item{c4_curvature}{
    The empirical curvature parameter of the hyperbola (\code{dimensionless}).
    If \code{c4_curvature} is not a number, then there must be a column in
    \code{exdf_obj} called \code{c4_curvature} with appropriate units. A numeric
    value supplied here will overwrite the values in the \code{c4_curvature}
    column of \code{exdf_obj} if it exists.
  }

  \item{c4_slope}{
    The empirical slope parameter of the hyperbola (\code{mol m^(-2) s^(-1)}).
    If \code{c4_slope} is not a number, then there must be a column in
    \code{exdf_obj} called \code{c4_slope} with appropriate units. A numeric
    value supplied here will overwrite the values in the \code{c4_slope}
    column of \code{exdf_obj} if it exists.
  }

  \item{rL}{
    The respiration rate, expressed in \code{micromol m^(-2) s^(-1)}. If
    \code{rL} is not a number, then there must be a column in \code{exdf_obj}
    called \code{rL} with appropriate units. A numeric value supplied here will
    overwrite the values in the \code{rL} column of \code{exdf_obj} if it
    exists.
  }

  \item{Vmax}{
    The maximum gross assimilation rate, expressed in
    \code{micromol m^(-2) s^(-1)}. If \code{Vmax} is not a number, then there
    must be a column in \code{exdf_obj} called \code{Vmax} with appropriate
    units. A numeric value supplied here will overwrite the values in the
    \code{Vmax} column of \code{exdf_obj} if it exists.
  }

  \item{ci_column_name}{
    The name of the column in \code{exdf_obj} that contains the intercellular
    CO2 concentration, expressed in \code{micromol mol^(-1)}.
  }

  \item{hard_constraints}{
    An integer numerical value indicating which types of hard constraints to
    place on the values of input parameters; see below for more details.
  }

  \item{perform_checks}{
    A logical value indicating whether to check units for the required columns.
    This should almost always be \code{TRUE}. The option to disable these checks
    is only intended to be used when \code{\link{fit_c4_aci_hyperbola}} calls
    this function, since performing these checks many times repeatedly slows
    down the fitting procedure.
  }

  \item{return_exdf}{
    A logical value indicating whether to return an \code{exdf} object. This
    should almost always be \code{TRUE}. The option to return a vector is mainly
    intended to be used when  \code{\link{fit_c4_aci_hyperbola}} calls this
    function, since creating an \code{exdf} object to return will slow down the
    fitting procedure.
  }
}

\details{
  \strong{General Description of the Model}

  In contrast to the mechanistic model implemented in
  \code{\link{calculate_c4_assimilation}}, this is a simple empirical model for
  C4 assimilation based on a four-parameter hyperbola. In this model, the net
  CO2 assimilation rate (\code{An}) is given by

  \code{An = Ag - rL},

  where \code{Ag} is the gross assimilation rate and \code{rL} is the
  respiration rate. In turn, \code{Ag} is given by the smaller root of the
  following quadratic equation:

  \code{curvature * Ag^2 - (Vinitial + Vmax) * Ag + Vinitial * Vmax = 0},

  where 0 <= \code{curvature} <= 1 is an empirical curvature factor, \code{Vmax}
  is the maximum gross assimilation rate, and \code{Vinitial} represents the
  initial response of \code{Ag} to increases in the intercellular CO2
  concentration (\code{Ci}):

  \code{Vinitial = slope * Ci}.

  Here the \code{slope} is another empirical factor.

  By including the respiration offset, it is also possible to define two other
  quantities: the maximum net CO2 assimilation rate (\code{Amax}) and the
  initial net CO2 assimilation rate (\code{Ainitial}). These are given by

  \code{Amax = Vmax - rL}

  and

  \code{Ainitial = Vinitial - rL}.

  Overall, this model exhibits a linear response of \code{An} to \code{Ci} at
  low \code{Ci}, a flat plateau of \code{An} at high \code{Ci}, and a smooth
  transition between these regions. The sharpess of the transition is set by the
  \code{curvature}. When \code{curvature = 1}, the model simplifies to

  \code{An = min{Vinitial, Vmax} - rL = min{Ainitial, Amax}}.

  As the \code{curvature} increases to 1, the transition becomes smoother. When
  the \code{curvature} is not zero, \code{An} approaches \code{Amax}
  asymptotically, and may not reach \code{Amax} at a reasonable value of
  \code{Ci}.

  \strong{Code implementation}

  In this function, \code{curvature} and \code{slope} above are referred to as
  \code{c4_curvature} and \code{c4_slope} to avoid any potential ambiguity with
  other models that may also have curvature and slope parameters.

  \strong{Temperature response}

  Because this model does not represent any photosynthetic mechanisms,
  temperature response functions are not applied.

  \strong{Hard constraints}

  Most input parameters to the this model have hard constraints on their values
  which are set by their interpretation; for example, \code{Vmax} cannot be
  negative and \code{c4_curvature} must lie between 0 and 1. Yet, because of
  measurement noise, sometimes it is necessary to use values outside these
  ranges when fitting an A-Ci curve with \code{\link{fit_c4_aci_hyperbola}}. To
  accomodate different potential use cases, it is possible to selectively apply
  these hard constraints by specifying different values of the
  \code{hard_constraints} input argument:

  \itemize{
    \item \code{hard_constraints = 0}: No constraints are applied.

    \item \code{hard_constraints = 1}: Checks whether all \code{Ci} values are
          non-negative.

    \item \code{hard_constraints = 2}: Includes the same constraints as when
          \code{hard_constraints} is 1, which additional constraints on the
          parameters that can be fitted. For example, \code{Vmax} must be
          non-negative and \code{c4_curvature} must lie between 0 and 1.
  }

  If any input values violate any of the specified constraints, an error message
  will be thrown.
}

\value{
  The return value depends on the value of \code{return_exdf}:
  \itemize{
    \item If \code{return_exdf} is \code{TRUE}, the return value is an
          \code{exdf} object with the following columns: \code{Ag},
          \code{Ainitial}, \code{Amax}, \code{An}, \code{c4_curvature},
          \code{c4_slope}, \code{rL}, \code{Vinitial}, \code{Vmax}, and
          \code{c4_assimilation_hyperbola_msg}. Most of these are calculated as
          described above, while several are copies of the input arguments with
          the same name. The \code{c4_assimilation_hyperbola_msg} is usually
          blank but may contain information about any issues with the inputs.
          The category for each of these new columns is
          \code{calculate_c4_assimilation_hyperbola} to indicate that they were
          created using this function.
    \item If \code{return_exdf} is \code{FALSE}, the return value is a numeric
          vector containing the calculated values of \code{An}.
  }
}

\examples{
# Simulate a C4 A-Ci curve and plot the net assimilation rate.
npts <- 101

inputs <- exdf(data.frame(
  Ci = seq(0, 1000, length.out = npts),
  total_pressure = 1
))

inputs <- document_variables(
  inputs,
  c('', 'Ci',             'micromol mol^(-1)'),
  c('', 'total_pressure', 'bar')
)

assim <- calculate_c4_assimilation_hyperbola(inputs, 0.8, 0.5, 1.0, 55)

lattice::xyplot(
  Ainitial + Amax + An ~ Ci,
  data = cbind(inputs, assim)$main_data,
  type = 'l',
  grid = TRUE,
  auto = TRUE,
  ylim = c(-5, 65),
  xlab = paste0('Intercellular CO2 concentration (', inputs$units$Ci, ')'),
  ylab = paste0('Net CO2 assimilation rate (', assim$units$An, ')')
)
}

\concept{exdf}
