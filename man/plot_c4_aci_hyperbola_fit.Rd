\name{plot_c4_aci_hyperbola_fit}

\alias{plot_c4_aci_hyperbola_fit}

\title{Plot the results of a hyperbolic C4 CO2 response curve fit}

\description{
  Plots the output from \code{\link{fit_c4_aci_hyperbola}}.
}

\usage{
  plot_c4_aci_hyperbola_fit(
    fit_results,
    identifier_column_name,
    a_column_name = 'A',
    ci_column_name = 'Ci',
    \dots
  )
}

\arguments{
  \item{fit_results}{
    A list of three \code{exdf} objects named \code{fits}, \code{parameters},
    and \code{fits_interpolated}, as calculated by
    \code{\link{fit_c4_aci_hyperbola}}.
  }

  \item{identifier_column_name}{
    The name of a column in each element of \code{fit_results} whose value can
    be used to identify each response curve within the data set; often, this is
    \code{'curve_identifier'}.
  }

  \item{a_column_name}{
    The name of the columns in the elements of \code{fit_results} that contain
    the net assimilation in \code{micromol m^(-2) s^(-1)}; should be the same
    value that was passed to \code{fit_c4_aci_hyperbola}.
  }

  \item{ci_column_name}{
    The name of the columns in the elements of \code{fit_results} that contain
    the intercellular CO2 concentration in \code{micromol mol^(-1)}; should be
    the same value that was passed to \code{fit_c4_aci_hyperbola}.
  }

  \item{\dots}{
    Additional arguments to be passed to \code{\link[lattice]{xyplot}}.
  }
}

\details{
  This is a convenience function for plotting the results of a C4 A-Ci curve
  fit. It is typically used for displaying several fits at once, in which case
  \code{fit_results} is actually the output from calling
  \code{\link{consolidate}} on a list created by calling \code{\link{by.exdf}}
  with \code{FUN = fit_c4_aci_hyperbola}.

  The resulting plot will show curves for the fitted rates \code{An},
  \code{Ainitial}, and \code{Amax}, along with points for the measured values of
  \code{A}.

  Internally, this function uses \code{\link[lattice]{xyplot}} to perform the
  plotting. Optionally, additional arguments can be passed to \code{xyplot}.
  These should typically be limited to things like \code{xlim}, \code{ylim},
  \code{main}, and \code{grid}, since many other \code{xyplot} arguments will be
  set internally (such as \code{xlab}, \code{ylab}, \code{auto}, and others).

  See the help file for \code{\link{fit_c4_aci_hyperbola}} for an example using
  this function.
}


\value{
  A \code{trellis} object created by \code{lattice::xyplot}.
}

\examples{
# See the help file for fit_c4_aci_hyperbola for an example using this function.
}

\concept{exdf}
