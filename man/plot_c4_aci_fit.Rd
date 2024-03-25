\name{plot_c4_aci_fit}

\alias{plot_c4_aci_fit}

\title{Plot the results of a C4 CO2 response curve fit}

\description{
  Plots the output from \code{\link{fit_c4_aci}}.
}

\usage{
  plot_c4_aci_fit(
    fit_results,
    identifier_column_name,
    x_name,
    plot_operating_point = TRUE,
    a_column_name = 'A',
    ci_column_name = 'Ci',
    pcm_column_name = 'PCm',
    \dots
  )
}

\arguments{
  \item{fit_results}{
    A list of three \code{exdf} objects names \code{fits}, \code{parameters},
    and \code{fits_interpolated}, as calculated by \code{\link{fit_c4_aci}}.
  }

  \item{identifier_column_name}{
    The name of a column in each element of \code{fit_results} whose value can
    be used to identify each response curve within the data set; often, this is
    \code{'curve_identifier'}.
  }

  \item{x_name}{
    The name of the column that should be used for the x-axis in the plot. This
    should refer to either \code{'Ci'} or \code{'Cc'}, and it must be the same
    as \code{ci_column_name} or \code{cc_column_name}.
  }

  \item{plot_operating_point}{
    A logical value indicating whether to plot the operating point.
  }

  \item{a_column_name}{
    The name of the column in \code{replicate_exdf} that contains the net
    assimilation in \code{micromol m^(-2) s^(-1)}; should be the same value that
    was passed to \code{fit_c4_aci}.
  }

  \item{ci_column_name}{
    The name of the column in \code{replicate_exdf} that contains the
    intercellular CO2 concentration in \code{micromol mol^(-1)}; should be the
    same value that was passed to \code{fit_c4_aci}.
  }

  \item{pcm_column_name}{
    The name of the column in \code{exdf_obj} that contains the partial pressure
    of CO2 in the mesophyll, expressed in \code{microbar}.
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
  with \code{FUN = fit_c4_aci}.

  The resulting plot will show curves for the fitted rates \code{An},
  \code{Apr}, \code{Apc}, and \code{Ar}, along with points for the measured
  values of \code{A}, and (optionally) the estimated operating point. The x-axis
  can be set to either \code{Ci} or \code{PCm}.

  Internally, this function uses \code{\link[lattice]{xyplot}} to perform the
  plotting. Optionally, additional arguments can be passed to \code{xyplot}.
  These should typically be limited to things like \code{xlim}, \code{ylim},
  \code{main}, and \code{grid}, since many other \code{xyplot} arguments will be
  set internally (such as \code{xlab}, \code{ylab}, \code{auto}, and others).

  See the help file for \code{\link{fit_c4_aci}} for an example using this
  function.
}


\value{
  A \code{trellis} object created by \code{lattice::xyplot}.
}

\examples{
# See the help file for fit_c4_aci for an example using this function.
}

\concept{exdf}
