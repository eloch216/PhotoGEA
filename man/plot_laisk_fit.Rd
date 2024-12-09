\name{plot_laisk_fit}

\alias{plot_laisk_fit}

\title{Plot the results of a C3 CO2 response curve fit}

\description{
  Plots the output from \code{\link{fit_laisk}}.
}

\usage{
  plot_laisk_fit(
    fit_results,
    identifier_column_name,
    plot_type,
    cols = multi_curve_colors(),
    a_column_name = 'A',
    ci_column_name = 'Ci',
    ppfd_column_name = 'PPFD',
    \dots
  )
}

\arguments{
  \item{fit_results}{
    A list of four \code{exdf} objects named \code{first_fits},
    \code{first_fit_parameters}, \code{second_fits},
    and \code{second_fit_parameters}, as calculated by
    \code{\link{fit_laisk}}.
  }

  \item{identifier_column_name}{
    The name of a column in each element of \code{fit_results} whose value can
    be used to identify each replicate within the data set; often, this is
    \code{'curve_identifier'}.
  }

  \item{plot_type}{
    Must be either \code{'first'} or \code{'second'} (case-insensitive);
    determines which type of plot to create (see below for details).
  }

  \item{cols}{
    A vector of color specifications to use for each light level when plotting.
  }

  \item{a_column_name}{
    The name of the columns in the elements of \code{fit_results} that contain
    the net assimilation in \code{micromol m^(-2) s^(-1)}; should be the same
    value that was passed to \code{fit_laisk}.
  }

  \item{ci_column_name}{
    The name of the column in the elements of \code{fit_results} that contain
    the intercellular CO2 concentration in \code{micromol mol^(-1)}; should be
    the same value that was passed to \code{fit_laisk}.
  }

  \item{ppfd_column_name}{
    The name of the column in the elements of \code{fit_results} that can be
    used to split the data into individual response curves; should be the same
    value that was passed to \code{fit_laisk}.
  }

  \item{\dots}{
    Additional arguments to be passed to \code{\link[lattice]{xyplot}}.
  }
}

\details{
  This is a convenience function for plotting the results of a Laisk curve
  fit. It is typically used for displaying several fits at once, in which case
  \code{fit_results} is actually the output from calling
  \code{\link{consolidate}} on a list created by calling \code{\link{by.exdf}}
  with \code{FUN = fit_laisk}.

  Because the Laisk fitting process involves two sets of linear fits, there are
  two possible graphs that can be created. When \code{plot_type} is
  \code{'first'}, this function will plot the individual A-Ci curves at each
  PPFD, along with the linear fits and the estimated intersection point. When
  \code{plot_type} is \code{'second'}, this function will plot the Laisk
  intercept vs. Laisk slope from the results of the first fits, along with a
  linear fit of Laisk intercept vs. Laisk slope. See
  \code{\link{fit_laisk}} for more details.

  Internally, this function uses \code{\link[lattice]{xyplot}} to perform the
  plotting. Optionally, additional arguments can be passed to \code{xyplot}.
  These should typically be limited to things like \code{xlim}, \code{ylim},
  \code{main}, and \code{grid}, since many other \code{xyplot} arguments will be
  set internally (such as \code{xlab}, \code{ylab}, \code{auto}, and others).

  See the help file for \code{\link{fit_laisk}} for an example using this
  function.
}


\value{
  A \code{trellis} object created by \code{lattice::xyplot}.
}

\examples{
# See the help file for fit_laisk for an example using this function.
}

\concept{exdf}
