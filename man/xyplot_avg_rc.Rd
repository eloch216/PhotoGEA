\name{xyplot_avg_rc}

\alias{xyplot_avg_rc}

\title{xyplot_avg_rc}

\description{
  A wrapper for \code{lattice::xyplot} that plots average response curves with
  error bars.
}

\usage{
  xyplot_avg_rc(
    Y,
    X,
    point_identifier,
    group_identifier,
    x_error_bars = FALSE,
    cols = multi_curve_colors(),
    eb_length = 0.05,
    eb_lwd = 1,
    ...
  )
}

\arguments{
  \item{Y}{A numeric vector of y-values.}

  \item{X}{A numeric vector of x-values with the same length as \code{Y}}

  \item{point_identifier}{
    A vector with the same length as \code{Y} that indicates the location of
    each (x, y) pair along the response curve; typically this is the
    \code{seq_num} column of an \code{exdf} object.
  }

  \item{group_identifier}{
    A vector with the same length as \code{Y} that indicates the "group" of each
    response curve.
  }

  \item{x_error_bars}{
    A logical value indicating whether to plot x-axis error bars.
  }

  \item{cols}{A vector of color specifications.}

  \item{eb_length}{The width of the error bars.}

  \item{eb_lwd}{The line width (thickness) of the error bars.}

  \item{\dots}{Additional arguments to be passed to \code{lattice::xyplot}.}
}

\details{
  This function calculates average values of \code{X} and \code{Y} at each value
  of the \code{point_identifier} across groups defined by
  \code{group_identifier}, and then uses these values to plot average response
  curves for each group. Error bars are determined by calculating the standard
  errors of \code{X} and \code{Y} at each value of the \code{point_identifier}
  across groups defined by \code{group_identifier}.
}

\value{
  A \code{trellis} object created by \code{lattice::xyplot}.
}

\examples{
# Read an example Licor file included in the PhotoGEA package
licor_file <- read_licor_file(
  system.file("extdata", "ball_berry_1.xlsx", package = "PhotoGEA"),
  c(3, 5, 7, 9, 11, 13), 14, 15, 16, 17, 'time'
)

# Organize the response curve data
licor_file <- organize_response_curve_data(
  licor_file,
  c('species', 'plot'),
  c(),
  'Qin'
)

# Plot the average light response curve for each species (here there is only one
# curve for tobacco, so there are no tobacco error bars)
xyplot_avg_rc(
  licor_file[, 'A'],
  licor_file[, 'Qin'],
  licor_file[, 'seq_num'],
  licor_file[, 'species'],
  ylim = c(0, 50),
  xlab = paste0('Incident PPFD (', licor_file$units$Qin, ')'),
  ylab = paste0('Average net assimilation (', licor_file$units$A, ')'),
  auto = TRUE,
  grid = TRUE
)
}
