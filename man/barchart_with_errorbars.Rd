\name{barchart_with_errorbars}

\alias{barchart_with_errorbars}
\alias{bwplot_wrapper}

\title{Barcharts with error bars}

\description{
  \code{barchart_with_errorbars} is a wrapper for \code{lattice::barchart} that
  includes error bars on the chart, while \code{bwplot_wrapper} is a simple
  wrapper for \code{lattice::bwplot} that gives it the same function signature
  as \code{barchart_with_errorbars}.
}

\usage{
  barchart_with_errorbars(
    Y,
    X,
    eb_width = 0.2,
    eb_lwd = 1,
    eb_col = 'black',
    na.rm = TRUE,
    remove_outliers = FALSE,
    ...
  )

  bwplot_wrapper(Y, X, ...)
}

\arguments{
  \item{Y}{A numeric vector.}

  \item{X}{
    A vector with the same length as \code{Y} that can be used as a factor to
    split \code{Y} into one or more distinct subsets.
  }

  \item{eb_width}{The width of the error bars.}

  \item{eb_lwd}{The line width (thickness) of the error bars.}

  \item{eb_col}{The color of the error bars.}

  \item{na.rm}{
    A logical value indicating whether or not to remove NA values before
    calculating means and standard errors.
  }

  \item{remove_outliers}{
    A logical value indicating whether or not to remove outliers using
    \code{\link{exclude_outliers}} before calculating means and standard errors.
  }

  \item{\dots}{
    Additional arguments to be passed to \code{lattice::barchart} or
    \code{lattice::bwplot}.
  }
}

\details{
  The \code{barchart_with_errorbars} function uses \code{\link{tapply}} to
  calculate the mean and standard error for each subset of \code{Y} as
  determined by the values of \code{X}. In other words, \code{means <-
  tapply(Y, X, mean)}, and similar for the standard errors. The mean values are
  represented as bars in the final plot, while the standard error is used to
  create error bars located at \code{mean +/- standard_error}.

  The \code{bwplot_wrapper} function is a simple wrapper for
  \code{lattice::bwplot} that gives it the same input arguments as
  \code{barchart_with_errorbars}. In other words, the same \code{X} and \code{Y}
  vectors can be used to create a barchart using \code{barchart_with_errorbars}
  or a box-whisker plot with \code{bwplot_wrapper}.
}

\value{
  A \code{trellis} object created by \code{lattice::barchart} or
  \code{lattice::bwplot}.
}

\examples{
# Read an example Licor file included in the PhotoGEA package
licor_file <- read_gasex_file(
  PhotoGEA_example_file_path('ball_berry_1.xlsx')
)

# Plot the average assimilation value for each species. (Note: this is not a
# meaningful calculation since we are combining assimilation values measured
# at different PPFD.)
barchart_with_errorbars(
  licor_file[, 'A'],
  licor_file[, 'species'],
  ylim = c(0, 50),
  xlab = 'Species',
  ylab = paste0('Net assimilation (', licor_file$units$A, ')')
)

# Make a box-whisker plot using the same data. (Note: this is not a meaningful
# plot since we are combining assimilation values measured at different PPFD.)
bwplot_wrapper(
  licor_file[, 'A'],
  licor_file[, 'species'],
  ylim = c(0, 50),
  xlab = 'Species',
  ylab = paste0('Net assimilation (', licor_file$units$A, ')')
)

# Another way to create the plots. This method illustrates the utility of the
# bwplot_wrapper function.
plot_parameters <- list(
  Y = licor_file[, 'A'],
  X = licor_file[, 'species'],
  ylim = c(0, 50),
  xlab = 'Species',
  ylab = paste0('Net assimilation (', licor_file$units$A, ')')
)
do.call(barchart_with_errorbars, plot_parameters)
do.call(bwplot_wrapper, plot_parameters)
}
