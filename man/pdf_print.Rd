\name{pdf_print}

\alias{pdf_print}

\title{Print a plot object or save it to a PDF}

\description{
  A convenience function that either displays a plot object in an R graphics
  window or saves it to a PDF.
}

\usage{
  pdf_print(
    plot_obj,
    width = 7,
    height = 7,
    save_to_pdf = FALSE,
    file = NULL,
    new_window = TRUE,
    \dots
  )
}

\arguments{
  \item{plot_obj}{
    A plotting object that can be printed, such as a trellis object returned by
    a call to \code{\link[lattice]{xyplot}}.
  }

  \item{width}{The width of the figure in inches.}

  \item{height}{The width of the figure in inches.}

  \item{save_to_pdf}{
    When \code{save_to_pdf} is \code{TRUE}, \code{plot_obj} will be saved as a
    PDF; otherwise it will be printed to an R graphics window.
  }

  \item{file}{
    A file name to use when \code{save_to_pdf} is \code{TRUE}. If \code{file} is
    \code{NULL}, then the default value will be determined by the
    \code{\link{pdf}} function.
  }

  \item{new_window}{
    When printing \code{plot_obj} to an R graphics window, a new window will be
    created if \code{new_window} is \code{TRUE}. Otherwise, the plot will
    replace the currently active plot window (if one exists).
  }

  \item{\dots}{Additional arguments to be passed to \code{\link{pdf}}.}
}

\details{
  This function is helpful when developing and using analysis scripts. In this
  context, it is recommended to define a boolean called \code{SAVE_TO_PDF} early
  in the script and to always use \code{pdf_print} when creating figures,
  passing the boolean as the \code{save_to_pdf} input argument. Figures can be
  initially displayed in R (setting \code{SAVE_TO_PDF = FALSE}), and then saved
  as PDFs once graphing parameters have been optimized (setting
  \code{SAVE_TO_PDF = TRUE}).

  Note that calling \code{\link{pdf}} from the command line (as is done
  internally by \code{pdf_print}) is different than exporting an R graphics
  object as a PDF from RGui or RStudio. For some reason, RGui and RStudio
  override some of the \code{\link{pdf}} defaults and set \code{useDingbats} to
  \code{TRUE}. This setting almost always causes problems when opening the PDFs
  in software like Adobe Illustrator or Inkscape.
}

\value{The \code{pdf_print} function does not return anything.}

\examples{
SAVE_TO_PDF = FALSE # change this to TRUE to save to a PDF

pdf_print(
  lattice::xyplot(
    1:4 ~ 11:14,
    xlab = 'X',
    ylab = 'Y',
    type = 'b'
  ),
  save_to_pdf = SAVE_TO_PDF,
  file = 'example.pdf', # this name will only be used when saving to a PDF
  new_window = FALSE    # necessary for rendering the documentation examples
)
}
