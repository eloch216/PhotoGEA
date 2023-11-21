\name{multi_curve_colors}

\alias{multi_curve_colors}
\alias{multi_curve_line_colors}
\alias{multi_curve_point_colors}

\title{Set of colors for plotting multiple curves}

\description{
  \code{multi_curve_colors} returns a vector of color specifications that work
  reasonably well for plotting multiple curves on the same axes.

  \code{multi_curve_line_colors} returns the same vector, but with the
  first color set to be transparent. \code{multi_curve_point_colors} also
  returns the same vector, but with all colors except the first set to
  transparent. These color specifications can be helpful when plotting measured
  data along with fits, allowing the data to be displayed as points and the fits
  as lines.
}

\usage{
  multi_curve_colors()

  multi_curve_line_colors()

  multi_curve_point_colors()
}

\details{
  The color set was originally formed by calling the following:

  \code{multi_curve_colors <- c(
      "#000000",
      RColorBrewer::brewer.pal(8, "Set2"),
      RColorBrewer::brewer.pal(12, "Paired")[c(1:10,12)],
      RColorBrewer::brewer.pal(8, "Dark2")
  )}
}

\value{
  A character vector with 28 elements, each of which is a hexadecimal color
  specification.
}
