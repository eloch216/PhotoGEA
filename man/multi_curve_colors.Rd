\name{multi_curve_colors}

\alias{multi_curve_colors}

\title{Set of colors for plotting multiple curves}

\description{
  A set of color specifications that work reasonably well for plotting multiple
  curves on the same axes.
}

\usage{multi_curve_colors()}

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
