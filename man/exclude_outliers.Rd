\name{exclude_outliers}

\alias{exclude_outliers}
\alias{exclude_outliers.data.frame}
\alias{exclude_outliers.exdf}

\title{Exclude outliers from a data set}

\description{
  Excludes outliers from a data set using the "1.5 interquartile range" rule.
}

\usage{
  exclude_outliers(x, col_for_analysis, INDICES)

  \method{exclude_outliers}{data.frame}(x, col_for_analysis, INDICES)

  \method{exclude_outliers}{exdf}(x, col_for_analysis, INDICES)
}

\arguments{
  \item{x}{A data table}

  \item{col_for_analysis}{
    The name of a column of \code{x} that should be used to determine outliers.
  }

  \item{INDICES}{
    A factor or list of factors that each \code{nrow(x)} elements.
  }
}

\value{
  A copy of \code{X} where outliers have been removed.
}

\details{
  \code{exclude_outliers} is generic, with methods defined for data frames and
  \code{exdf} objects. This function uses a simple rule to detect outliers,
  where any point that deviates from the mean by more than \code{1.5 * IQR},
  where \code{IQR} is the interquartile range, is said to be an outlier. This
  method is also sometimes referred to as "Tukey's Fences," as seen in the
  \href{https://en.wikipedia.org/wiki/Outlier}{Wikipedia page about outliers}.

  For data sets with extreme outliers, it may be necessary to exclude outliers
  more than once to actually remove them all.
}

\seealso{\code{\link{exdf}}}

\examples{
# Read a Licor file included with the PhotoGEA package; this file includes
# several light response curves that can be identified by the 'species' and
# 'plot' columns.
licor_file <- read_gasex_file(
  PhotoGEA_example_file_path('ball_berry_1.xlsx')
)

# Exclude points from each response curve in the data where the leaf temperature
# is determined to be an outlier
licor_file_clean <- exclude_outliers(
  licor_file,
  'TleafCnd',
  list(licor_file[, 'species'], licor_file[, 'plot'])
)

# Check to see how many points remain after removing outliers
str(list('original' = nrow(licor_file), 'clean' = nrow(licor_file_clean)))
}

\concept{exdf}
