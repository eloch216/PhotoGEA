\name{residual_stats}

\alias{residual_stats}

\title{Calculate statistics that describe the residuals of a fit}

\description{
  Calculates several key statistics from the residuals of of a fit: the residual
  sum of squares (\code{RSS}), the mean squared error (\code{MSE}), the root
  mean squared error (\code{RMSE}), and the residual standard error
  (\code{RSE}).
}

\usage{
  residual_stats(fit_residuals, units, nparam)
}

\arguments{
  \item{fit_residuals}{
    A numeric vector representing the residuals from a fit, i.e., the
    differences between the measured and fitted values.
  }

  \item{units}{
    A string expressing the units of the residuals.
  }

  \item{nparam}{
    The number of free parameters that were varied when performing the fit.
  }
}

\details{
  This function calculates several model-independent measures of the quality of
  a fit. The basis for these statistics are the \code{residuals} (also known as
  the \code{errors}). If the measured values of a quantity \code{y} are given by
  \code{y_measured} and the fitted values are \code{y_fitted}, then the
  residuals are defined to be \code{residual = y_measured - y_fitted}. The key
  statistics that can be calculated from the residuals are as follows:
  \itemize{
    \item The residual sum of squares (\code{RSS}) is also known as the sum of
          squared errors (\code{SSE}). As its name implies, it is simply the sum
          of all the squared residuals: \code{RSS = sum(residuals^2)}.
    \item The mean squared error (\code{MSE}) is the mean value of the squared
          residuals: \code{MSE = sum(residuals^2) / n = RSS / n}, where \code{n}
          is the number of residuals.
    \item The root mean squared error (\code{RMSE}) is the square root of the
          mean squared error: \code{RMSE = sqrt(MSE) = sqrt(RSS / n)}.
    \item The residual standard error \code{RSE} is given by \code{RSE =
          sqrt(RSS / dof)}, where \code{dof = n - nparam} is the number of
          degrees of freedom involved in the fit.
  }

  For a given model, the \code{RMSE} is usually a good way to compare the
  quality of different fits. When trying to decide which model best fits the
  measured data, the \code{RSE} may be a more appropriate metric since it
  controls for the number of parameters in the model.
}

\value{
  An \code{exdf} object with one row and the following columns: \code{npts} (the
  number of residual values), \code{nparam}, \code{dof}, \code{RSS}, \code{MSE},
  \code{RMSE}, and \code{RSE}.
}

\examples{
# Generate some random residuals
residuals <- runif(10, -1, 1)

# Calculate residual stats as if these values had units of `kg` and were related
# to a model with 3 free parameters
residual_stats(residuals, 'kg', 3)
}

\concept{exdf}
