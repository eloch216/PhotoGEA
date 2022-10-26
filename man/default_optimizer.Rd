\name{default_optimizer}

\alias{default_optimizer}

\title{A default optimizer}

\description{
  This function returns an optimizer that meets the requirements for the
  \code{OPTIM_FUN} input argument of \code{\link{fit_c3_aci}} and other fitting
  functions. Essentially, this is a wrapper for \code{dfoptim::nmkb}.
}

\usage{default_optimizer(tol = 1e-7, maxfeval = 2000, restarts.max = 10)}

\arguments{
  \item{tol}{
    A convergence tolerance value; to be passed to \code{dfoptim::nmkb} via its
    \code{control} input argument.
  }

  \item{maxfeval}{
    A maximum value for the number of function evaluations to allow during
    optimization; to be passed to \code{dfoptim::nmkb} via its \code{control}
    input argument.
  }

  \item{restarts.max}{
    A maximum value for the number of restarts allowed during optimization;
    to be passed to \code{dfoptim::nmkb} via its \code{control} input argument.
  }
}
