\name{optimizers}

\alias{optimizers}
\alias{optimizer_nmkb}
\alias{optimizer_deoptim}

\title{Optimizers}

\description{
  These functions return optimizers that meet requirements for the
  \code{OPTIM_FUN} input argument of \code{\link{fit_c3_aci}},
  \code{\link{fit_c3_variable_j}}, and \code{\link{fit_c4_aci}}. Essentially,
  they are wrappers for optimizers from other libraries that serve to
  standardize their inputs and outputs.
}

\usage{
  optimizer_nmkb(tol, maxfeval = 2000, restarts.max = 10)

  optimizer_deoptim(itermax, VTR = -Inf)
}

\arguments{
  \item{tol}{
    A convergence tolerance value; to be passed to \code{\link[dfoptim]{nmkb}}
    via its \code{control} input argument. A typical value is \code{1e-7}.
  }

  \item{maxfeval}{
    A maximum value for the number of function evaluations to allow during
    optimization; to be passed to \code{\link[dfoptim]{nmkb}} via its
    \code{control} input argument.
  }

  \item{restarts.max}{
    A maximum value for the number of restarts allowed during optimization;
    to be passed to \code{\link[dfoptim]{nmkb}} via its \code{control} input
    argument.
  }

  \item{itermax}{
    The maximum number of generations to be used; to be passed to
    \code{\link[DEoptim]{DEoptim}} via its \code{control} input argument. Note
    that when \code{VTR} is \code{-Inf}, the optimizer will always use the
    maximum number of generations. A typical value is \code{200}.
  }

  \item{VTR}{
    The value to be reached; to be passed to \code{\link[DEoptim]{DEoptim}} via
    its \code{control} input argument.
  }
}

\details{
  \code{optimizer_nmkb} is a wrapper for \code{\link[dfoptim]{nmkb}}.

  \code{optimizer_deoptim} is a wrapper for \code{\link[DEoptim]{DEoptim}}.
}
