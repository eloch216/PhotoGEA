\name{optimizers}

\alias{optimizers}
\alias{optimizer_nmkb}
\alias{optimizer_hjkb}
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

  optimizer_hjkb(tol, maxfeval = Inf, target = Inf)

  optimizer_deoptim(itermax, VTR = -Inf)
}

\arguments{
  \item{tol}{
    A convergence tolerance value; to be passed to \code{\link[dfoptim]{nmkb}}
    or \code{\link[dfoptim]{hjkb}} via their \code{control} input arguments. A
    typical value is \code{1e-7}.
  }

  \item{maxfeval}{
    A maximum value for the number of function evaluations to allow during
    optimization; to be passed to \code{\link[dfoptim]{nmkb}} or
    \code{\link[dfoptim]{hjkb}} via their \code{control} input arguments.
  }

  \item{target}{
    A real number restricting the absolute function value; to be passed to
    \code{\link[dfoptim]{hjkb}} via its \code{control} input argument.
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

  \code{optimizer_hjkb} is a wrapper for \code{\link[dfoptim]{hjkb}}.

  \code{optimizer_deoptim} is a wrapper for \code{\link[DEoptim]{DEoptim}}.

  See the documentation for those functions for more information about how the
  optimizers work.
}
