\name{identify_c3_limiting_processes}

\alias{identify_c3_limiting_processes}

\title{Identify C3 Limiting Processes}

\description{
  Identify limiting processes in a C3 curve, typically the result of a fit.
  It is rate for users to call this function directly because it is used
  internally by \code{\link{fit_c3_aci}} and \code{\link{fit_c3_variable_j}}.
}

\usage{
  identify_c3_limiting_processes(
    data_table,
    a_column_name = 'A_fit',
    ac_column_name = 'Ac',
    aj_column_name = 'Aj',
    ap_column_name = 'Ap',
    tol = 1e-3
)
}

\arguments{
  \item{data_table}{
    A table-like R object such as a data frame or an \code{exdf}.
  }

  \item{a_column_name}{
    The name of the column in \code{data_table} that contains the modeled net CO2
    assimilation rate in \code{micromol m^(-2) s^(-1)}.
  }

  \item{ac_column_name}{
    The name of the column in \code{data_table} that contains the modeled
    Rubisco-limited net CO2 assimilation rate in \code{micromol m^(-2) s^(-1)}.
  }

  \item{aj_column_name}{
    The name of the column in \code{data_table} that contains the modeled
    RuBP-regeneration-limited net CO2 assimilation rate in
    \code{micromol m^(-2) s^(-1)}.
  }

  \item{ap_column_name}{
    The name of the column in \code{data_table} that contains the modeled
    TPU-limited net CO2 assimilation rate in \code{micromol m^(-2) s^(-1)}.
  }

  \item{tol}{
    A relative tolerance factor used to identify when two rates are equal.
  }
}

\details{
  For a C3 leaf, \code{An} is given by either \code{Ac}, \code{Aj}, or
  \code{Ap}. See the documentation for \code{\link{calculate_c3_assimilation}}
  for more information.

  This function first identifies points where \code{An = Ac}, \code{An = Aj},
  and \code{An = Ap}. The results are stored in columns called
  \code{Ac_limiting}, \code{Aj_limiting}, and \code{Ap_limiting}, where a value
  of \code{TRUE} indicates that the corresponding process is limiting.

  Then, the overall limiting state is specified in the \code{limiting_process}
  column. For example, points where \code{An} equals \code{Ac} but not \code{Aj}
  or \code{Ap} are designated by \code{limiting_process = 'Ac'}, and likewise
  for the other potential limiting processes. If more than one process is
  limiting for a point, \code{limiting_process} is set to \code{'co-limited'}.

}

\value{
  An \code{exdf} object based on \code{licor_exdf} that includes new columns as
  described above: \code{Ac_limiting}, \code{Aj_limiting}, \code{Ap_limiting},
  and \code{limiting_process}. The categories of these new columns are
  set to \code{identify_c3_limiting_processes} to indicate that they were
  created using this function.
}

\examples{
# Identify limiting processes in an example curve
example_curve <- exdf(
  data.frame(
    A_fit = c(1.0, 2.0, 3.0, 4.0, 4.0),
    Ac    = c(1.0, 2.0, 5.0, 8.0, 9.0),
    Aj    = c(2.0, 2.5, 3.0, 4.0, 8.0),
    Ap    = c(NA,  NA,  4.0, 4.0, 4.0)
  ),
  units = data.frame(
    A_fit = 'micromol m^(-2) s^(-1)',
    Ac    = 'micromol m^(-2) s^(-1)',
    Aj    = 'micromol m^(-2) s^(-1)',
    Ap    = 'micromol m^(-2) s^(-1)',
    stringsAsFactors = FALSE
  )
)

identify_c3_limiting_processes(example_curve)

# This function also works for data frames
identify_c3_limiting_processes(example_curve$main_data)
}

\concept{exdf}
