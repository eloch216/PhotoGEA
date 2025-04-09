\name{deprecated}

\alias{deprecated}
\alias{read_tdl_file}
\alias{read_licor_file}
\alias{check_licor_data}
\alias{calculate_arrhenius}
\alias{calculate_peaked_gaussian}

\title{Deprecated functions}

\description{
  Deprecated functions that will be fully removed in future releases. Each of
  these functions will produce an error when called that will redirect the user
  to a suitable replacement.
}

\usage{
  read_tdl_file(\dots)

  read_licor_file(\dots)

  check_licor_data(\dots)

  calculate_arrhenius(\dots)

  calculate_peaked_gaussian(\dots)
}

\arguments{
  \item{\dots}{
    Additional arguments; currently unused.
  }
}

\value{None of the deprecated functions return anything.}

\examples{
# These functions all throw errors, so we will wrap them in `tryCatch` here

tryCatch(
  read_tdl_file(),
  error = function(e) {print(e)}
)

tryCatch(
  read_licor_file(),
  error = function(e) {print(e)}
)

tryCatch(
  check_licor_data(),
  error = function(e) {print(e)}
)

tryCatch(
  calculate_arrhenius(),
  error = function(e) {print(e)}
)

tryCatch(
  calculate_peaked_gaussian(),
  error = function(e) {print(e)}
)
}
