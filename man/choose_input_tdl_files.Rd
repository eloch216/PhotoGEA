\name{choose_input_tdl_files}

\alias{choose_input_tdl_files}

\title{Choosing TDL data files}

\description{
  Opens up a window to choose TDL data files. Only available in interactive
  R sessions running in Microsoft Windows. The output from this function is
  typically passed to \code{\link{read_tdl_file}} via \code{\link{lapply}}.
}

\usage{
  choose_input_tdl_files()
}

\value{
  A character vector of full file paths.
}

\examples{
# Interactively select one or more TDL data files and read each one to create a
# list of exdf objects

\dontrun{

lapply(choose_input_tdl_files(), function(fname) {
  read_tdl_file(
    fname,
    rows_to_skip = 1,
    variable_name_row = 2,
    variable_unit_row = 3,
    data_start_row = 5,
    timestamp_colname = 'TIMESTAMP'
  )
})
}
}
