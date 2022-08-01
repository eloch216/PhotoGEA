\name{choose_input_licor_files}

\alias{choose_input_licor_files}

\title{Choosing Licor Excel files}

\description{
  Opens up a window to choose Licor Excel files. Only available in interactive
  R sessions running in Microsoft Windows. The output from this function is
  typically passed to \code{\link{read_licor_file}} via \code{\link{lapply}}.
}

\usage{
  choose_input_licor_files()
}

\value{
  A character vector of full file paths.
}

\examples{
# Interactively select one or more Licor Excel files and read each one to create
# a list of exdf objects

\dontrun{

lapply(choose_input_licor_files(), function(fname) {
  read_licor_file(
    fname,
    preamble_data_rows = c(3, 5, 7, 9, 11, 13),
    variable_category_row = 14,
    variable_name_row = 15,
    variable_unit_row = 16,
    data_start_row = 17,
    timestamp_colname = 'time'
  )
})
}
}
