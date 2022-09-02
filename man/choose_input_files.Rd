\name{choose_input_files}

\alias{choose_input_files}
\alias{choose_input_licor_files}
\alias{choose_input_tdl_files}

\title{Choosing input files}

\description{
  Tools for choosing input files via dialog windows.
}

\usage{
  choose_input_files()

  choose_input_licor_files()

  choose_input_tdl_files()
}

\details{
  These functions are only available in interactive sessions; moreover,
  \code{choose_input_licor_files} and \code{choose_input_tdl_files} are only
  available in Microsoft Windows.

  \code{choose_input_files} will prompt the user to select a single file, and
  will return full file paths for all files in the same directory that have the
  same extension. The output from this function is typically passed to
  \code{\link{read_licor_file}} or \code{\link{read_tdl_file}}.

  \code{choose_input_licor_files} can be used to select one or more Microsoft
  Excel files (with extension \code{*.xlsx}). The output from this function is
  typically passed to \code{\link{read_licor_file}} via \code{\link{lapply}}.

  \code{choose_input_tdl_files} can be used to select one or more TDL data files
  (with extension \code{*.dat}). The output from this function is typically
  passed to \code{\link{read_tdl_file}} via \code{\link{lapply}}.
}

\value{
  A character vector of full file paths.
}

\examples{
# Interactively select a single file and get full file paths to all
# other files in the same directory that have the same extension

\dontrun{

file_paths <- choose_input_files()
}


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