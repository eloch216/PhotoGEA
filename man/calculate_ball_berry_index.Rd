\name{calculate_ball_berry_index}

\alias{calculate_ball_berry_index}

\title{calculate_ball_berry_index}

\description{
  Calculates the Ball-Berry index. This function can accomodate alternative
  column names for the variables taken from the Licor file in case they change
  at some point in the future. This function also checks the units of each
  required column and will produce an error if any units are incorrect.
}

\usage{
  calculate_ball_berry_index(
    licor_exdf,
    a_column_name,
    rhleaf_column_name,
    csurface_column_name
  )
}

\arguments{
  \item{licor_exdf}{
    An \code{exdf} object representing data from a Licor gas exchange
    measurement system.
  }

  \item{a_column_name}{
    The name of the column in \code{licor_exdf} that contains the net
    assimilation (typically \code{'A'}).
  }

  \item{rhleaf_column_name}{
    The name of the column in \code{licor_exdf} that contains the relative
    humidity at the leaf surface (typically \code{'RHleaf'}).
  }

  \item{csurface_column_name}{
    The name of the column in \code{licor_exdf} that contains the CO2
    concentration at the leaf surface (typically \code{'Csurface'}).
  }
}

\details{
  The Ball-Berry index is defined as \code{A * h_s / c_s}, where
  \code{A} is the net assimilation rate, \code{h_s} is the relative humidity
  at the leaf surface, and \code{c_s} is the CO2 concentration at the leaf
  surface. This variable is a key part of the Ball-Berry model, which assumes
  that stomatal conductance is linearly related to the Ball-Berry index. For
  more information, please see the original publication describing the model:
  Ball, J. T., Woodrow, I. E. and Berry, J. A. "A Model Predicting Stomatal
  Conductance and its Contribution to the Control of Photosynthesis under
  Different Environmental Conditions." in "Progress in Photosynthesis Research:
  Volume 4" (1986) [\doi{10.1007/978-94-017-0519-6_48}].

  Typically, the relative humidity and CO2 concentration at the leaf surface are
  not included in Licor output files. Instead, the output files only include the
  relative humidity and CO2 concentration in the sample chamber, and conditions
  at the leaf surface may be slightly different. These required inputs can be
  calculated using the \code{\link{calculate_gas_properties}} function.
}

\value{
  An \code{exdf} object based on \code{licor_exdf} that includes the Ball-Berry
  index as a new column called \code{bb_index}. The category of this new column
  is \code{calculate_ball_berry_index} to indicate that it was created using
  this function.
}

\examples{
# Read an example Licor file included in the PhotoGEA package, calculate
# additional gas properties, and then calculate the Ball-Berry index.
licor_file <- read_licor_file(
  system.file("extdata", "ball_berry_1.xlsx", package = "PhotoGEA"),
  preamble_data_rows = c(3, 5, 7, 9, 11, 13),
  variable_category_row = 14,
  variable_name_row = 15,
  variable_unit_row = 16,
  data_start_row = 17,
  timestamp_colname = 'time'
)

licor_file <- calculate_gas_properties(
  licor_file,
  'A',
  'Ca',
  'DeltaPcham',
  'E',
  'gbw',
  'gsw',
  'H2O_s',
  'Pa',
  'TleafCnd'
)

licor_file <- calculate_ball_berry_index(
    licor_file,
    'A',
    'RHleaf',
    'Csurface'
)

licor_file$units$bb_index      # View the units of the new `bb_index` column
licor_file$categories$bb_index # View the category of the new `bb_index` column
licor_file[,'bb_index']        # View the values of the new `bb_index` column
}

\concept{exdf}