\name{calculate_ball_berry_index}

\alias{calculate_ball_berry_index}

\title{Calculate the Ball-Berry index}

\description{
  Calculates the Ball-Berry index. This function can accomodate alternative
  column names for the variables taken from the Licor file in case they change
  at some point in the future. This function also checks the units of each
  required column and will produce an error if any units are incorrect.
}

\usage{
  calculate_ball_berry_index(
    licor_exdf,
    a_column_name = 'A',
    rhleaf_column_name = 'RHleaf',
    csurface_column_name = 'Csurface'
  )
}

\arguments{
  \item{licor_exdf}{
    An \code{exdf} object representing data from a Licor gas exchange
    measurement system.
  }

  \item{a_column_name}{
    The name of the column in \code{licor_exdf} that contains the net
    assimilation in \code{micromol m^(-2) s^(-1)}.
  }

  \item{rhleaf_column_name}{
    The name of the column in \code{licor_exdf} that contains the relative
    humidity at the leaf surface in \code{\%}.
  }

  \item{csurface_column_name}{
    The name of the column in \code{licor_exdf} that contains the CO2
    concentration at the leaf surface in \code{micromol mol^(-1)}.
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
# Read an example Licor file included in the PhotoGEA package, calculate the
# total pressure, calculate additional gas properties, and finally calculate the
# Ball-Berry index.
licor_file <- read_gasex_file(
  system.file('extdata', 'ball_berry_1.xlsx', package = 'PhotoGEA', mustWork = TRUE)
)

licor_file <- calculate_total_pressure(licor_file)

licor_file <- calculate_gas_properties(licor_file)

licor_file <- calculate_ball_berry_index(licor_file)

licor_file$units$bb_index      # View the units of the new `bb_index` column
licor_file$categories$bb_index # View the category of the new `bb_index` column
licor_file[,'bb_index']        # View the values of the new `bb_index` column
}

\concept{exdf}
