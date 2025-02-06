\name{initial_guess_c4_aci_hyperbola}

\alias{initial_guess_c4_aci_hyperbola}

\title{Make an initial guess of C4 hyperbola parameter values for one curve}

\description{
  Creates a function that makes an initial guess of C4 hyperbola model
  parameter values for one curve. This function is used internally by
  \code{\link{fit_c4_aci_hyperbola}}.

  Values estimated by this guessing function should be considered inaccurate,
  and should always be improved upon by an optimizer.
}

\usage{
  initial_guess_c4_aci_hyperbola(
    a_column_name = 'A'
  )
}

\arguments{
  \item{a_column_name}{
    The name of the column in \code{rc_exdf} that contains the net assimilation
    rate in \code{micromol m^(-2) s^(-1)}.
  }
}

\details{
  Here we estimate values of \code{c4_curvature}, \code{c4_slope}, \code{rL},
  and \code{Vmax} from a measured C4 CO2 response curve. For more information
  about these parameters, see the documentation for
  \code{\link{calculate_c4_assimilation_hyperbola}}.

  Here we take a very simple approach to forming the initial guess. We always
  choose \code{c4_curvature = 0.5}, \code{c4_slope = 1.0}, and \code{rL = 0.0}.
  For \code{Vmax}, we use \code{Vmax = max{A} - rL_guess}, where \code{max{A}}
  is the largest observed net CO2 assimilation rate and \code{rL_guess} is the
  guess for \code{rL}.
}

\value{
  A function with one input argument \code{rc_exdf}, which should be an
  \code{exdf} object representing one C4 CO2 response curve. The return value of
  this function will be a numeric vector with four elements, representing the
  values of \code{c4_curvature}, \code{c4_slope}, \code{rL}, and \code{Vmax}
  (in that order).
}

\examples{
# Read an example Licor file included in the PhotoGEA package
licor_file <- read_gasex_file(
  PhotoGEA_example_file_path('c4_aci_1.xlsx')
)

# Define a new column that uniquely identifies each curve
licor_file[, 'species_plot'] <-
  paste(licor_file[, 'species'], '-', licor_file[, 'plot'] )

# Organize the data
licor_file <- organize_response_curve_data(
    licor_file,
    'species_plot',
    c(9, 10, 16),
    'CO2_r_sp'
)

# Create the guessing function
guessing_func <- initial_guess_c4_aci_hyperbola()

# Apply it and see the initial guesses for each curve
str(by(licor_file, licor_file[, 'species_plot'], guessing_func))
}

\concept{exdf}
